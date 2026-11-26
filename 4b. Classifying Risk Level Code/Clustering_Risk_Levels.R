# =================================================================
# K-means clustering into 3 risk levels (Low / Medium / High)
# Fasting glucose and HbA1c are used on their original scale (not logged), so
# the diabetes High group is mostly people with clearly diabetic values.
# Hypertension, dyslipidaemia, diabetes - KNHANES 2024, adults aged 40-74
# =================================================================
library(dplyr)
library(data.table)
library(mice)
library(ggplot2)
library(tidyr)

# Show graphs only - never write them to a file (stops Rplots.pdf when run with Rscript)
if (!interactive()) pdf(NULL)

# =================================================================
# STEP 1. Read dt and keep ages 40 to 74
# =================================================================
dt <- fread("hn24_all.csv")
dt <- dt[age >= 40 & age <= 74, ]
cat("Step 1: rows aged 40-74 =", nrow(dt), "\n")

# =================================================================
# STEP 2. Select the required factors for each disease
# =================================================================
# English name = "original KNHANES name"
all_vars <- c(
  # Person ID (used only to merge the three datasets in Step 6)
  ID = "ID",
  # Demographics
  sex = "sex", age = "age", income_quartile_indiv = "incm", education = "edu",
  # Body measures & blood pressure
  bmi = "HE_BMI", waist_cm = "HE_wc", sbp = "HE_sbp", dbp = "HE_dbp",
  # Blood tests
  fasting_glucose = "HE_glu", hba1c = "HE_HbA1c", total_cholesterol = "HE_chol",
  hdl = "HE_HDL_st2", triglycerides = "HE_TG", ldl = "HE_LDL_drct",
  # Lifestyle
  drink_freq_yr = "BD1_11", cigs_per_day = "BS3_2", strength_days_week = "BE5_1",
  sodium_mg = "N_NA",
  # Helper columns: used ONLY to impute other variables in Step 6, never clustered
  height_cm = "HE_ht", weight_kg = "HE_wt", wbc = "HE_WBC", haemoglobin = "HE_HB",
  uric_acid = "HE_Uacid", ast = "HE_ast", alt = "HE_alt", mcv = "HE_MCV",
  urine_sodium = "HE_UNa", urine_creatinine = "HE_Ucrea",
  # Physical activity components (days, hours, minutes)
  vig_work_days = "BE3_72", vig_work_h = "BE3_73", vig_work_m = "BE3_74",
  vig_leis_days = "BE3_76", vig_leis_h = "BE3_77", vig_leis_m = "BE3_78",
  mod_work_days = "BE3_82", mod_work_h = "BE3_83", mod_work_m = "BE3_84",
  mod_leis_days = "BE3_86", mod_leis_h = "BE3_87", mod_leis_m = "BE3_88",
  sit_h = "BE8_1", sit_m = "BE8_2",
  # Sleep components (bedtime and wake-up time, weekday and weekend)
  bed_wd_h = "BP16_11", bed_wd_m = "BP16_12", wake_wd_h = "BP16_13", wake_wd_m = "BP16_14",
  bed_we_h = "BP16_21", bed_we_m = "BP16_22", wake_we_h = "BP16_23", wake_we_m = "BP16_24",
  # Family history (father, mother, siblings)
  fh_1 = NA, fh_2 = NA, fh_3 = NA,   # filled per disease below
  # Doctor diagnosis (filled per disease below)
  diagnosed = NA
)

# Helper columns for imputation only (see Step 6)
helper_cols <- c("height_cm", "weight_kg", "wbc", "haemoglobin", "uric_acid",
                 "ast", "alt", "mcv", "urine_sodium", "urine_creatinine")

# Columns every disease uses
shared <- c("ID", "sex", "age", "income_quartile_indiv", "education",
            "bmi", "waist_cm", "sbp", "dbp", "strength_days_week", helper_cols,
            grep("^(vig_|mod_|sit_|bed_|wake_)", names(all_vars), value = TRUE),
            "fh_1", "fh_2", "fh_3", "diagnosed")

# Extra columns for each disease
disease_setup <- list(
  hypertension = list(
    extra     = c("fasting_glucose", "total_cholesterol", "hdl", "triglycerides", "ldl",
                  "drink_freq_yr", "cigs_per_day", "sodium_mg"),
    family    = c("HE_HPfh1", "HE_HPfh2", "HE_HPfh3"),
    diagnosis = "DI1_dg"),
  dyslipidaemia = list(
    extra     = c("total_cholesterol", "hdl", "triglycerides", "ldl",
                  "cigs_per_day", "drink_freq_yr"),
    family    = c("HE_HLfh1", "HE_HLfh2", "HE_HLfh3"),
    diagnosis = "DI2_dg"),
  diabetes = list(
    extra     = c("fasting_glucose", "hba1c", "hdl", "triglycerides", "cigs_per_day"),
    family    = c("HE_DMfh1", "HE_DMfh2", "HE_DMfh3"),
    diagnosis = "DE1_dg")
)

dts <- list()
for (d in names(disease_setup)) {
  s    <- disease_setup[[d]]
  cols <- all_vars[c(shared, s$extra)]
  cols[c("fh_1", "fh_2", "fh_3")] <- s$family
  cols["diagnosed"]               <- s$diagnosis
  dts[[d]] <- dt %>% select(all_of(cols))   # selects AND renames
  cat("Step 2:", d, "->", ncol(dts[[d]]), "columns\n")
}


# =================================================================
# STEP 3. Add vigorous and moderate activity time (+ sitting and sleep time)
# =================================================================
# Activity hours per week = days x (hours + minutes / 60), work + leisure
#   days: 8 = not active -> 0, 9 = don't know -> NA
#   hours / minutes: 88 = not active -> 0, 99 = don't know -> NA
# Sitting hours per day = hours + minutes / 60   (99 -> NA)
# Sleep hours = wake-up time - bedtime, weighted 5 weekdays + 2 weekend days

pa_hours_week <- function(days, h, m) {
  days <- case_when(days == 8 ~ 0, days == 9 ~ NA_real_, TRUE ~ as.numeric(days))
  h    <- case_when(h == 88 ~ 0, h == 99 ~ NA_real_, TRUE ~ as.numeric(h))
  m    <- case_when(m == 88 ~ 0, m == 99 ~ NA_real_, TRUE ~ as.numeric(m))
  ifelse(days == 0, 0, days * (h + m / 60))
}

sleep_duration <- function(bed_h, bed_m, wake_h, wake_m) {
  bad <- bed_h == 99 | bed_m == 99 | wake_h == 99 | wake_m == 99
  dur <- ((wake_h + wake_m / 60) - (bed_h + bed_m / 60)) %% 24
  ifelse(bad | dur < 2 | dur > 16, NA, dur)     # implausible -> NA
}

for (d in names(dts)) {
  dts[[d]] <- dts[[d]] %>%
    mutate(
      vig_pa_hrs_wk   = pa_hours_week(vig_work_days, vig_work_h, vig_work_m) +
        pa_hours_week(vig_leis_days, vig_leis_h, vig_leis_m),
      mod_pa_hrs_wk   = pa_hours_week(mod_work_days, mod_work_h, mod_work_m) +
        pa_hours_week(mod_leis_days, mod_leis_h, mod_leis_m),
      sitting_hrs_day = ifelse(sit_h == 99 | sit_m == 99, NA, sit_h + sit_m / 60),
      sleep_hrs       = (5 * sleep_duration(bed_wd_h, bed_wd_m, wake_wd_h, wake_wd_m) +
                           2 * sleep_duration(bed_we_h, bed_we_m, wake_we_h, wake_we_m)) / 7
    ) %>%
    select(-matches("^(vig_work|vig_leis|mod_work|mod_leis|sit_|bed_|wake_)"))
}
cat("Step 3: activity, sitting and sleep time added\n")


# =================================================================
# STEP 4. Clean special codes (8, 9, 88, 99, 888, 999)
# =================================================================
#   drink_freq_yr (BD1_11): 1 / 8 (none in past year / never) -> 0,
#                           2-6 -> 1-5, 9 -> NA
#   cigs_per_day (BS3_2):   888 (non-smoker) -> 0, 999 -> NA
#   strength_days_week:     1 = none ... 6 = 5+ days -> 0-5 days, 9 -> NA
#   diagnosed:              9 -> NA
#   family history:         YES / NO / UNKNOWN (8 = no siblings, 9 = don't know)
#     family_history = 1 if father, mother or a sibling was diagnosed, else 0
#     fh_unknown     = 1 if the person did not know (kept as its own category,
#                      NOT imputed: age, sex and biomarkers predict it poorly,
#                      cross-validated AUC 0.54-0.65)

family_history <- function(f1, f2, f3) {
  any_yes   <- (f1 %in% 1) | (f2 %in% 1) | (f3 %in% 1)
  all_known <- !(is.na(f1) | f1 == 9) & !(is.na(f2) | f2 == 9)
  ifelse(any_yes, 1, ifelse(all_known, 0, NA))
}

for (d in names(dts)) {
  x <- dts[[d]]
  if ("drink_freq_yr" %in% names(x))
    x$drink_freq_yr <- case_when(x$drink_freq_yr %in% c(1, 8) ~ 0,
                                 x$drink_freq_yr == 9         ~ NA_real_,
                                 TRUE ~ as.numeric(x$drink_freq_yr) - 1)
  if ("cigs_per_day" %in% names(x))
    x$cigs_per_day <- case_when(x$cigs_per_day == 888 ~ 0,
                                x$cigs_per_day == 999 ~ NA_real_,
                                TRUE ~ as.numeric(x$cigs_per_day))
  x$strength_days_week <- ifelse(x$strength_days_week == 9, NA, x$strength_days_week - 1)
  x$diagnosed          <- ifelse(x$diagnosed == 9, NA, x$diagnosed)
  fh                   <- family_history(x$fh_1, x$fh_2, x$fh_3)
  x$fh_unknown         <- as.integer(is.na(fh))      # 1 = "don't know"
  x$family_history     <- ifelse(is.na(fh), 0, fh)   # 1 = yes, 0 = no or unknown
  x$male               <- ifelse(x$sex == 1, 1, 0)
  dts[[d]] <- x %>% select(-fh_1, -fh_2, -fh_3)
}
cat("Step 4: special codes cleaned\n")


# =================================================================
# STEP 5. Check missing values
# =================================================================
for (d in names(dts)) {
  cat("\n------------------------------\nStep 5: missing values -", d, "\n")
  na_count <- colSums(is.na(dts[[d]]))
  print(na_count[na_count > 0])
}
summary(dts$hypertension)   # full summary for one dataset (change the name to view others)


# -----------------------------------------------------------------
# STEP 5b. Look at WHICH rows are missing and WHY
# -----------------------------------------------------------------
# Missing values come in blocks (someone who skipped the blood test is missing
# every blood value), and the people with missing values tend to be older.
# That is why a single population mode / median is biased: it ignores the
# person's age, BMI, etc. Step 6 predicts each missing value from the rest
# of that person's data instead.
explore_missing <- function(x, d) {
  na <- is.na(x)
  cat("\n------------------------------\nStep 5b:", d, "\n")
  cat("Rows with at least one NA:", sum(rowSums(na) > 0), "of", nrow(x), "\n")
  blood <- intersect(c("fasting_glucose", "hba1c", "total_cholesterol", "hdl",
                       "triglycerides", "ldl"), names(x))
  pa    <- c("vig_pa_hrs_wk", "mod_pa_hrs_wk", "sitting_hrs_day", "strength_days_week")
  cat("Missing ALL blood tests:", sum(apply(na[, blood, drop = FALSE], 1, all)),
      "| missing ALL activity questions:", sum(apply(na[, pa], 1, all)), "\n")
  cat("Mean age of people with vs without a missing value:\n")
  for (v in names(x)[colSums(na) > 0]) {
    m <- na[, v]
    cat(sprintf("  %-22s n missing = %4d | age %.1f vs %.1f\n",
                v, sum(m), mean(x$age[m]), mean(x$age[!m])))
  }
}
for (d in names(dts)) explore_missing(dts[[d]], d)


# -----------------------------------------------------------------
# STEP 5c. Exclude rows where half or more of the clustering features are missing
# -----------------------------------------------------------------
# For these people often only age and sex are known, so any imputed profile
# would be mostly guesswork. They are excluded from that disease's analysis.
max_missing_share <- 0.5

# Original columns behind each disease's clustering features (age / sex are never missing)
feature_sources <- list(
  hypertension  = c("bmi", "waist_cm", "sleep_hrs", "sbp", "dbp", "fasting_glucose",
                    "total_cholesterol", "hdl", "triglycerides", "ldl", "drink_freq_yr",
                    "cigs_per_day", "strength_days_week", "vig_pa_hrs_wk", "mod_pa_hrs_wk",
                    "sitting_hrs_day", "sodium_mg", "family_history"),
  dyslipidaemia = c("bmi", "waist_cm", "sleep_hrs", "sbp", "dbp", "total_cholesterol", "hdl",
                    "triglycerides", "ldl", "cigs_per_day", "drink_freq_yr",
                    "strength_days_week", "vig_pa_hrs_wk", "mod_pa_hrs_wk",
                    "sitting_hrs_day", "family_history"),
  diabetes      = c("bmi", "waist_cm", "sleep_hrs", "sbp", "dbp", "fasting_glucose", "hba1c",
                    "hdl", "triglycerides", "strength_days_week", "vig_pa_hrs_wk",
                    "mod_pa_hrs_wk", "sitting_hrs_day", "family_history", "cigs_per_day")
)

for (d in names(dts)) {
  v          <- feature_sources[[d]]
  share_miss <- rowMeans(is.na(dts[[d]][, ..v]))
  drop       <- share_miss >= max_missing_share
  cat("Step 5c:", d, "- excluded", sum(drop), "rows with >=",
      100 * max_missing_share, "% of features missing; kept", sum(!drop), "\n")
  dts[[d]] <- dts[[d]][!drop, ]
}


# =====================================================================
# STEP 6 — Selected-predictor MICE + regression equation audit
# =====================================================================
# Same idea as before (predict each missing value from the person's other
# variables), but each variable is predicted ONLY from a hand-picked list of
# clinically related predictors (predictor_choices), instead of from every
# other column. The doctor diagnosis is never used as a predictor, so disease
# labels cannot leak into the risk features.
#
# It also fits one ordinary regression per imputed variable on the people who
# have no missing values (the "audit"), and writes its equation out, so the
# imputation models can be reported. These audit equations are diagnostic
# only; MICE itself adds random draws when it fills in values.
#
# R-squared check: every LINEAR model is used, including weak ones
# (R-squared < min_r2 = 0.15). This is safe with PMM: the regression only picks
# donors, and each missing value is copied from a REAL person with a similar
# prediction, so no impossible or "average" values are created. Weak models are
# flagged in the output.
# Set use_hotdeck_if_weak <- TRUE to fill weak-model variables by an age-sex
# hot-deck instead (random donor of the same sex and age band).
#
# Helper columns (height, weight, white blood cells, haemoglobin, uric acid,
# AST, ALT, MCV, urine sodium, urine creatinine) are predictors only. They were
# chosen because they are measured for most people whose target is missing and
# are biologically linked to it (smoking raises WBC / haemoglobin, alcohol raises
# uric acid / liver enzymes / MCV, body size and urine sodium track salt intake).
if (!requireNamespace("mice", quietly = TRUE)) install.packages("mice")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")

# How continuous variables are imputed:
#   "pmm"  = predictive mean matching (DEFAULT): fits a LINEAR REGRESSION, then copies
#            a real value from one of the 5 people with the closest prediction.
#            Never gives impossible values (negative cigarettes) and keeps skewed shapes.
#   "norm" = plain (Bayesian) linear regression: imputed value = regression prediction
#            + random noise. Can give impossible values (in this data: negative activity
#            hours, cigarettes and triglycerides), so they are capped to the observed
#            range after imputation.
continuous_method <- "pmm"
min_r2            <- 0.15     # linear models below this R-squared are flagged as weak
use_hotdeck_if_weak <- FALSE  # FALSE = still use the (weak) regression via PMM

binary_cols  <- c("diagnosed")   # family_history is complete (unknown = own category)
ordinal_cols <- c("drink_freq_yr", "strength_days_week")
not_imputed  <- c("income_quartile_indiv", "education")

# Predictors for each variable. Only those present in a disease dataset are used.
predictor_choices <- list(
  bmi                = c("age", "sex", "waist_cm"),
  waist_cm           = c("age", "sex", "bmi"),
  sbp                = c("age", "sex", "bmi", "dbp"),
  dbp                = c("age", "sex", "bmi", "sbp"),
  fasting_glucose    = c("age", "sex", "bmi", "waist_cm", "hba1c"),
  hba1c              = c("age", "sex", "bmi", "fasting_glucose"),
  total_cholesterol  = c("age", "sex", "ldl", "hdl", "triglycerides"),
  ldl                = c("age", "sex", "total_cholesterol", "hdl", "triglycerides"),
  hdl                = c("age", "sex", "bmi", "triglycerides"),
  triglycerides      = c("age", "sex", "bmi", "waist_cm", "hdl", "fasting_glucose"),
  # salt intake: body size + urine sodium (adjusted for urine creatinine)
  sodium_mg          = c("age", "sex", "bmi", "height_cm", "weight_kg",
                         "urine_sodium", "urine_creatinine"),
  # smoking: raises white blood cells and haemoglobin, lowers HDL, raises TG.
  # (drinking is NOT used: everyone missing smoking also skipped drinking)
  cigs_per_day       = c("age", "sex", "bmi", "wbc", "haemoglobin", "hdl", "triglycerides"),
  # alcohol: raises uric acid, liver enzymes (AST, ALT), red-cell size (MCV), TG.
  # (smoking is NOT used: everyone missing drinking also skipped smoking)
  drink_freq_yr      = c("age", "sex", "bmi", "uric_acid", "ast", "alt", "mcv",
                         "triglycerides", "haemoglobin"),
  strength_days_week = c("age", "sex", "bmi", "waist_cm"),
  vig_pa_hrs_wk      = c("age", "sex", "bmi", "waist_cm"),
  mod_pa_hrs_wk      = c("age", "sex", "bmi", "waist_cm"),
  sitting_hrs_day    = c("age", "sex", "bmi", "waist_cm"),
  # sleep: no available variable predicts it (mental-health items are missing
  # for the same people), so it will fall under the R-squared rule
  sleep_hrs          = c("age", "sex", "bmi"),
  # Diagnosis is imputed ONLY for the post-clustering % diagnosed summary.
  # It is never a k-means input and never a predictor for other variables.
  diagnosed          = c("age", "sex", "bmi", "sbp", "dbp",
                         "fasting_glucose", "hba1c", "ldl", "hdl"),
  # Helper columns (missing for people who skipped the blood / urine test)
  height_cm          = c("age", "sex", "weight_kg"),
  weight_kg          = c("age", "sex", "height_cm", "waist_cm"),
  wbc                = c("age", "sex"),
  haemoglobin        = c("age", "sex"),
  uric_acid          = c("age", "sex", "bmi"),
  ast                = c("age", "sex", "alt"),
  alt                = c("age", "sex", "ast"),
  mcv                = c("age", "sex", "haemoglobin"),
  urine_sodium       = c("age", "sex", "urine_creatinine"),
  urine_creatinine   = c("age", "sex", "weight_kg")
)

# Age-sex hot-deck: each missing value is copied from a random person with the
# same sex and age band who has the value. Used when the linear model is weak.
# (mice calls this when method = "hotdeck"; its predictors must be age and sex)
mice.impute.hotdeck <- function(y, ry, x, wy = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  band <- cut(x[, "age"], c(-Inf, 49, 59, 69, Inf))
  key  <- paste(x[, "sex"], band)
  obs  <- which(ry)
  vapply(which(wy), function(i) {
    pool <- obs[key[obs] == key[i]]
    if (!length(pool)) pool <- obs                      # fallback: anyone
    y[pool[sample.int(length(pool), 1)]]
  }, numeric(1))
}

# Text of a fitted regression equation
format_equation <- function(fit, type, target) {
  b <- stats::coef(fit)
  b <- b[is.finite(b)]
  if (!length(b)) return(NA_character_)
  terms <- paste(ifelse(names(b) == "(Intercept)",
                        sprintf("%+.6g", unname(b)),                 # intercept: number only
                        sprintf("%+.6g*(%s)", unname(b), names(b))), collapse = " ")
  if (type == "logreg") return(paste0("logit(P(", target, "=1)) = ", terms))
  if (type == "polr") {
    # MASS::polr uses logit(P(Y <= k)) = zeta_k - X*beta
    zeta <- paste(sprintf("%s:%.5g", names(fit$zeta), fit$zeta), collapse = "; ")
    return(paste0("logit[P(", target, " <= k)] = zeta_k - (", terms,
                  "); cutpoints {", zeta, "}"))
  }
  paste0(target, " = ", terms)
}

# Diagnostic regression on people with no missing values in target + predictors
fit_audit <- function(x, target, predictors, method, disease) {
  model_type <- c(pmm = "Linear regression", norm = "Linear regression",
                  logreg = "Logistic regression",
                  polr = "Ordered logistic regression")[[method]]
  result <- data.table(disease = disease, variable = target, method = method,
                       model_type = model_type,
                       predictors = paste(predictors, collapse = ", "),
                       formula = paste(target, "~", paste(predictors, collapse = " + ")),
                       n_complete = NA_integer_, equation = NA_character_,
                       status = "not fitted")
  if (!length(predictors)) return(list(row = result, fit = NULL))
  cc <- x[stats::complete.cases(x[, c(target, predictors), drop = FALSE]),
          c(target, predictors), drop = FALSE]
  result$n_complete <- nrow(cc)
  if (nrow(cc) < max(30L, length(predictors) + 10L)) {
    result$status <- "too few complete observations"
    return(list(row = result, fit = NULL))
  }
  fit <- tryCatch({
    f <- stats::reformulate(predictors, response = target)
    if (method %in% c("pmm", "norm")) stats::lm(f, data = cc)
    else if (method == "logreg") stats::glm(f, data = cc, family = stats::binomial())
    else if (method == "polr")   MASS::polr(f, data = cc, method = "logistic", Hess = TRUE)
    else stop("unhandled method")
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    result$status <- paste("fit failed:", conditionMessage(fit))
    return(list(row = result, fit = NULL))
  }
  result$equation <- format_equation(fit, method, target)
  result$status   <- "OK (diagnostic complete-case fit)"
  list(row = result, fit = fit)
}

# -----------------------------------------------------------------
# Regression summaries for the audit models
# -----------------------------------------------------------------
# Coefficient table for any of the three model types
coef_table <- function(fit, disease, target, model_type) {
  if (inherits(fit, "polr")) {
    ct <- stats::coef(summary(fit))                  # slopes AND cutpoints, no p-value
    p  <- 2 * stats::pnorm(-abs(ct[, "t value"]))    # Wald p-value
    out <- data.table(term = rownames(ct), estimate = ct[, "Value"],
                      std_error = ct[, "Std. Error"], statistic = ct[, "t value"],
                      p_value = p,
                      term_type = ifelse(rownames(ct) %in% names(fit$zeta),
                                         "cutpoint", "slope"))
  } else {
    ct <- summary(fit)$coefficients
    out <- data.table(term = rownames(ct), estimate = ct[, 1], std_error = ct[, 2],
                      statistic = ct[, 3], p_value = ct[, 4],
                      term_type = ifelse(rownames(ct) == "(Intercept)", "intercept", "slope"))
  }
  if (inherits(fit, "glm") || inherits(fit, "polr"))
    out[term_type == "slope", odds_ratio := exp(estimate)]  # logistic: odds ratio per unit
  out[, signif := cut(p_value, c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                      labels = c("***", "**", "*", ".", ""))]
  cbind(data.table(disease = disease, variable = target, model_type = model_type), out)
}

# One row of fit statistics per model
fit_stats <- function(fit, disease, target, model_type) {
  if (inherits(fit, "lm") && !inherits(fit, "glm")) {
    sm <- summary(fit)
    fst <- sm$fstatistic
    return(data.table(disease = disease, variable = target, model_type = model_type,
                      n = length(fit$residuals),
                      r_squared = sm$r.squared, adj_r_squared = sm$adj.r.squared,
                      residual_se = sm$sigma,
                      model_p_value = stats::pf(fst[1], fst[2], fst[3], lower.tail = FALSE),
                      aic = stats::AIC(fit), pseudo_r2_mcfadden = NA_real_))
  }
  if (inherits(fit, "glm")) {
    return(data.table(disease = disease, variable = target, model_type = model_type,
                      n = stats::nobs(fit),
                      r_squared = NA_real_, adj_r_squared = NA_real_, residual_se = NA_real_,
                      model_p_value = stats::pchisq(fit$null.deviance - fit$deviance,
                                                    fit$df.null - fit$df.residual,
                                                    lower.tail = FALSE),
                      aic = stats::AIC(fit),
                      pseudo_r2_mcfadden = 1 - fit$deviance / fit$null.deviance))
  }
  # polr: compare with an intercept-only (cutpoints only) model
  null_fit <- MASS::polr(stats::reformulate("1", response = target),
                         data = fit$model, method = "logistic")
  lr <- null_fit$deviance - fit$deviance
  data.table(disease = disease, variable = target, model_type = model_type,
             n = nrow(fit$model),
             r_squared = NA_real_, adj_r_squared = NA_real_, residual_se = NA_real_,
             model_p_value = stats::pchisq(lr, fit$edf - null_fit$edf, lower.tail = FALSE),
             aic = stats::AIC(fit),
             pseudo_r2_mcfadden = 1 - fit$deviance / null_fit$deviance)
}

impute_dataset <- function(x, disease, m = 5L, seed = 123L) {
  x <- as.data.frame(x)
  
  # Factor conversion for categorical outcome models
  for (v in intersect(binary_cols,  names(x))) x[[v]] <- factor(x[[v]])
  for (v in intersect(ordinal_cols, names(x))) x[[v]] <- factor(x[[v]], ordered = TRUE)
  
  # ---- (completed from here) ----
  # Method for each column ("" = nothing to impute)
  meth <- mice::make.method(x)
  meth[intersect(names(meth), binary_cols)]  <- "logreg"
  meth[intersect(names(meth), ordinal_cols)] <- "polr"
  meth[!(names(meth) %in% c(binary_cols, ordinal_cols)) & meth != ""] <- continuous_method
  meth[intersect(names(meth), not_imputed)]  <- ""     # income / education stay NA
  
  # Predictor matrix: start with all zeros, switch on only the chosen predictors
  pred <- matrix(0, ncol(x), ncol(x), dimnames = list(names(x), names(x)))
  used <- list()
  for (v in names(meth)[meth != ""]) {
    p <- predictor_choices[[v]]
    if (is.null(p)) p <- c("age", "sex")                # fallback if a variable is not listed
    p <- setdiff(intersect(p, names(x)), c(v, not_imputed, "male", "diagnosed"))
    pred[v, p] <- 1
    used[[v]]  <- p
  }
  
  # Audit regressions (on the original, un-imputed data)
  audit_res <- lapply(names(used), function(v)
    fit_audit(x, v, used[[v]], meth[[v]], disease))
  names(audit_res) <- names(used)
  audit  <- rbindlist(lapply(audit_res, `[[`, "row"))
  fits   <- lapply(audit_res, `[[`, "fit")
  
  # R-squared check for LINEAR models (weak = R2 < min_r2)
  audit[, r_squared := sapply(fits[variable], function(f)
    if (inherits(f, "lm") && !inherits(f, "glm")) summary(f)$r.squared else NA_real_)]
  audit[, imputed_with := fifelse(method == "logreg", "logistic regression",
                                  fifelse(method == "polr", "ordered logistic regression",
                                          fifelse(!is.na(r_squared) & r_squared >= min_r2,
                                                  paste0("linear regression (", toupper(method), ")"),
                                                  if (use_hotdeck_if_weak)
                                                    "age-sex hot-deck (weak: R2 < 0.15, equation NOT used)"
                                                  else paste0("linear regression (", toupper(method),
                                                              ") - weak: R2 < 0.15"))))]
  for (v in audit[grepl("hot-deck", imputed_with), variable]) {
    meth[v]    <- "hotdeck"
    pred[v, ]  <- 0
    pred[v, c("age", "sex")] <- 1
  }
  
  imp <- mice::mice(x, m = m, method = meth, predictorMatrix = pred,
                    maxit = 10, seed = seed, printFlag = FALSE)
  
  out <- mice::complete(imp, 1)                       # first completed dataset for clustering
  # "norm" can produce impossible values (e.g. negative hours / triglycerides, which
  # would break the log transforms in Step 7). Keep imputed values inside the
  # range actually observed for that variable. Not needed for "pmm".
  if (continuous_method == "norm") {
    for (v in names(meth)[meth == "norm"]) {
      lo <- min(x[[v]], na.rm = TRUE); hi <- max(x[[v]], na.rm = TRUE)
      out[[v]] <- pmin(pmax(out[[v]], lo), hi)
    }
  }
  for (v in intersect(c(binary_cols, ordinal_cols), names(out)))
    out[[v]] <- as.numeric(as.character(out[[v]]))
  list(data = as.data.table(out), imp = imp, audit = audit, fits = fits)
}

imps <- list(); audits <- list(); audit_fits <- list(); was_missing <- list()
for (d in names(dts)) {
  # Guard: Step 6 must run on the un-imputed data from Step 5c, exactly once.
  if (!anyNA(dts[[d]][, !(names(dts[[d]]) %in% not_imputed), with = FALSE]))
    stop("Step 6 has already been run on dts$", d, ". Delete any older copy of ",
         "Step 6 above this one, then re-run the script from Step 1.")
  
  # Which values were missing BEFORE imputing (1 = filled in); Step 8 uses the
  # diagnosis column to compare with people who actually answered it
  was_missing[[d]] <- cbind(ID = dts[[d]]$ID,
                            dts[[d]][, lapply(.SD, function(v) as.integer(is.na(v))),
                                     .SDcols = !"ID"])
  
  ids         <- dts[[d]]$ID                          # ID is not used in the imputation
  res         <- impute_dataset(dts[[d]][, !"ID"], disease = d)
  dts[[d]]    <- cbind(ID = ids, res$data)
  imps[[d]]   <- res$imp                              # keep all 5 imputations
  audits[[d]] <- res$audit
  audit_fits[[d]] <- res$fits
  na_left     <- sum(is.na(dts[[d]][, !(names(dts[[d]]) %in%
                                          not_imputed), with = FALSE]))
  cat("Step 6:", d, "- NA remaining in clustering data:", na_left,
      "| income / education left as NA:", sum(is.na(dts[[d]]$income_quartile_indiv)),
      "/", sum(is.na(dts[[d]]$education)), "\n")
}

# -----------------------------------------------------------------
# Regression audit output
# -----------------------------------------------------------------
# 1) One row per imputed variable per disease: model type, predictors, equation
audit_table <- rbindlist(audits)
cat("\nStep 6: imputation models\n")
print(audit_table[, .(disease, variable, model_type, r_squared = round(r_squared, 3),
                      imputed_with, n_complete)])
fwrite(audit_table, "hn24_imputation_equations.csv")

# 2) Coefficients (estimate, SE, t/z, p-value, odds ratio) and model fit for every model
coef_all  <- list(); stats_all <- list()
for (d in names(audit_fits)) for (v in names(audit_fits[[d]])) {
  fit <- audit_fits[[d]][[v]]
  if (is.null(fit)) next
  mt  <- audit_table[disease == d & variable == v, model_type]
  coef_all[[paste(d, v)]]  <- coef_table(fit, d, v, mt)
  stats_all[[paste(d, v)]] <- fit_stats(fit, d, v, mt)
}
coef_all  <- rbindlist(coef_all,  fill = TRUE)
stats_all <- rbindlist(stats_all, fill = TRUE)
fwrite(coef_all,  "hn24_imputation_coefficients.csv")
fwrite(stats_all, "hn24_imputation_model_fit.csv")

cat("\nStep 6: model fit - LINEAR regressions (continuous variables)\n")
print(stats_all[model_type == "Linear regression",
                .(disease, variable, n, r_squared = round(r_squared, 3),
                  adj_r_squared = round(adj_r_squared, 3),
                  residual_se = round(residual_se, 2), model_p_value = signif(model_p_value, 3))])
cat("\nStep 6: model fit - LOGISTIC regressions (binary and ordinal variables)\n")
print(stats_all[model_type != "Linear regression",
                .(disease, variable, model_type, n,
                  mcfadden_r2 = round(pseudo_r2_mcfadden, 3), aic = round(aic, 1),
                  model_p_value = signif(model_p_value, 3))])

# 3) Full R summary() of every model, grouped by disease and model type
show_regression_summaries <- function(disease = names(audit_fits),
                                      type = c("Linear regression", "Logistic regression",
                                               "Ordered logistic regression")) {
  for (d in disease) for (v in names(audit_fits[[d]])) {
    fit <- audit_fits[[d]][[v]]
    mt  <- audit_table[disease == d & variable == v, model_type]
    if (is.null(fit) || !(mt %in% type)) next
    cat("\n=================================================================\n")
    cat(d, "|", v, "|", mt, "\n")
    cat("Equation:", audit_table[disease == d & variable == v, equation], "\n")
    cat("Used for imputation:", audit_table[disease == d & variable == v, imputed_with], "\n")
    cat("-----------------------------------------------------------------\n")
    print(summary(fit))
    if (inherits(fit, "polr")) {               # polr summary has no p-values
      cat("Wald p-values:\n")
      print(coef_all[disease == d & variable == v, .(term, p_value = signif(p_value, 3), signif)])
    }
  }
}
show_regression_summaries()          # all models
# Examples for a shorter view:
# show_regression_summaries("diabetes")
# show_regression_summaries(type = "Logistic regression")

# -----------------------------------------------------------------
# Merge the three imputed datasets into one CSV (one row per person)
# -----------------------------------------------------------------
# Joined on ID. Columns that are never imputed (ID, sex, age, male, income,
# education) appear once. Every other column gets a suffix for its disease
# dataset (_hyp, _dys, _dia), because each dataset was imputed separately, so
# a shared variable (e.g. bmi) can have slightly different imputed values in
# each. People excluded in Step 5c for one disease have NA in that disease's
# columns.
key_cols <- c("ID", "sex", "age", "male", "income_quartile_indiv", "education")
suffix   <- c(hypertension = "hyp", dyslipidaemia = "dys", diabetes = "dia")

merged <- NULL
for (d in names(dts)) {
  x     <- copy(dts[[d]])
  other <- setdiff(names(x), key_cols)
  setnames(x, other, paste0(other, "_", suffix[[d]]))
  merged <- if (is.null(merged)) x else merge(merged, x, by = key_cols, all = TRUE)
}
setorder(merged, ID)
fwrite(merged, "hn24_imputed_merged.csv")
cat("Step 6: merged dataset saved -", nrow(merged), "people x", ncol(merged), "columns\n")

# Check: imputed values (red) should look like observed values (blue)
# densityplot(imps$hypertension, ~ fasting_glucose + sitting_hrs_day)
# stripplot(imps$hypertension, family_history ~ .imp)


# -----------------------------------------------------------------
# STEP 6b. Normal Q-Q plots: which continuous variables are skewed or heavy-tailed?
# -----------------------------------------------------------------
# Each variable is standardised (z-score), so every panel shares the same
# grey reference line y = x. Points ON the line = normal. Points curving
# ABOVE the line at the right end = long right tail (right-skewed); points
# BELOW at the left end = long left tail; an S-shape = heavy tails on both sides.
# Skewness: |skew| < 0.5 roughly symmetric, 0.5-1 moderate, > 1 strong.
# (drink_freq_yr and strength_days_week are 0-5 categories, so not shown.)

qq_vars <- c("age", "bmi", "waist_cm", "sbp", "dbp", "fasting_glucose", "hba1c",
             "total_cholesterol", "hdl", "triglycerides", "ldl", "sodium_mg",
             "cigs_per_day", "vig_pa_hrs_wk", "mod_pa_hrs_wk", "sitting_hrs_day",
             "sleep_hrs")

# Take each variable from the first disease dataset that contains it
qq_long <- rbindlist(lapply(qq_vars, function(v) {
  d <- names(dts)[sapply(dts, function(x) v %in% names(x))][1]
  data.table(variable = v, value = dts[[d]][[v]])
}))

skewness <- function(x) { z <- (x - mean(x)) / sd(x); mean(z^3) }
kurtosis <- function(x) { z <- (x - mean(x)) / sd(x); mean(z^4) - 3 }   # excess kurtosis

shape <- qq_long[, .(skewness = round(skewness(value), 2),
                     excess_kurtosis = round(kurtosis(value), 2)), by = variable]
shape[, shape := fifelse(abs(skewness) > 1, "Strongly skewed / long tail",
                         fifelse(abs(skewness) > 0.5, "Moderately skewed", "Roughly symmetric"))]
setorder(shape, -skewness)
cat("\nStep 6b: shape of each continuous variable (after imputation)\n")
print(shape)

qq_long <- merge(qq_long, shape[, .(variable, skewness, shape)], by = "variable")
qq_long[, z := (value - mean(value)) / sd(value), by = variable]
qq_long[, panel := factor(paste0(variable, "  (skew ", sprintf("%.1f", skewness), ")"))]
qq_long[, panel := factor(panel, levels = unique(panel[order(-skewness)]))]   # most skewed first
qq_long[, shape := factor(shape, levels = c("Strongly skewed / long tail",
                                            "Moderately skewed", "Roughly symmetric"))]

shape_cols <- c("Strongly skewed / long tail" = "#eb6834",
                "Moderately skewed"           = "#eda100",
                "Roughly symmetric"           = "#2a78d6")

qq_theme <- theme_minimal(base_size = 11) +
  theme(plot.background  = element_rect(fill = "#fcfcfb", colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e4e3df", linewidth = 0.3),
        strip.text       = element_text(face = "bold", size = 9, colour = "#0b0b0b"),
        axis.text        = element_text(colour = "#52514e", size = 8),
        axis.title       = element_text(colour = "#52514e"),
        legend.position  = "top", legend.justification = "left",
        plot.title       = element_text(face = "bold", colour = "#0b0b0b"),
        plot.subtitle    = element_text(colour = "#52514e"))

p_qq <- ggplot(qq_long, aes(sample = z, colour = shape)) +
  geom_abline(slope = 1, intercept = 0, colour = "#8a8984", linewidth = 0.5) +
  stat_qq(size = 0.6, alpha = 0.6) +
  facet_wrap(~ panel, ncol = 6) +
  scale_colour_manual(values = shape_cols, name = NULL,
                      guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  coord_cartesian(ylim = c(-4, 8)) +
  labs(title = "Normal Q-Q plots of continuous variables (after imputation)",
       subtitle = "Standardised values; points on the grey line = normally distributed. Sorted by skewness.",
       x = "Theoretical normal quantile", y = "Observed (standardised)") +
  qq_theme
print(p_qq)

# Same check AFTER the transformations used in Step 7
qq_tr <- rbindlist(list(
  data.table(variable = "log(triglycerides)",  value = log(dts$hypertension$triglycerides)),
  data.table(variable = "log(sodium_mg)",      value = log(dts$hypertension$sodium_mg)),
  data.table(variable = "log1p(vig_pa_hrs_wk)", value = log1p(dts$hypertension$vig_pa_hrs_wk)),
  data.table(variable = "log1p(mod_pa_hrs_wk)", value = log1p(dts$hypertension$mod_pa_hrs_wk))
))
qq_tr[, z := (value - mean(value)) / sd(value), by = variable]
tr_skew <- qq_tr[, .(skewness = round(skewness(value), 2)), by = variable]
cat("\nStep 6b: skewness after the Step 7 transformations\n")
print(tr_skew)
qq_tr <- merge(qq_tr, tr_skew, by = "variable")
qq_tr[, panel := paste0(variable, "  (skew ", sprintf("%.1f", skewness), ")")]

p_qq_tr <- ggplot(qq_tr, aes(sample = z)) +
  geom_abline(slope = 1, intercept = 0, colour = "#8a8984", linewidth = 0.5) +
  stat_qq(size = 0.6, alpha = 0.6, colour = "#2a78d6") +
  facet_wrap(~ panel, ncol = 4) +
  coord_cartesian(ylim = c(-4, 8)) +
  labs(title = "After transformation (as used in Step 7)",
       subtitle = "Log fixes triglycerides and sodium; activity hours stay non-normal because most people report 0",
       x = "Theoretical normal quantile", y = "Observed (standardised)") +
  qq_theme
print(p_qq_tr)

# =================================================================
# STEP 7. K-means clustering (k = 3)
# =================================================================
# Transformations before clustering
#   log_triglycerides, log_sodium  -> right-skewed
#   log_vig_pa, log_mod_pa         -> mostly zeros with a long tail (log1p)
#   sleep_dev_7h                   -> distance from 7 hours (short AND long sleep are risky)
add_transforms <- function(x) {
  x$log_vig_pa   <- log1p(x$vig_pa_hrs_wk)
  x$log_mod_pa   <- log1p(x$mod_pa_hrs_wk)
  x$sleep_dev_7h <- abs(x$sleep_hrs - 7)
  if ("triglycerides" %in% names(x)) x$log_triglycerides <- log(x$triglycerides)
  if ("sodium_mg"     %in% names(x)) x$log_sodium        <- log(x$sodium_mg)
  x
}
for (d in names(dts)) dts[[d]] <- add_transforms(dts[[d]])

# -----------------------------------------------------------------
# Why these log transforms? Distributions before and after
# -----------------------------------------------------------------
# Top row: original scale - long right tails (a few very large values).
# Bottom row: after log / log1p - much closer to symmetric.
# Skewness: 0 = symmetric, > 1 = strongly right-skewed.
# Uses the hypertension dataset (it contains all four variables).
# skewness() is defined in Step 6b

skew_vars <- data.table(
  var   = c("sodium_mg", "triglycerides", "vig_pa_hrs_wk", "mod_pa_hrs_wk",
            "log_sodium", "log_triglycerides", "log_vig_pa", "log_mod_pa"),
  label = c("Sodium (mg/day)", "Triglycerides (mg/dL)",
            "Vigorous activity (hrs/week)", "Moderate activity (hrs/week)",
            "log(sodium)", "log(triglycerides)",
            "log(1 + vigorous hrs)", "log(1 + moderate hrs)"),
  stage = rep(c("Original scale", "After log transform"), each = 4)
)

dist_long <- rbindlist(lapply(seq_len(nrow(skew_vars)), function(i) {
  v <- dts$hypertension[[skew_vars$var[i]]]
  data.table(value = v,
             panel = sprintf("%s\nskewness = %.1f", skew_vars$label[i], skewness(v)),
             stage = skew_vars$stage[i])
}))
dist_long[, panel := factor(panel, levels = unique(panel))]
dist_long[, stage := factor(stage, levels = c("Original scale", "After log transform"))]

p_skew <- ggplot(dist_long, aes(x = value, fill = stage)) +
  geom_histogram(bins = 40, colour = "#fcfcfb", linewidth = 0.2) +
  facet_wrap(~ panel, nrow = 2, scales = "free") +
  scale_fill_manual(values = c("Original scale" = "#2a78d6",
                               "After log transform" = "#eb6834"), name = NULL) +
  labs(title = "Right-skewed variables before and after log transformation",
       subtitle = "Hypertension dataset after imputation, KNHANES 2024, aged 40-74",
       x = NULL, y = "Number of people") +
  theme_minimal(base_size = 11) +
  theme(plot.background    = element_rect(fill = "#fcfcfb", colour = NA),
        legend.position    = "top", legend.justification = "left",
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e4e3df", linewidth = 0.4),
        strip.text         = element_text(face = "bold", colour = "#0b0b0b"),
        axis.text          = element_text(colour = "#52514e"),
        axis.title         = element_text(colour = "#52514e"),
        plot.title         = element_text(face = "bold", colour = "#0b0b0b"),
        plot.subtitle      = element_text(colour = "#52514e"))

print(p_skew)

# Share of people with zero activity (a spike at 0 that no transform can remove)
cat("Zero vigorous activity:", round(100 * mean(dts$hypertension$vig_pa_hrs_wk == 0), 1), "% |",
    "zero moderate activity:", round(100 * mean(dts$hypertension$mod_pa_hrs_wk == 0), 1), "%\n")

# Clustering inputs ("# added" = not in your original list)
cluster_inputs <- list(
  hypertension  = c("age", "male", "bmi", "sleep_dev_7h", "sbp", "dbp",
                    "fasting_glucose", "total_cholesterol", "hdl", "log_triglycerides", "ldl",
                    "drink_freq_yr", "cigs_per_day", "strength_days_week",
                    "log_vig_pa", "log_mod_pa", "sitting_hrs_day",
                    "waist_cm",        # added: central obesity
                    "log_sodium",      # added: salt intake
                    "family_history",  # added
                    "fh_unknown"),     # added: family history "don't know"
  dyslipidaemia = c("age", "male", "bmi", "sleep_dev_7h", "sbp", "dbp",
                    "total_cholesterol", "hdl", "log_triglycerides", "ldl",
                    "cigs_per_day", "drink_freq_yr", "strength_days_week",
                    "log_vig_pa", "log_mod_pa", "sitting_hrs_day",
                    "waist_cm",        # added
                    "family_history",  # added
                    "fh_unknown"),     # added: family history "don't know"
  diabetes      = c("age", "male", "bmi", "sleep_dev_7h", "sbp", "dbp",
                    "fasting_glucose", "hba1c", "strength_days_week",
                    "log_vig_pa", "log_mod_pa", "sitting_hrs_day",
                    "waist_cm",           # added
                    "hdl",                # added: insulin resistance
                    "log_triglycerides",  # added: insulin resistance
                    "family_history",     # added
                    "fh_unknown",         # added: family history "don't know"
                    "cigs_per_day")       # added: smoking raises type 2 diabetes risk
)

# Direction of risk (+1 = higher is worse, -1 = higher is better, 0 = not used for ranking)
risk_direction <- c(
  age = 1, male = 0, bmi = 1, waist_cm = 1, sleep_dev_7h = 1, sbp = 1, dbp = 1,
  fasting_glucose = 1, hba1c = 1, total_cholesterol = 1, hdl = -1, log_triglycerides = 1,
  ldl = 1, drink_freq_yr = 1, cigs_per_day = 1, strength_days_week = -1,
  log_vig_pa = -1, log_mod_pa = -1, sitting_hrs_day = 1, log_sodium = 1,
  family_history = 1,
  fh_unknown = 0          # "don't know" is neither better nor worse for ranking
)

# -----------------------------------------------------------------
# Feature weights for each disease
# -----------------------------------------------------------------
# Without weights, every variable counts equally, so 6 lifestyle variables
# outvote the 2 blood-pressure readings that actually define hypertension.
#
#   3   = defining clinical marker of the disease (diagnostic criterion)
#   2-2.5 = strong, well-established risk factor / core component
#   1.5 = established risk factor
#   1   = moderate / indirect risk factor
#   0.5 = weak, noisy self-report, or overlaps with another input
#
# A weight of w makes the variable count like w variables in the distance
# (columns are multiplied by sqrt(w), because k-means uses squared distance).
feature_weights <- list(
  hypertension = c(
    sbp = 3, dbp = 3,                         # diagnostic criterion (>= 140/90 mmHg)
    age = 2,                                  # arteries stiffen with age
    waist_cm = 1.5, bmi = 1.5,                # obesity: strongest modifiable cause
    log_sodium = 1.5,                         # salt raises blood pressure
    drink_freq_yr = 1.5,                      # alcohol raises blood pressure
    family_history = 1.5, fh_unknown = 0.5,
    fasting_glucose = 1, hdl = 1, log_triglycerides = 1, cigs_per_day = 1,
    ldl = 0.5, total_cholesterol = 0.5,       # not a cause of high BP; overlap with each other
    strength_days_week = 0.5, log_vig_pa = 0.5, log_mod_pa = 0.5,   # self-reported activity
    sitting_hrs_day = 0.5, sleep_dev_7h = 0.5,
    male = 0.5                                # stops sex from splitting the clusters
  ),
  dyslipidaemia = c(
    ldl = 3,                                  # primary treatment target
    hdl = 2.5, log_triglycerides = 2.5,       # core lipid components
    total_cholesterol = 1,                    # mostly LDL + HDL, so overlaps
    age = 1.5, waist_cm = 1.5, family_history = 1.5, fh_unknown = 0.5,
    bmi = 1,
    cigs_per_day = 1,                         # smoking lowers HDL
    drink_freq_yr = 1,                        # alcohol raises triglycerides
    sbp = 0.5, dbp = 0.5,                     # related, but not a cause
    strength_days_week = 0.5, log_vig_pa = 0.5, log_mod_pa = 0.5,
    sitting_hrs_day = 0.5, sleep_dev_7h = 0.5,
    male = 0.5
  ),
  diabetes = c(
    hba1c = 3, fasting_glucose = 3,           # diagnostic criteria (>= 6.5%, >= 126 mg/dL)
    waist_cm = 2,                             # central obesity drives insulin resistance
    bmi = 1.5, age = 1.5, family_history = 1.5, fh_unknown = 0.5,
    log_triglycerides = 1.5, hdl = 1.5,       # markers of insulin resistance
    cigs_per_day = 1,                         # smoking: ~30-40% higher diabetes risk
    sbp = 0.5, dbp = 0.5,
    strength_days_week = 0.5, log_vig_pa = 0.5, log_mod_pa = 0.5,
    sitting_hrs_day = 0.5, sleep_dev_7h = 0.5,
    male = 0.5
  )
)

# Weighted k-means (k = 3) on one dataset; returns cluster ids and Low / Medium / High labels
run_cluster <- function(x, d, seed = 123) {
  vars <- cluster_inputs[[d]]
  w    <- feature_weights[[d]][vars]
  if (any(is.na(w))) stop("Missing weight for: ", paste(vars[is.na(w)], collapse = ", "))
  
  X  <- scale(as.matrix(x[, ..vars]))          # standardise: k-means is distance-based
  Xw <- sweep(X, 2, sqrt(w), "*")               # apply feature weights
  
  set.seed(seed)
  km <- kmeans(Xw, centers = 3, nstart = 50, iter.max = 100)
  
  # Label clusters Low / Medium / High by a WEIGHTED risk score of the cluster centres
  # (centres taken from the standardised, unweighted data)
  centres    <- apply(X, 2, function(col) tapply(col, km$cluster, mean))
  wd         <- w * risk_direction[vars]
  score      <- as.vector(centres[, vars] %*% wd) / sum(abs(wd))
  levels_map <- c("Low", "Medium", "High")[rank(score)]
  
  list(km = km, cluster = km$cluster,
       risk_level = factor(levels_map[km$cluster], levels = c("Low", "Medium", "High")))
}

results <- list()
for (d in names(dts)) {
  cl <- run_cluster(dts[[d]], d)
  dts[[d]]$cluster    <- cl$cluster
  dts[[d]]$risk_level <- cl$risk_level
  results[[d]]        <- cl$km
}
cat("\nStep 7: k-means done\n")


# =================================================================
# STEP 8. Show the results
# =================================================================
# Original units (not log / standardised) for readability
display_name <- c(log_triglycerides = "triglycerides", log_sodium = "sodium_mg",
                  log_vig_pa = "vig_pa_hrs_wk", log_mod_pa = "mod_pa_hrs_wk",
                  sleep_dev_7h = "sleep_hrs")

for (d in names(dts)) {
  vars <- cluster_inputs[[d]]
  show <- unname(ifelse(vars %in% names(display_name), display_name[vars], vars))
  
  tab <- dts[[d]] %>%
    group_by(risk_level) %>%
    summarise(n             = n(),
              pct           = round(100 * n() / nrow(dts[[d]]), 1),
              pct_diagnosed = round(100 * mean(diagnosed), 1),   # check: not an input
              across(all_of(show), ~ round(mean(.x), 2)),
              .groups = "drop") %>%
    as.data.frame()
  
  cat("\n==================================================\n")
  cat("Step 8:", d, "\n")
  cat("Between-cluster SS / total SS (weighted):",
      round(results[[d]]$betweenss / results[[d]]$totss, 3), "\n")
  out <- t(tab[, -1]); colnames(out) <- as.character(tab$risk_level)
  print(out)
  
  # Compare the risk levels with people who ARE / are NOT doctor-diagnosed.
  # Uses only people who actually answered the diagnosis question (not imputed).
  # Note: diagnosed people are often on medication, so their BP / lipids /
  # glucose can be LOWER than their untreated level.
  answered_dx <- was_missing[[d]]$diagnosed[match(dts[[d]]$ID, was_missing[[d]]$ID)] == 0
  dx <- dts[[d]][answered_dx, ] %>%
    mutate(group = ifelse(diagnosed == 1, "Diagnosed", "Not diagnosed")) %>%
    group_by(group) %>%
    summarise(n             = n(),
              pct           = round(100 * n() / sum(answered_dx), 1),
              pct_diagnosed = round(100 * mean(diagnosed), 1),
              across(all_of(show), ~ round(mean(.x), 2)),
              .groups = "drop") %>%
    as.data.frame()
  dx_out <- t(dx[, -1]); colnames(dx_out) <- dx$group
  comp   <- cbind(out, dx_out[rownames(out), c("Diagnosed", "Not diagnosed")])
  cat("\nMeans by risk level vs diagnosed / not diagnosed (",
      sum(answered_dx), "people with a known diagnosis ):\n")
  print(comp)
  fwrite(data.table(variable = rownames(comp), comp),
         paste0("hn24_", d, "_profile_vs_diagnosed.csv"))
  
  fwrite(dts[[d]], paste0("hn24_", d, "_risk_levels.csv"))
}

clust_hypertension  <- dts$hypertension
clust_dyslipidaemia <- dts$dyslipidaemia
clust_diabetes      <- dts$diabetes