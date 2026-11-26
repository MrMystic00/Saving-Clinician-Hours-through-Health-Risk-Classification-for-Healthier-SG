library(data.table)
library(rpart)
library(rpart.plot)

raw <- fread("hn24_imputed_merged_with_tiers.csv")
raw <- raw[age <= 74]                       # team scope: residents aged 40-74

# ---- 1. Modelling table ----------------------------------------------
# Inputs a nurse can get on the spot: age, sex, BP machine, tape measure,
# scale (BMI), smoking question, finger-prick HbA1c
d <- data.table(age = raw$age, male = raw$male,
                sbp = raw$sbp_hyp, dbp = raw$dbp_hyp, bmi = raw$bmi_hyp,
                waist = raw$waist_cm_hyp, cigs = raw$cigs_per_day_hyp,
                hba1c = raw$hba1c_dia,
                diag_hyp = raw$diagnosed_hyp, diag_dys = raw$diagnosed_dys,
                diag_dia = raw$diagnosed_dia,
                tier_hyp = raw$tier_hyp, tier_dia = raw$tier_dia)
d <- d[!is.na(tier_hyp) & !is.na(tier_dia) & !is.na(hba1c) & !is.na(sbp) &
       !is.na(diag_hyp) & !is.na(diag_dys) & !is.na(diag_dia)]

# Target: HIGH if High for hypertension or diabetes, otherwise NOT HIGH
d[, triage := "NotHigh"]
d[tier_hyp == "High" | tier_dia == "High", triage := "High"]
d[, triage := factor(triage, levels = c("NotHigh", "High"))]
d[, diag_any := diag_hyp == 1 | diag_dys == 1 | diag_dia == 1]
print(nrow(d))
print(table(d$triage))

# ---- 2. Sense check: does HIGH contain more people with disease? --------
print(d[, .(n = .N, share = round(100 * .N / nrow(d), 1),
            diag_hyp = round(100 * mean(diag_hyp), 1),
            diag_dia = round(100 * mean(diag_dia), 1),
            diag_any = round(100 * mean(diag_any), 1)), keyby = triage])

# ---- 3. Train / test split (70 / 30) -------------------------------------
set.seed(2406)
train_rows <- sample(nrow(d), round(0.7 * nrow(d)))
trainset <- d[train_rows]
testset  <- d[-train_rows]

# ---- 4. Grow the maximal tree on the trainset -----------------------------
set.seed(2406)                                   # for randomisation in 10-fold CV
m1 <- rpart(triage ~ age + male + sbp + dbp + bmi + waist + cigs + hba1c,
            data = trainset, method = "class", control = rpart.control(cp = 0))
printcp(m1)
plotcp(m1)

# ---- 5. Optimal tree: 1 SE rule --------------------------------------------
# CV error cap = minimum CV error + its standard error

CVerror.cap <- 0.31061 + 0.017315
CVerror.cap

cp1 <- sqrt(0.00757576 * 0.00865801)

cart_opt <- prune(m1, cp = cp1)
printcp(cart_opt)
cart_opt

# ---- 6. Nurse tree: largest tree with at most 5 splits ----
# From printcp(m1): row 4 has 4 splits (row 5 has 7, too many for a checklist)
# cp = geometric mean of the CP of row 4 and row 3
cp_nurse <- sqrt(0.00865801 * 0.04329004)
cart_nurse <- prune(m1, cp = cp_nurse)
printcp(cart_nurse)                 # should end at 4 splits

# ---- 7. Test both trees on the testset ---------------------------------------
# Sensitivity = share of HIGH people the tree catches (the safety measure)
testset[, pred_opt   := predict(cart_opt,   newdata = testset, type = "class")]
testset[, pred_nurse := predict(cart_nurse, newdata = testset, type = "class")]

print(table(Actual = testset$triage, Predicted = testset$pred_opt))
print(paste("Optimal tree |",
            "accuracy",    round(mean(testset$pred_opt == testset$triage), 3),
            "| sensitivity", round(mean(testset[triage == "High", pred_opt] == "High"), 3),
            "| specificity", round(mean(testset[triage == "NotHigh", pred_opt] == "NotHigh"), 3),
            "| referred",    round(mean(testset$pred_opt == "High"), 3)))

print(table(Actual = testset$triage, Predicted = testset$pred_nurse))
print(paste("Nurse tree |",
            "accuracy",    round(mean(testset$pred_nurse == testset$triage), 3),
            "| sensitivity", round(mean(testset[triage == "High", pred_nurse] == "High"), 3),
            "| specificity", round(mean(testset[triage == "NotHigh", pred_nurse] == "NotHigh"), 3),
            "| referred",    round(mean(testset$pred_nurse == "High"), 3)))

# Baseline: call everyone NOT HIGH
print(paste("Baseline accuracy:", round(mean(testset$triage == "NotHigh"), 3)))

# ---- 8. The nurse checklist ---------------------------------------------------
print(cart_opt)
print(cart_nurse)
cart_opt$variable.importance
rpart.plot(cart_opt, nn = T, digits = 3, box.palette = c("#D9D9D9", "#C0392B"), shadow.col = 0, main = "Stage 1 triage: optimal tree (1 SE rule)")

rpart.plot(cart_nurse, nn = T, digits = 3,
           box.palette = c("#D9D9D9", "#C0392B"),  
           shadow.col = 0,
           main = "Stage 1 triage: nurse tree")

# ---- 9. Who gets referred? (nurse tree, testset) -------------------------------
print(testset[, .(n = .N, share = round(100 * .N / nrow(testset), 1),
                  diag_any = round(100 * mean(diag_any), 1),
                  hyp_or_dia_diagnosed = round(100 * mean(diag_hyp == 1 | diag_dia == 1), 1)),
              keyby = pred_nurse])
# Truly HIGH people that the nurse tree refers, never diagnosed with hypertension or diabetes
referred_high <- testset[pred_nurse == "High" & triage == "High"]
print(paste("Referred HIGH people never diagnosed with hypertension or diabetes:",
            round(100 * mean(referred_high$diag_hyp == 0 & referred_high$diag_dia == 0), 1), "%"))

# ---- 9B. Breakdown of the people the nurse tree refers (testset) ---------------
referred <- testset[pred_nurse == "High"]
print(paste("Referred to a doctor:", nrow(referred)))
print(paste("Model errors (actually Not High):", nrow(referred[triage == "NotHigh"])))
print(paste("Correctly referred (actually High):", nrow(referred[triage == "High"])))

correct <- referred[triage == "High"]
print(paste("Already diagnosed with hypertension or diabetes:",
            nrow(correct[diag_hyp == 1 | diag_dia == 1])))
undiagnosed <- correct[diag_hyp == 0 & diag_dia == 0]
print(paste("Never diagnosed:", nrow(undiagnosed)))

# Readings already in the diagnostic range: BP >= 140/90 mmHg or HbA1c >= 6.5%
print(paste("Readings in the disease range (likely undiagnosed disease):",
            nrow(undiagnosed[sbp >= 140 | dbp >= 90 | hba1c >= 6.5])))
print(paste("Readings not yet in the disease range (high risk only):",
            nrow(undiagnosed[sbp < 140 & dbp < 90 & hba1c < 6.5])))

# Diagnosed OR readings in the disease range: referred vs not referred
print(testset[, .(n = .N,
                  diagnosed_or_disease_range = round(100 * mean(diag_hyp == 1 | diag_dia == 1 |
                                                                  sbp >= 140 | dbp >= 90 | hba1c >= 6.5), 1)),
              keyby = pred_nurse])

# ---- 10. Referral rate by age band (to scale up with Singapore population) ------
d[, pred_nurse := predict(cart_nurse, newdata = d, type = "class")]
d[, band := cut(age, breaks = c(39, 49, 59, 69, 74), labels = c("40-49", "50-59", "60-69", "70-74"))]
print(d[, .(n = .N, referred = round(100 * mean(pred_nurse == "High"), 1), actual_high = round(100 * mean(triage == "High"), 1)), keyby = band])