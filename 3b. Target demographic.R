# Goal: Singapore residents by age band, aligned to the NPHS reporting bands
## So prevalence numbers can be applied

# Data: SingaporePopulationDataGovSG.csv

library(data.table)
library(ggplot2)

choose.dir()
setwd()

# Year choices
BASE_YEAR <- 2024                    # NPHS 2024
YEARS <- c(2023, 2024, 2025, 2026)

# ---- 1. Read the file ----------------------------------------
dt <- fread("SingaporePopulationDataGovSG.csv", header = TRUE,
            select = c("DataSeries", as.character(2010:2026)))

# Label the three blocks by row number, and number the rows inside each block
dt[1:97,    sex := "Total"]
dt[98:194,  sex := "Male"]
dt[195:291, sex := "Female"]
dt[, row := 1:.N, by = sex]

# Layout check: the first row of each block must be the expected header
print(dt[row == 1, DataSeries])
## Output should be:
## "Singapore Residents" "Singapore Male Residents" "Singapore Female Residents"

# ---- 2. Wide (one column per year) -> long (one row per year) --
long <- melt(dt, id.vars = c("DataSeries", "sex", "row"),
             measure.vars = as.character(2010:2026),
             variable.name = "year", value.name = "pop")
long[, year := as.numeric(as.character(year))]

# Age of each row: rows 2-91 are ages 0-89, row 97 is "90 and over" (coded 90).
# Row 1 (block total) and rows 92-96 (subtotals) have no age, so they are
# left out of the band sums and never double counted.
long[row >= 2 & row <= 91, age := row - 2]
long[row == 97, age := 90]

# ---- 3. Single ages -> NPHS-aligned bands ----------------------
long[age < 18, band := "Under 18"]
long[age >= 18 & age < 30, band := "18-29"]
long[age >= 30 & age < 40, band := "30-39"]
long[age >= 40 & age < 50, band := "40-49"]
long[age >= 50 & age < 60, band := "50-59"]
long[age >= 60 & age < 70, band := "60-69"]
long[age >= 70 & age < 75, band := "70-74"]
long[age >= 75, band := "75+"]

pop_band <- long[!is.na(age), .(pop = sum(pop)), by = .(year, sex, band)]

# ---- 4. Tables for the report ----------------------------------
band_levels <- c("Under 18", "18-29", "30-39", "40-49", "50-59",
                 "60-69", "70-74", "75+")
pop_band[, band := factor(band, levels = band_levels)]
inscope <- c("18-29", "30-39", "40-49", "50-59", "60-69", "70-74")   # NPHS covers 18-74

# Table 1: base year, by sex
tab_base <- dcast(pop_band[year == BASE_YEAR], band ~ sex, value.var = "pop")
setorder(tab_base, band)
denom  <- tab_base[band %in% inscope, sum(Total)]   # residents aged 18-74
allres <- tab_base[, sum(Total)]                    # all residents
tab_base[band %in% inscope, pct_of_18_74 := round(100 * Total / denom, 1)]
tab_base[, pct_of_residents := round(100 * Total / allres, 1)]
tab_base[, in_NPHS_scope := "No"]
tab_base[band %in% inscope, in_NPHS_scope := "Yes"]
print(tab_base)

# Table 2: total residents, one column per year
tab_years <- dcast(pop_band[sex == "Total" & year %in% YEARS], band ~ year, value.var = "pop")
setorder(tab_years, band)
print(tab_years)

# ---- 6. Ageing trend and linear regression (Category A) --------
# Rows 1, 92 and 94 of the Total block = all residents, 65+, 75+
tr <- dcast(long[sex == "Total" & row %in% c(1, 92, 94)], year ~ row, value.var = "pop")
setnames(tr, c("year", "residents", "p65", "p75"))
tr[, share65 := 100 * p65 / residents]
tr[, share75 := 100 * p75 / residents]
print(tr[year %in% c(2010, 2015, 2020, 2023, 2024, 2025, 2026)])

# Growth is accelerating, so a straight line understates the recent pace
# Use it to describe the trend, not to forecast far ahead.
fit <- lm(share65 ~ year, data = tr)
summary(fit)

# ---- 7. Bar chart: residents by age band (all residents) ---------
plot_df <- pop_band[year == BASE_YEAR & sex == "Total" &
                      band %in% c("Under 18", "18-29", "30-39", "40-49", "50-59",
                                  "60-69", "70-74", "75+")]
plot_df[, band := factor(band, levels = rev(c("Under 18", "18-29", "30-39", "40-49",
                                              "50-59", "60-69", "70-74", "75+")))]
plot_df[, scope := "Inside NPHS (18-74)"]
plot_df[band %in% c("Under 18", "75+"), scope := "Outside NPHS (not surveyed)"]

p1 <- ggplot(plot_df, aes(x = band, y = pop / 1000, fill = scope)) +
  geom_col() +
  geom_text(aes(label = round(pop / 1000)), nudge_y = 40) +   # value ('000) beside each bar
  expand_limits(y = 800) +                                    # room for the labels
  coord_flip() +
  scale_fill_manual(values = c("Inside NPHS (18-74)" = "steelblue",
                               "Outside NPHS (not surveyed)" = "grey")) +
  labs(title = paste("Singapore residents by age band, end-June", BASE_YEAR),
       x = "Age band (years)", y = "Residents ('000)", fill = "NPHS coverage") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p1)
ggsave(paste0("fig_population_by_band_", BASE_YEAR, ".png"), width = 8, height = 4.5)

tr_long <- melt(tr, id.vars = "year", measure.vars = c("share65", "share75"),
                variable.name = "group", value.name = "share")
tr_long[, group := as.character(group)]
tr_long[group == "share65", group := "Aged 65+"]
tr_long[group == "share75", group := "Aged 75+"]

# Lines = actual shares; straight line = the linear regression for the 75+ share
p2 <- ggplot(tr_long, aes(x = year, y = share, color = group)) +
  geom_line() +
  geom_smooth(data = tr_long[group == "Aged 65+"], method = lm) +
  labs(title = "Share of Singapore residents aged 65+ and 75+, end-June",
       subtitle = "Straight line with shaded band = linear regression fit for the 65+ share",
       x = "Year", y = "% of residents", color = "Age group") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p2)
ggsave("fig_ageing_trend.png", width = 7, height = 4)