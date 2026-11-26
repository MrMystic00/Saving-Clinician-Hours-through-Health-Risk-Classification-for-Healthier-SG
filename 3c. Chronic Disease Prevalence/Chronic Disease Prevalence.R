# Sources: NPHS 2024 and SingStat 2026

library(data.table)
library(ggplot2)

# FACET CHART 1
# Chronic Disease Prevalence by Age Group (%)

disease.dt <- data.frame(
  AgeGroup = rep(c("40-49", "50-59", "60-69", "70-74"),3),
  Disease = rep(c("Hypertension","High Cholesterol","Diabetes"),each=4),
  Prevalence = c(29.0, 48.1, 59.2, 69.4,
               24.1, 40.5, 54.2, 58.8,
               5.4, 12.9, 21.9, 22.0))

disease.dt$AgeGroup <- factor(disease.dt$AgeGroup, ordered = T, levels = c("40-49","50-59","60-69","70-74"))
disease.dt$Disease <- factor(disease.dt$Disease)

ggplot(data=disease.dt, aes(x = AgeGroup, y = Prevalence, fill = AgeGroup)) +
  geom_bar(stat ="identity") +
  geom_text(aes(label = paste0(Prevalence, "%")), vjust = -0.5) +
  facet_grid(.~Disease) +
  scale_fill_manual(values = c("coral", "lightsalmon" , "sandybrown", "indianred")) +
  labs(title = "Chronic Disease Prevalence by Age Group (%)",
       x = "Age Group", y = "Prevalence (%)", fill = "Age Group") +
  ylim(0, 80)

# FACET CHART 2
# Distribution of Estimated Chronic Disease Cases by Age Group (%)

# Import SingStat population data
setwd("~/Downloads")
population.dt <- fread("SingaporePopulationDataGovSG.csv")

# Identify column names
names(population.dt)

# Calculate 2024 population by age group
population <- c(
  sum(population.dt[[4]][43:52]),  # 40 to 49
  sum(population.dt[[4]][53:62]),  # 50 to 59
  sum(population.dt[[4]][63:72]),  # 60 to 69
  sum(population.dt[[4]][73:77])   # 70 to 74
)

print(population)

# Add population figures
disease2.dt <- disease.dt

# Calculate estimated cases
disease2.dt$Population <- rep(population, 3)
disease2.dt$EstimatedCases <- disease2.dt$Population * disease2.dt$Prevalence / 100

View(disease2.dt)

# Calculate percentage of estimated cases by age group (adding the percentage column into disease2.dt)
disease2.dt$CasePercent <- c(
  disease2.dt$EstimatedCases[1:4] / sum(disease2.dt$EstimatedCases[1:4]) * 100,
  disease2.dt$EstimatedCases[5:8] / sum(disease2.dt$EstimatedCases[5:8]) * 100,
  disease2.dt$EstimatedCases[9:12] / sum(disease2.dt$EstimatedCases[9:12]) * 100
)

View(disease2.dt)

ggplot(data = disease2.dt, aes(x = AgeGroup, y = CasePercent, fill = AgeGroup)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(CasePercent, 1), "%")), vjust = -0.5) +
  facet_grid(.~Disease) +
  scale_fill_manual(values = c("coral", "lightsalmon" , "sandybrown", "indianred")) +
  labs(title = "Distribution of Estimated Chronic Disease Cases by Age Group (%)",
       x = "Age Group", y = "Share of Estimated Cases (%)", fill = "Age Group") +
  ylim(0, 50)

# Total cases by each disease
# hypertension
sum(disease2.dt$EstimatedCases[1:4])

# high cholesterol
sum(disease2.dt$EstimatedCases[5:8])

# diabetes
sum(disease2.dt$EstimatedCases[9:12])



