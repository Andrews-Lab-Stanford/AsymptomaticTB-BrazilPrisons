# ASYMPTOMATIC TRANSMISSION ANALAYSES             
# User: Esther Jung                                                    
# Last Updated: April 25, 2025                                  

rm(list = ls())

# DATA LOADING ================================================================

library(tidyverse)
library(lme4)
library(rstanarm)
library(broom.mixed)
library(ggplot2)
library(mgcv)
library(corrplot)
library(gridExtra)
library(data.table)

data <- read.csv('~/Google Drive/My Drive/Subclinical TB/Data/asymptomatic_data.csv')
cases <- read.csv('~/Google Drive/My Drive/Subclinical TB/Data/cases.csv')

# DATA CLEANING ================================================================

data <- data %>%
  mutate(
    QFT_max = pmax(qft_tb1, qft_tb2),
    QFT_result_over_.35 = as.integer(QFT_max >= 0.35),
    QFT_result_over_.70 = as.integer(QFT_max >= 0.70),
    QFT_result_over_2.00 = as.integer(QFT_max >= 2.00),
    QFT_result_over_4.00 = as.integer(QFT_max >= 4.00),
    prison = gsub("[^a-zA-Z]", "", cell),
    GeneXpert.semiquant = case_when(
      xpert_status_case %in% 1:2 ~ "med/high",
      xpert_status_case %in% 3:4 ~ "low/verylow",
      TRUE ~ "negative"
    ),
    cell_entry = as.Date(cell_entry, origin = "1899-12-30"),
    cell_entry_case = as.Date(cell_entry_case, origin = "1899-12-30"),
    xpert_date_case = as.Date(xpert_date_case, origin = "1899-12-30"),
    date_qft_t1 = as.Date(date_qft, origin = "1899-12-30"),
    contact_duration = as.numeric(xpert_date_case - cell_entry),  
    case_duration = as.numeric(xpert_date_case - cell_entry_case),
    max_exposure = pmin(contact_duration, case_duration),
    cell_density = total_cellmates / m2
  ) %>%
  mutate(
    race_group = case_when(race == '1' ~ 'white',
                           race == '2' ~ 'black',
                           race == '3' ~ 'brown',
                           race == '4' ~ 'yellow',
                           race == '5' ~ 'indigenous',
                           TRUE ~ 'white'),
    schooling_group = case_when(
      schooling %in% 1:7 ~ 'up to 8 years',
      schooling >= 8 ~ 'More than 8 years',
      TRUE ~ 'NA'
    ),
    n_pdls_cell = case_when(
      total_cellmates <= 10 ~ '10 or less',
      total_cellmates <= 20 ~ '11-20',
      total_cellmates <= 40 ~ '21-40',
      total_cellmates >= 41 ~ '41 to 60',
      TRUE ~ 'NA'
    ),
    n_pdls_dic = case_when(
      total_cellmates <= 20 ~ '20 or less',
      total_cellmates > 20 ~ '21 or more',
      TRUE ~ 'NA'
    ),
    drugs = as.integer(rowSums(across(c(marijuana, cocaine, crack, heroin, glue_sniffing, coca_paste, hashish, injectables)) == 1, na.rm = TRUE) > 0)
  ) 

no_prev <- data %>% filter(prev_prison == 0) # Restrict to those without previous incarceration

cases <- cases %>%
  mutate(schooling_group = case_when(
    n11 %in% 1:7 ~ 'up to 8 years',
    n11 >= 8 ~ 'More than 8 years',
    TRUE ~ 'NA'
  ), .after = n11)

# Table 1 ======================================================================

pairwise_fisher <- function(var, group, data) {
  combos <- combn(unique(na.omit(data[[group]])), 2, simplify = FALSE)
  
  map_dfr(combos, function(groups) {
    subset_df <- filter(data, .data[[group]] %in% groups, !is.na(.data[[var]]))
    tbl <- table(subset_df[[var]], subset_df[[group]])
    p_val <- ifelse(all(dim(tbl) == c(2, 2)), fisher.test(tbl)$p.value, NA)
    tibble(variable = var, group1 = groups[1], group2 = groups[2], p_value = p_val)
  })
}

binary_vars <- c("schooling_group", "smoking_status", "alcohol", "drugs", "bcg", "previous_ctt_tb", 
                 "n_pdls_cell", "n_pdls_dic")
table1_results_binary <- map_dfr(binary_vars, ~ pairwise_fisher(.x, "tb_exposition", no_prev))

CrossTable(no_prev$tb_exposition) # n

aggregate(no_prev$age, list(no_prev$tb_exposition), FUN = summary) # Median age
wilcox.test(no_prev$age[no_prev$tb_exposition == 1], no_prev$age[no_prev$tb_exposition == 2])
wilcox.test(no_prev$age[no_prev$tb_exposition == 1], no_prev$age[no_prev$tb_exposition == 0])
wilcox.test(no_prev$age[no_prev$tb_exposition == 2], no_prev$age[no_prev$tb_exposition == 0])

CrossTable(x = no_prev$race_group, y = no_prev$tb_exposition) # Race
fisher.test(matrix(c(79, 121, 132 - 79, 224 - 121), nrow = 2))$p.value # mixed
fisher.test(matrix(c(38, 66, 132 - 38, 224 - 66), nrow = 2))$p.value   # white
fisher.test(matrix(c(13, 27, 132 - 13, 224 - 27), nrow = 2))$p.value   # black
fisher.test(matrix(c(2, 7, 132 - 2, 224 - 7), nrow = 2))$p.value       # indigenous
fisher.test(matrix(c(0, 3, 132 - 0, 224 - 3), nrow = 2))$p.value       # asian
fisher.test(matrix(c(79, 187, 132 - 79, 330 - 187), nrow = 2))$p.value  # mixed
fisher.test(matrix(c(38, 116, 132 - 38, 330 - 116), nrow = 2))$p.value  # white
fisher.test(matrix(c(13, 20, 132 - 13, 330 - 20), nrow = 2))$p.value    # black
fisher.test(matrix(c(2, 5, 132 - 2, 330 - 5), nrow = 2))$p.value        # indigenous
fisher.test(matrix(c(0, 2, 132 - 0, 330 - 2), nrow = 2))$p.value.       # asian
fisher.test(matrix(c(121, 187, 224 - 121, 330 - 187), nrow = 2))$p.value  # mixed
fisher.test(matrix(c(66, 116, 224 - 66, 330 - 116), nrow = 2))$p.value    # white
fisher.test(matrix(c(27, 20, 224 - 27, 330 - 20), nrow = 2))$p.value.     # black
fisher.test(matrix(c(7, 5, 224 - 7, 330 - 5), nrow = 2))$p.value.         # indigenous
fisher.test(matrix(c(3, 2, 224 - 3, 330 - 2), nrow = 2))$p.value          # asian

CrossTable(x = no_prev$schooling_group, y = no_prev$tb_exposition) # Education
CrossTable(no_prev$smoking_status,no_prev$tb_exposition) # Smoking
CrossTable(no_prev$alcohol,no_prev$tb_exposition) # Alcohol
CrossTable(no_prev$drugs,no_prev$tb_exposition) # Drug use
CrossTable(no_prev$bcg, no_prev$tb_exposition) # BCG vaccine
CrossTable(no_prev$previous_ctt_tb, no_prev$tb_exposition) # Knows someone with TB
CrossTable(no_prev$n_pdls_dic, no_prev$tb_exposition) # Number of individuals in cell

aggregate(no_prev$months_prison, list(no_prev$tb_exposition), FUN = summary) # Months incarcerated
wilcox.test(no_prev$months_prison[no_prev$tb_exposition == 1], no_prev$months_prison[no_prev$tb_exposition == 2])
wilcox.test(no_prev$months_prison[no_prev$tb_exposition == 1], no_prev$months_prison[no_prev$tb_exposition == 0])
wilcox.test(no_prev$months_prison[no_prev$tb_exposition == 2], no_prev$months_prison[no_prev$tb_exposition == 0])

CrossTable(no_prev$n_pdls_cell, no_prev$tb_exposition) # Number of individuals in cell
fisher.test(matrix(c(12, 33, 132 - 12, 224 - 33), nrow = 2))$p.value
fisher.test(matrix(c(12, 68, 132 - 12, 330 - 68), nrow = 2))$p.value
fisher.test(matrix(c(33, 68, 224-33, 330-68), nrow = 2))$p.value

fisher.test(matrix(c(56, 46, 132 - 56, 224 - 46), nrow = 2))$p.value
fisher.test(matrix(c(56, 203, 132 - 56, 330 - 203), nrow = 2))$p.value
fisher.test(matrix(c(46, 203, 224-46, 330-203), nrow = 2))$p.value

fisher.test(matrix(c(23, 83, 132 - 23, 224 - 83), nrow = 2))$p.value
fisher.test(matrix(c(23, 43, 132 - 23, 330 - 43), nrow = 2))$p.value
fisher.test(matrix(c(83, 43, 224-83, 330-43), nrow = 2))$p.value

fisher.test(matrix(c(23, 83, 132 - 23, 224 - 83), nrow = 2))$p.value
fisher.test(matrix(c(23, 43, 132 - 23, 330 - 43), nrow = 2))$p.value
fisher.test(matrix(c(83, 43, 224-83, 330-43), nrow = 2))$p.value

fisher.test(matrix(c(41, 62, 132 - 41, 224 - 62), nrow = 2))$p.value
fisher.test(matrix(c(41, 16, 132 - 41, 330 - 16), nrow = 2))$p.value
fisher.test(matrix(c(64, 16, 224-62, 330-16), nrow = 2))$p.value

CrossTable(no_prev$prison, no_prev$tb_exposition) # Prison
fisher.test(matrix(c(43, 40, 132 - 43, 224 - 40), nrow = 2))$p.value
fisher.test(matrix(c(43, 51, 132 - 43, 330 - 51), nrow = 2))$p.value
fisher.test(matrix(c(40, 51, 224-40, 330-51), nrow = 2))$p.value

fisher.test(matrix(c(36,57, 132 - 36, 224 - 57), nrow = 2))$p.value
fisher.test(matrix(c(36,201, 132 - 36, 330 - 201), nrow = 2))$p.value
fisher.test(matrix(c(57, 201, 224-57, 330-201), nrow = 2))$p.value

fisher.test(matrix(c(53,127, 132 - 53, 224 - 127), nrow = 2))$p.value
fisher.test(matrix(c(53, 78, 132 - 53, 330 - 78), nrow = 2))$p.value
fisher.test(matrix(c(40, 78, 224-127, 330-78), nrow = 2))$p.value

# Table 2 ======================================================================
CrossTable(cases$tb_type) # n 
fisher.test(matrix(c(104, 186, 186, 104), nrow = 2))$p.value

aggregate(cases$age, list(cases$tb_type), FUN = summary) # Median age
wilcox.test(cases$age[cases$tb_type == 1], cases$age[cases$tb_type == 0])

CrossTable(x = cases$race, y = cases$tb_type) # Race
fisher.test(matrix(c(55, 130, 104-55, 186-130), nrow = 2))$p.value
fisher.test(matrix(c(26, 26, 104-26, 186-26), nrow = 2))$p.value
fisher.test(matrix(c(20, 24, 104-20, 186-24), nrow = 2))$p.value
fisher.test(matrix(c(1, 6, 104-1, 186-6), nrow = 2))$p.value
fisher.test(matrix(c(2, 0, 104-2, 186-0), nrow = 2))$p.value

CrossTable(x = cases$schooling_group, y = cases$tb_type, fisher = TRUE) # Education
CrossTable(x = cases$smoking_status, y = cases$tb_type, fisher = TRUE) # Smoking
CrossTable(x = cases$alcohol, y = cases$tb_type, fisher = TRUE) # Alcohol
CrossTable(x = cases$drugs, y = cases$tb_type, fisher = TRUE) # Drug use
CrossTable(x = cases$bcg, y = cases$tb_type, fisher = TRUE) # BCG Vaccine
CrossTable(x = cases$n26, y = cases$tb_type, fisher = TRUE) # Knows someone with TB
CrossTable(x = cases$prev_prison, y = cases$tb_type, fisher = TRUE) # Previous incarceration

CrossTable(x = cases$cell_prison, y = cases$tb_type) # Prison
fisher.test(matrix(c(56, 104, 104-56, 186-104), nrow = 2))$p.value
fisher.test(matrix(c(21, 34, 104-21, 186-34), nrow = 2))$p.value
fisher.test(matrix(c(27, 48, 104-27, 186-48), nrow = 2))$p.value

# GLMMs: QFT Positivity ========================================================

### TB exposure models
run_stan_model <- function(outcome, exposure, data) {
  formula <- as.formula(paste0(outcome, " ~ relevel(factor(tb_exposition), ref = '0') + prison + cell_density + age + months_prison + (1 | cell)"))
  fit <- stan_glmer(formula, family = binomial, data = data, 
                    chains = 4, iter = 2000, warmup = 1000, seed = 123)
  tidy(fit, effects = "fixed") %>%
    mutate(
      p_value = 2 * (1 - pnorm(abs(estimate / std.error))),
      OR = exp(estimate),
      lower_CrI = exp(estimate - 1.96 * std.error),
      upper_CrI = exp(estimate + 1.96 * std.error)
    )
}

model_qft_035 <- run_stan_model("QFT_result_over_.35", "tb_exposition", no_prev)
model_qft_070 <- run_stan_model("QFT_result_over_.70", "tb_exposition", no_prev)
model_qft_400 <- run_stan_model("QFT_result_over_4.00", "tb_exposition", no_prev)

print(model_qft_035)
print(model_qft_070)
print(model_qft_400)

### GeneXpert models
run_gene_model <- function(outcome, data) {
  data <- data %>% mutate(GeneXpert.semiquant = relevel(factor(GeneXpert.semiquant), ref = "negative"))
  formula <- as.formula(paste0(outcome, " ~ GeneXpert.semiquant + prison + cell_density + age + months_prison + (1 | cell)"))
  fit <- stan_glmer(formula, family = binomial, data = data, 
                    chains = 4, iter = 2000, warmup = 1000, seed = 123)
  tidy(fit, effects = "fixed") %>%
    mutate(
      p_value = 2 * (1 - pnorm(abs(estimate / std.error))),
      OR = exp(estimate),
      lower_CrI = exp(estimate - 1.96 * std.error),
      upper_CrI = exp(estimate + 1.96 * std.error)
    )
}

model_gene_035 <- run_gene_model("QFT_result_over_.35", no_prev)
model_gene_070 <- run_gene_model("QFT_result_over_.70", no_prev)
model_gene_400 <- run_gene_model("QFT_result_over_4.00", no_prev)

print(model_gene_035)
print(model_gene_070)
print(model_gene_400)

# FIGURE 2: BIVARIATE GAMS =====================================================

# Duration of Incarceration vs QFT Positivity
prop_incarc <- no_prev %>%
  group_by(prison, months_prison) %>%
  summarise(QFT_proportion = mean(QFT_result_over_.35, na.rm = TRUE), .groups = "drop")

prison_labels <- c("MAXIMA" = "Prison A", "PED" = "Prison B", "PENAL" = "Prison C")

p1 <- ggplot(prop_incarc, aes(x = months_prison, y = QFT_proportion)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20), 
              method.args = list(family = "quasibinomial", gamma = 1),
              se = TRUE, color = "blue", fill = "blue", alpha = 0.2) +
  labs(x = "Duration of Incarceration (Months)", y = "QuantiFERON Positivity") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ prison, labeller = labeller(prison = prison_labels)) +
  theme_minimal(base_size = 12)

# Cell Density vs QFT Positivity
prop_density <- no_prev %>%
  group_by(prison, cell_density) %>%
  summarise(QFT_proportion = mean(QFT_result_over_.35, na.rm = TRUE), .groups = "drop")

p2 <- ggplot(prop_density, aes(x = cell_density, y = QFT_proportion)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 25),
              method.args = list(family = "quasibinomial", gamma = 5),
              se = TRUE, color = "blue", fill = "blue", alpha = 0.2) +
  labs(x = expression("Cell Density Index (# of people per m"^2*")"), y = "QuantiFERON Positivity") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ prison, labeller = labeller(prison = prison_labels)) +
  theme_minimal(base_size = 12)

# Save Panel
combined_plot <- arrangeGrob(p1, p2, nrow = 1)
combined_plot <- grobTree(
  combined_plot,
  textGrob("A", x = 0.02, y = 0.94, gp = gpar(fontsize = 18, fontface = "bold")),
  textGrob("B", x = 0.52, y = 0.94, gp = gpar(fontsize = 18, fontface = "bold"))
)

ggsave("~/Google Drive/My Drive/Subclinical TB/Results/combined_gams.png", 
       plot = combined_plot, width = 16, height = 4, dpi = 600)

# SAMPLE SIZE ==================================================================

var_cluster <- as.numeric(VarCorr(mod1.stan)$cell) # Extract random effect variance
icc <- var_cluster / (var_cluster + (pi^2 / 3)) # Compute ICC
m <- mean(table(no_prev$cell))
design_effect <- 1 + (m - 1) * icc # Design effect

