
# rm(list=ls())

# ---- analysis ----

library(dplyr)
library(brms)
library(tidybayes)
library(ggplot2)
library(report)
library(posterior)
library(boot)
library(tibble)
library(broom.mixed)
library(broman)
library(tidyverse)
library(smatr)
library(factoextra)
library(officer)
library(flextable)
library(tidyverse)
library(rptR)

# =========================================
# Load data
# =========================================
earwig_body<-read.csv("data/raw/earwig_body_measurements.csv")
list_of_files <- list.files(path = "data/raw/egg_data_files",
                            recursive = TRUE,
                            pattern = "\\.csv$",
                            full.names = TRUE)
earwig_egg <- readr::read_csv(list_of_files)
earwig_egg$num <- as.factor(earwig_egg$num)
earwig_egg$id <- as.factor(earwig_egg$id)

# =========================================
# Trait repeatability
# =========================================
# body morphology
# make wide body data long
data_long.body <- gather(earwig_body, trait, measurement, pronotum_length_1:body_mass_3, factor_key=TRUE)
data_long.body$id <- as.factor(data_long.body$id)

# repeatability of pronotum length
# data_long.body.pronotum <- data_long.body %>%
#   filter(str_detect(trait, "pronotum_")) %>%
#   separate(trait, into = c("trait", "part2"), sep = "_(?=[^_]*_)", extra = "merge")
# 
# rep1 <- rpt(measurement ~ 1 +  (1| id), grname = c("id"),
#             data = data_long.body.pronotum, datatype = "Gaussian", nboot = 1000, npermut = 0)
# saveRDS(rep1, file = "data/processed/rep1.rds")
rep1 <- readRDS(file = "data/processed/rep1.rds")

# repeatability of body mass
# data_long.body.mass <- data_long.body %>%
#   filter(str_detect(trait, "body_")) %>%
#   separate(trait, into = c("trait", "part2"), sep = "_(?=[^_]*_)", extra = "merge")
# 
# rep2 <- rpt(measurement ~ 1 +  (1| id), grname = c("id"),
#             data = data_long.body.mass, datatype = "Gaussian", nboot = 1000, npermut = 0)
# saveRDS(rep2, file = "data/processed/rep2.rds")
rep2 <- readRDS(file = "data/processed/rep2.rds")

# calculate average pronotum length and body mass
earwig_body_2 <- earwig_body %>%
  rowwise() %>%
  mutate(mean.pronotum = mean(c_across(pronotum_length_1:pronotum_length_3)), mean.mass = mean(c_across(body_mass_1:body_mass_3))) %>%
  dplyr::select(-c(pronotum_length_1, pronotum_length_2, pronotum_length_3,body_mass_1, body_mass_2, body_mass_3)) %>%
  ungroup()

# =========================================
# Principal components
# =========================================
pc <- prcomp(earwig_body_2[,c(-1, -2)],
             center = TRUE,
             scale. = TRUE) # scale = TRUE causes PCA on correlation matrix
attributes(pc)
summary(pc)
fviz_pca_var(pc, col.var = "black")

pc1.variables <- get_pca_var(pc)
pc1 <- get_pca_ind(pc)
earwig_body_2$pc1 <- pc1$coord[,1]

res.var <- get_pca_var(pc)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# =========================================
# Scaled mass index
# =========================================
fit<-sma(log(earwig_body_2$mean.mass)~log(earwig_body_2$mean.pronotum))
summary(fit)
plot(fit)
plot(fit, which="qq")#assumptions are met

#reference population
L0<-mean(earwig_body_2$mean.pronotum)
L0
fit<-sma(log(earwig_body_2$mean.mass)~log(earwig_body_2$mean.pronotum))
fit
b.mass<-1.833134
Mi_hat<-(earwig_body_2$mean.mass)*(L0/earwig_body_2$mean.pronotum)^b.mass
earwig_body_2$Mi_hat <- Mi_hat

# =========================================
# Egg size repeatability
# =========================================

# egg perimeter
# make wide egg data long
# data_long <- gather(earwig_egg, perimeter, measurement, perim:perim3, factor_key=TRUE)
# data_long
# data_long$brood <- as.numeric(data_long$brood)
# data_long$num <- as.factor(data_long$num)
# data_long$id <- as.factor(data_long$id)

# analyze Rpt for each brood separately
# data_long_2 <- data_long %>%
#   unite(peregg, c("num", "id"))
# 
# data_long_3 <- data_long_2 %>%
#   filter(brood==1)
# 
# rep3 <- rpt(measurement ~ 1 +  (1| peregg), grname = c("peregg"),
#             data = data_long_3, datatype = "Gaussian", nboot = 1000, npermut = 0)
# saveRDS(rep3, file = "data/processed/rep3.rds")
rep3 <- readRDS(file = "data/processed/rep3.rds")

# data_long_4 <- data_long_2 %>%
#   filter(brood==2)
# 
# rep4 <- rpt(measurement ~ 1 +  (1| peregg), grname = c("peregg"),
#             data = data_long_4, datatype = "Gaussian", nboot = 1000, npermut = 0)
# saveRDS(rep4, file = "data/processed/rep4.rds")
rep4 <- readRDS(file = "data/processed/rep4.rds")

# calculate average egg perimeter
earwig_egg_2 <- earwig_egg %>%
  rowwise() %>%
  mutate(mean.perim = mean(c_across(perim:perim3))) %>%
  ungroup()

# =========================================
# Final dataset
# =========================================
# total.data <- earwig_egg_2 %>%
#   inner_join(earwig_body_2, by="id") %>%
#   mutate(brood = recode(brood, "1" = "one", "2" = "two")) %>%
#   ungroup() %>%
#   dplyr::select(-c(num.x, num.y)) %>%
#   group_by(id,brood) %>%
#   mutate(num=n()) %>%
#   group_by(id, brood) %>%
#   summarise(egg.number=n(), mean.egg.size = mean(mean.perim),pc1=mean(pc1), mean.pro=mean(mean.pronotum), mean.mass=mean(mean.mass), mean.sci=mean(Mi_hat)) 
# 
# df <- total.data %>%
#   group_by(id) %>%
#   mutate(
#     egg_number_between = mean(egg.number),
#     egg_number_within  = egg.number - egg_number_between
#   ) %>%
#   ungroup()

# saveRDS(df, file = "data/processed/df.rds")
df <- readRDS(file = "data/processed/df.rds")

# =========================================
# Summary stats
# =========================================

summary.stat <- df %>%
  group_by(brood) %>%
  summarize(mean.size = mean(mean.egg.size), 
            sd.size   = sd(mean.egg.size, na.rm = TRUE),
            n.size    = n(),
            se.size   = sd.size / sqrt(n.size),
            mean.num = mean(egg.number), 
            sd.num   = sd(egg.number, na.rm = TRUE),
            n.num    = n(),
            se.num   = sd.num / sqrt(n.num),)

# =========================================
# Effect of brood number
# =========================================
# priors_mv2 <- c(
#   
#   # Gaussian response: mean.egg.size
#   
#   # Intercept
#   prior(normal(0, 5), class = "Intercept", resp = "meaneggsize"),
#   
#   # Brood effect
#   prior(normal(0, 2), class = "b", resp = "meaneggsize"),
#   
#   # Random intercept SD
#   prior(student_t(3, 0, 2.5), class = "sd", resp = "meaneggsize"),
#   
#   # Residual SD
#   prior(student_t(3, 0, 2.5), class = "sigma", resp = "meaneggsize"),
#   
#   
#   # Poisson response: egg.number
#   # (log link scale!)
#   
#   # Intercept (on log scale)
#   prior(normal(2, 1.5), class = "Intercept", resp = "eggnumber"),
#   
#   # Brood effect (log scale)
#   prior(normal(0, 0.5), class = "b", resp = "eggnumber"),
#   
#   # Random intercept SD (log scale)
#   prior(student_t(3, 0, 1), class = "sd", resp = "eggnumber")
# )
# 
# m_mv_2 <- brm(data = df,
#               family = list(gaussian(), poisson()),
#               bf(mvbind(mean.egg.size, egg.number) ~ 1 + brood + (1|id)) + set_rescor(FALSE),
#               prior = priors_mv2,
#               file = "data/processed/m_mv_2",
#               backend = "cmdstanr",
#               chains=4,cores=4,iter = 4000, warmup= 1000)
# print(summary(m_mv_2), digits = 4)
m_mv_2 <- readRDS(file = "data/processed/m_mv_2.rds")
m_mv_2.sum = tidy(m_mv_2, effects = "fixed")

# =========================================
# Multivariate model for residual correlation
# =========================================
priors_mv <- c(
  
  # FIXED EFFECTS (egg size)
  prior(normal(0, 2), class = "b", resp = "meaneggsize"),
  
  # INTERCEPTS (both models)
  prior(normal(0, 5), class = "Intercept", resp = "meaneggsize"),
  prior(normal(0, 5), class = "Intercept", resp = "eggnumber"),
  
  # RANDOM INTERCEPT SDs
  prior(student_t(3, 0, 2.5), class = "sd", resp = "meaneggsize"),
  prior(student_t(3, 0, 2.5), class = "sd", resp = "eggnumber"),
  
  # RESIDUAL SDs
  prior(student_t(3, 0, 2.5), class = "sigma", resp = "meaneggsize"),
  prior(student_t(3, 0, 2.5), class = "sigma", resp = "eggnumber"),
  
  # RESIDUAL CORRELATION
  prior(lkj(2), class = "rescor")
)

f_size <- bf(
  mean.egg.size ~ egg_number_within + egg_number_between + (1 | id)
)

f_number <- bf(
  egg.number ~ (1 | id)
)

m_mv <- brm(
  f_size + f_number + set_rescor(TRUE),
  data = df,
  family = gaussian(),
  prior = priors_mv,
  file = "data/processed/m_mv",
  backend = "cmdstanr",
  control = list(adapt_delta = 0.999),
  chains = 4,
  cores = 4,
  iter = 4000
)

print(summary(m_mv), digits = 4)
m_mv <- readRDS(file = "data/processed/m_mv.rds")

# Extract posterior draws
post <- as_draws_df(m_mv)

# Extract residual correlation draws
rescor_draws <- post$rescor__meaneggsize__eggnumber

# Posterior mean
r_within <- mean(rescor_draws)

# 95% credible interval
within_credible <- quantile(rescor_draws, probs = c(0.025, 0.975))

# =========================================
# Raw phenotypic correlation (all broods)
# =========================================

# Function to compute Pearson correlation
corr_raw_fun <- function(data, indices) {
  d <- data[indices, ]
  cor(d$mean.egg.size, d$egg.number, use = "complete.obs")
}

# Run bootstrap (R = 5000)
set.seed(123)
boot_raw <- boot(df, statistic = corr_raw_fun, R = 5000)

# 95% CI using percentile method
ci_raw_95 <- boot.ci(boot_raw, type = "perc", conf = 0.95)

# Extract point estimate and CI
r_raw <- boot_raw$t0
r_raw_ci95 <- ci_raw_95$percent[4:5]  # 2.5% and 97.5% percentiles

cat("Raw phenotypic correlation (r_raw) = ", round(r_raw, 3),
    ", 95% CI = [", round(r_raw_ci95[1], 3), ", ", round(r_raw_ci95[2], 3), "]\n")

# =========================================
# Between-female correlation
# =========================================
# Compute female means
df_female <- df %>%
  group_by(id) %>%
  summarise(
    mean_egg_size  = mean(mean.egg.size, na.rm = TRUE),
    mean_egg_number = mean(egg.number, na.rm = TRUE),
    .groups = "drop"
  )

# Function for bootstrapping correlation across female means
corr_between_fun <- function(data, indices) {
  d <- data[indices, ]
  cor(d$mean_egg_size, d$mean_egg_number)
}

# Run bootstrap
set.seed(123)
boot_between <- boot(df_female, statistic = corr_between_fun, R = 5000)

# 95% CI using percentile method
ci_between_95 <- boot.ci(boot_between, type = "perc", conf = 0.95)

# Extract point estimate and CI
r_between <- boot_between$t0
r_between_ci95 <- ci_between_95$percent[4:5]

cat("Between-female correlation (r_between) = ", round(r_between, 3),
    ", 95% CI = [", round(r_between_ci95[1], 3), ", ", round(r_between_ci95[2], 3), "]\n")


# ## from multivariate model: (this model needs (1|p|id) as randome effect)
# VarCorr(m_mv)
# 
# # residual between-female correlation (correlation between intercepts)
# post <- as_draws_df(m_mv)
# r_between <- post$`cor_id__meaneggsize_Intercept__eggnumber_Intercept`
# posterior_summary(r_between)
# mean(r_between)
# quantile(r_between, c(.025, .975))
# 
# # total between-female correlation
# post <- as_draws_df(m_mv)
# beta_B  <- post$b_meaneggsize_egg_number_between
# var_size_id   <- post$sd_id__meaneggsize_Intercept^2
# var_number_id <- post$sd_id__eggnumber_Intercept^2
# cor_id <- post$cor_id__meaneggsize_Intercept__eggnumber_Intercept
# cov_id <- cor_id * sqrt(var_size_id * var_number_id)
# 
# # Between covariance
# cov_between <- beta_B * var_number_id + cov_id
# 
# # Between variance in egg size
# var_size_between <- beta_B^2 * var_number_id +
#   var_size_id +
#   2 * beta_B * cov_id
# 
# # Total between correlation
# r_between_total <- cov_between /
#   sqrt(var_size_between * var_number_id)
# mean(r_between_total)
# quantile(r_between_total, c(.025, .975))

# =========================================
# Fixed-effects model
# =========================================

# priors <- c(
#   prior(normal(0, 0.01), class = "b", coef = "egg_number_within"),
#   prior(normal(0, 0.01), class = "b", coef = "egg_number_between"),
#   prior(normal(0, 1), class = "Intercept"),
#   prior(exponential(1), class = "sd"),
#   prior(exponential(1), class = "sigma")
# )
# 
# fit.null <- brm(
#   mean.egg.size ~
#     egg_number_within +
#     egg_number_between +
#     (1 | id),
#   data = df,
#   family = gaussian(),
#   prior = priors,
#   save_pars = save_pars(all = TRUE),
#   file = "data/processed/fit.null7",
#   backend = "cmdstanr",
#   chains = 4, cores = 4, thin = 6, iter = 4000
# )
# print(summary(fit.null), digits = 4)
# fit.null <- add_criterion(fit.null, "loo", moment_match=TRUE) # best model
fit.null <- readRDS(file = "data/processed/fit.null.rds")

# fit.pro <- brm(
#   mean.egg.size ~
#     egg_number_within +
#     egg_number_between + scale(mean.pro) +
#     (1 | id),
#   data = df,
#   control = list(adapt_delta = 0.999),
#   save_pars = save_pars(all = TRUE),
#   family = gaussian(),
#   file = "data/processed/fit.pro",
#   backend = "cmdstanr",
#   chains = 4, cores = 4,iter = 4000
# )
# print(summary(fit.pro), digits = 4)
# fit.pro <- add_criterion(fit.pro, "loo", moment_match=TRUE)
# 
# fit.mass <- brm(
#   mean.egg.size ~
#     egg_number_within +
#     egg_number_between + scale(mean.mass) +
#     (1 | id),
#   data = df,
#   control = list(adapt_delta = 0.999),
#   save_pars = save_pars(all = TRUE),
#   family = gaussian(),
#   file = "data/processed/fit.mass",
#   backend = "cmdstanr",
#   chains = 4, cores = 4,iter = 4000
# )
# print(summary(fit.mass), digits = 4)
# fit.mass <- add_criterion(fit.mass, "loo", moment_match=TRUE)
# 
# fit.pc1 <- brm(
#   mean.egg.size ~
#     egg_number_within +
#     egg_number_between + scale(pc1) +
#     (1 | id),
#   data = df,
#   control = list(adapt_delta = 0.999),
#   save_pars = save_pars(all = TRUE),
#   family = gaussian(),
#   file = "data/processed/fit.pc1",
#   backend = "cmdstanr",
#   chains = 4, cores = 4,iter = 4000
# )
# print(summary(fit.pc1), digits = 4)
# fit.pc1 <- add_criterion(fit.pc1, "loo", moment_match=TRUE)
# 
# fit.smi <- brm(
#   mean.egg.size ~
#     egg_number_within +
#     egg_number_between + scale(mean.sci) +
#     (1 | id),
#   data = df,
#   control = list(adapt_delta = 0.999),
#   save_pars = save_pars(all = TRUE),
#   family = gaussian(),
#   file = "data/processed/fit.smi",
#   backend = "cmdstanr",
#   chains = 4, cores = 4,iter = 4000
# )
# print(summary(fit.smi), digits = 4)
# fit.smi <- add_criterion(fit.smi, "loo", moment_match=TRUE)
# 
# model.comparison <- loo_compare(fit.null, fit.pro, fit.mass, fit.pc1, fit.smi, criterion = "loo")
# saveRDS(model.comparison, file = "data/processed/model.comparison.rds")
model.comparison <- readRDS(file = "data/processed/model.comparison.rds")

# model.comparison.report <-report(model.comparison) # fit.null is best model
# saveRDS(model.comparison.report, file = "data/processed/model.comparison.report.rds")
model.comparison.report <- readRDS(file = "data/processed/model.comparison.report.rds")

ran_int <- posterior_summary(fit.null, variable = "sd_id__Intercept")

# =========================================
# Standardized effects
# =========================================

post <- as_draws_df(fit.null)
sd_y <- sd(df$mean.egg.size)

sd_within  <- sd(df$egg_number_within)
sd_between <- sd(df$egg_number_between)

beta_within_std  <- post$b_egg_number_within  * sd_within  / sd_y
beta_between_std <- post$b_egg_number_between * sd_between / sd_y

std_effects_table <- data.frame(
  Predictor = c("Egg number (within female)",
                "Egg number (between females)"),
  Estimate = c(mean(beta_within_std),
               mean(beta_between_std)),
  LCI_95 = c(quantile(beta_within_std,  0.025),
             quantile(beta_between_std, 0.025)),
  UCI_95 = c(quantile(beta_within_std,  0.975),
             quantile(beta_between_std, 0.975))
)
# saveRDS(std_effects_table, file = "data/processed/std_effects_table.rds")

# =========================================
# Brood adjusted slope
# =========================================

# m_wb_brood <- brm(
#   mean.egg.size ~ 
#     egg_number_within +
#     egg_number_between +
#     brood +
#     (1 | id),
#   data = df,
#   family = gaussian(),
#   backend = "cmdstanr",
#   chains = 4, cores = 4,
#   file = "data/processed/m_wb_brood",
#   iter = 4000, warmup = 1000
# )
# print(summary(m_wb_brood), digits = 4)

m_wb_brood <- readRDS(file = "data/processed/m_wb_brood.rds")

clutch_order_within <- brms::fixef(m_wb_brood)["egg_number_within", ]

# =========================================
# Random slopes model
# =========================================
# m_rs <- brm(
#   mean.egg.size ~ 
#     egg_number_within + 
#     egg_number_between +
#     (1 + egg_number_within | id),
#   data = df,
#   family = gaussian(),
#   prior = c(
#     prior(normal(0, 0.5), class = "b"),
#     prior(normal(0, 1), class = "Intercept"),
#     prior(exponential(1), class = "sd"),
#     prior(lkj(2), class = "cor")
#   ),
#   chains = 4,
#   cores = 4,
#   iter = 5000,
#   control = list(adapt_delta = 0.99),
#   file = "data/processed/m_rs",
#   
# )
# print(summary(m_rs), digits = 4)

m_rs <- readRDS(file = "data/processed/m_rs.rds")

ran_int.slopes <- posterior_summary(m_rs, variable = "sd_id__egg_number_within")
ran_int.slopes.cor <- posterior_summary(m_rs, variable = "cor_id__Intercept__egg_number_within")


# # Convert brms model to draws_df
# post <- as_draws_df(m_rs)
# 
# # Keep only random effects
# post_re <- post %>%
#   select(starts_with("r_id["))
# 
# # Add a draw ID
# post_re <- post_re %>%
#   mutate(draw = row_number())
# post_long <- post_re %>%
#   pivot_longer(
#     cols = -draw,
#     names_to = "variable",
#     values_to = "value"
#   ) %>%
#   # Extract female ID and effect
#   mutate(
#     id = sub("r_id\\[([^,]+),.*\\]", "\\1", variable),
#     effect = sub("r_id\\[[^,]+,(.*)\\]", "\\1", variable)
#   )
# post_wide <- post_long %>%
#   pivot_wider(
#     id_cols = c(draw, id),   # include draw so each row is unique
#     names_from = effect,
#     values_from = value
#   )
# 
# ggplot(post_wide, aes(x = egg_number_between, y = Intercept)) +
#   geom_point(alpha = 0.2) +
#   labs(
#     x = "Female-specific slope (egg_number_between)",
#     y = "Female-specific intercept (mean egg size)"
#   )

# # Posterior median + 95% credible intervals per female
# post_summary <- post_wide %>%
#   group_by(id) %>%
#   summarise(
#     Intercept = median(Intercept),
#     Slope = median(egg_number_between),
#     Slope_LCI = quantile(egg_number_between, 0.025),
#     Slope_UCI = quantile(egg_number_between, 0.975)
#   )
# 
# ggplot(post_summary, aes(x = Slope, y = Intercept)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin = Slope_LCI, xmax = Slope_UCI), alpha = 0.5) +
#   labs(
#     x = "Female-specific slope (egg_number_between)",
#     y = "Female-specific intercept (mean egg size)"
#   ) +
#   theme_classic()

## ---- end

# =========================================
# 3-panel plot
# =========================================
# between-female correlation

# df_between <- df %>%
#   group_by(id) %>%
#   summarise(
#     mean_egg_number = mean(egg.number, na.rm = TRUE),
#     mean_egg_size   = mean(mean.egg.size, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# p_between <- ggplot(df_between, aes(x = mean_egg_number, y = mean_egg_size)) +
#   geom_point(size = 2) +
#   geom_smooth(method = "lm", se = TRUE) +
#   labs(
#     x = "Mean egg number per female",
#     y = "Mean egg size per female",
#     title = "Between-female correlation"
#   ) +
#   theme_classic()
# ggsave(p_between, filename="p_between.jpg", width=10.83, height=10.83, dpi=300,antialias="default")    
# 
# # within-female correlation with reaction norms
# 
# newdata_within <- tibble(
#   egg_number_within  = seq(
#     min(df$egg_number_within, na.rm = TRUE),
#     max(df$egg_number_within, na.rm = TRUE),
#     length.out = 100
#   ),
#   egg_number_between = 0
# )
# 
# pred_within <- m_rs %>%
#   add_epred_draws(
#     newdata = newdata_within,
#     re_formula = NA   # population-level only
#   ) %>%
#   group_by(egg_number_within) %>%
#   summarise(
#     mean = mean(.epred),
#     lwr  = quantile(.epred, 0.025),
#     upr  = quantile(.epred, 0.975),
#     .groups = "drop"
#   ) %>%
#   rename(egg_n_within = egg_number_within)
# 
# p_within <- ggplot(df_within,
#                    aes(x = egg_n_within,
#                        y = mean.egg.size,
#                        group = id)) +
#   geom_line(alpha = 0.3) +
#   geom_point(alpha = 0.4, size = 1) +
#   
#   # Posterior ribbon (population-level)
#   geom_ribbon(
#     data = pred_within,
#     aes(x = egg_n_within, ymin = lwr, ymax = upr),
#     inherit.aes = FALSE,
#     alpha = 0.25
#   ) +
#   
#   # Posterior mean line (population-level)
#   geom_line(
#     data = pred_within,
#     aes(x = egg_n_within, y = mean),
#     inherit.aes = FALSE,
#     linewidth = 1
#   ) +
#   
#   labs(
#     x = "Within-female deviation in egg number",
#     y = "Mean egg size",
#     title = "B) Within-female plasticity (posterior mean)"
#   ) +
#   theme_classic(base_size = 12)
# ggsave(p_within, filename="p_within.jpg", width=10.83, height=10.83, dpi=300,antialias="default")    
# 
# # residual correlation plot
# 
# # m_mv = your multivariate brms model with responses: mean.egg.size and egg.number
# 
# # Get posterior residuals (array: draws x observations x response)
# resid_array <- residuals(m_mv, summary = FALSE)  # do NOT summarize yet
# 
# # Compute posterior mean residual for each observation
# resid_mean <- apply(resid_array, c(2, 3), mean)  # dims: obs x response
# 
# # Convert to data frame, combine with original df
# df_resid <- df %>%
#   mutate(
#     resid_egg_size   = as.numeric(resid_mean[, "meaneggsize"]),
#     resid_egg_number = as.numeric(resid_mean[, "eggnumber"])
#   )
# 
# # Posterior draws of the residual correlation
# rescor_draws <- posterior_samples(m_mv, pars = "rescor__meaneggsize__eggnumber")
# 
# # Extract numeric vector
# rescor_vec <- rescor_draws[[1]]
# 
# # Summarize posterior: mean + 95% credible interval
# rescor_summary <- c(
#   mean = mean(rescor_vec),
#   lwr  = quantile(rescor_vec, 0.025),
#   upr  = quantile(rescor_vec, 0.975)
# )
# 
# rescor_summary
# 
# p_resid <- ggplot(df_resid, aes(x = resid_egg_number, y = resid_egg_size)) +
#   geom_point(alpha = 0.5, size = 2) +  # residual points
#   geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +  # trend line
#   labs(
#     x = "Residual egg number",
#     y = "Residual egg size",
#     title = paste0(
#       "Residual correlation (r = ",
#       round(rescor_summary["mean"], 2),
#       ", 95% CrI [",
#       round(rescor_summary["lwr.2.5%"], 2), ", ",
#       round(rescor_summary["upr.97.5%"], 2), "])"
#     )
#   ) +
#   theme_classic(base_size = 12)
# ggsave(p_resid, filename="p_resid.jpg", width=10.83, height=10.83, dpi=300,antialias="default")    
# 





