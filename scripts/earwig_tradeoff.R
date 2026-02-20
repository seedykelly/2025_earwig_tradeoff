
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
library(rptR)
library(modelsummary)
library(parameters)


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
# Raw phenotypic correlation: fecundity vs size
# =========================================

df.one <- df %>%
  filter(brood=="one")

## Pronotum

# Function to compute Pearson correlation
corr_raw_fun.pro <- function(data, indices) {
  d <- data[indices, ]
  cor(d$mean.pro, d$egg.number, use = "complete.obs")
}

# Run bootstrap (R = 5000)
set.seed(123)
boot_raw.pro <- boot(df.one, statistic = corr_raw_fun.pro, R = 5000)

# 95% CI using percentile method
ci_raw_95.pro <- boot.ci(boot_raw.pro, type = "perc", conf = 0.95)

# Extract point estimate and CI
r_raw.pro <- boot_raw.pro$t0
r_raw_ci95.pro <- ci_raw_95.pro$percent[4:5]  # 2.5% and 97.5% percentiles

cat("Raw phenotypic correlation (r_raw) = ", round(r_raw.pro, 3),
    ", 95% CI = [", round(r_raw_ci95.pro[1], 3), ", ", round(r_raw_ci95.pro[2], 3), "]\n")

## Body mass

# Function to compute Pearson correlation
corr_raw_fun.mass <- function(data, indices) {
  d <- data[indices, ]
  cor(d$mean.mass, d$egg.number, use = "complete.obs")
}

# Run bootstrap (R = 5000)
set.seed(123)
boot_raw.mass <- boot(df.one, statistic = corr_raw_fun.mass, R = 5000)

# 95% CI using percentile method
ci_raw_95.mass <- boot.ci(boot_raw.mass, type = "perc", conf = 0.95)

# Extract point estimate and CI
r_raw.mass <- boot_raw.mass$t0
r_raw_ci95.mass <- ci_raw_95.mass$percent[4:5]  # 2.5% and 97.5% percentiles

cat("Raw phenotypic correlation (r_raw) = ", round(r_raw.mass, 3),
    ", 95% CI = [", round(r_raw_ci95.mass[1], 3), ", ", round(r_raw_ci95.mass[2], 3), "]\n")

## PC1

# Function to compute Pearson correlation
corr_raw_fun.pc1 <- function(data, indices) {
  d <- data[indices, ]
  cor(d$pc1, d$egg.number, use = "complete.obs")
}

# Run bootstrap (R = 5000)
set.seed(123)
boot_raw.pc1 <- boot(df.one, statistic = corr_raw_fun.pc1, R = 5000)

# 95% CI using percentile method
ci_raw_95.pc1 <- boot.ci(boot_raw.pc1, type = "perc", conf = 0.95)

# Extract point estimate and CI
r_raw.pc1 <- boot_raw.pc1$t0
r_raw_ci95.pc1 <- ci_raw_95.pc1$percent[4:5]  # 2.5% and 97.5% percentiles

cat("Raw phenotypic correlation (r_raw) = ", round(r_raw.pc1, 3),
    ", 95% CI = [", round(r_raw_ci95.pc1[1], 3), ", ", round(r_raw_ci95.pc1[2], 3), "]\n")

## Scaled mass index (SCI)

# Function to compute Pearson correlation
corr_raw_fun.sci <- function(data, indices) {
  d <- data[indices, ]
  cor(d$mean.sci, d$egg.number, use = "complete.obs")
}

# Run bootstrap (R = 5000)
set.seed(123)
boot_raw.sci <- boot(df.one, statistic = corr_raw_fun.sci, R = 5000)

# 95% CI using percentile method
ci_raw_95.sci <- boot.ci(boot_raw.sci, type = "perc", conf = 0.95)

# Extract point estimate and CI
r_raw.sci <- boot_raw.sci$t0
r_raw_ci95.sci <- ci_raw_95.sci$percent[4:5]  # 2.5% and 97.5% percentiles

cat("Raw phenotypic correlation (r_raw) = ", round(r_raw.sci, 3),
    ", 95% CI = [", round(r_raw_ci95.sci[1], 3), ", ", round(r_raw_ci95.sci[2], 3), "]\n")

# =========================================
# Multivariate model
# =========================================
df <- df %>%
  mutate(
    egg_number_within_z  = scale(egg_number_within)[,1],
    egg_number_between_z    = scale(egg_number_between)[,1]
  )

df <- df %>%
  mutate(
    id = as.factor(id),
    brood = as.factor(brood)
  )

# Make "one" the reference level
df$brood <- relevel(df$brood, ref = "one")

levels(df$brood)

bf_number <- bf(
  egg.number ~ brood + (1 |p| id),
  family = negbinomial()
)
exp(-0.412)
bf_size <- bf(
  mean.egg.size ~ brood +
    egg_number_within_z +
    egg_number_between_z +
    (1 |p| id),
  family = gaussian()
)

mv_model <- bf_number + bf_size

get_prior(mv_model, data = df)

priors <- c(
  # -----------------------
  # FIXED EFFECTS
  # -----------------------
  prior(normal(0, 1), class = "b", resp = "eggnumber"),
  prior(normal(0, 1), class = "b", resp = "meaneggsize"),
  
  # -----------------------
  # INTERCEPTS
  # -----------------------
  prior(normal(0, 5), class = "Intercept", resp = "eggnumber"),
  prior(normal(3.5, 1), class = "Intercept", resp = "meaneggsize"),
  
  # -----------------------
  # RANDOM EFFECT SDs
  # -----------------------
  prior(exponential(1), class = "sd", group = "id", resp = "eggnumber"),
  prior(exponential(1), class = "sd", group = "id", resp = "meaneggsize"),
  
  # -----------------------
  # RESIDUAL SD (Gaussian only)
  # -----------------------
  prior(exponential(1), class = "sigma", resp = "meaneggsize"),
  
  # -----------------------
  # NEGATIVE BINOMIAL SHAPE
  # -----------------------
  prior(exponential(1), class = "shape", resp = "eggnumber"),
  
  # -----------------------
  # CORRELATION BETWEEN FEMALE INTERCEPTS
  # -----------------------
  prior(lkj(2), class = "cor", group = "id")
)

# fit_mv <- brm(
#   mv_model,
#   data = df,
#   prior = priors,
#   chains = 6,
#   cores = 4,
#   iter = 4000,
#   warmup = 1000,
#   control = list(adapt_delta = 0.95),save_pars = save_pars(all = TRUE),
#   file = "data/processed/fit_mv",
#   backend = "cmdstanr",
#   seed = 123
# )

fit_mv <- readRDS(file = "data/processed/fit_mv.rds")

print(summary(fit_mv), digits = 4)
summary(fit_mv)$fixed

## brood effect on egg number
fixef(fit_mv)["eggnumber_broodtwo", ]
exp(fixef(fit_mv)["eggnumber_broodtwo", "Estimate"])

## brood effect on egg size
fixef(fit_mv)["meaneggsize_broodtwo", ]

## Within-Female Trade-Off (Plasticity)
fixef(fit_mv)["meaneggsize_egg_number_within_z", ]

## Between-Female Effect (Strategic Allocation Axis)
fixef(fit_mv)["meaneggsize_egg_number_between_z", ]

## correlation between female intercepts
draws <- as_draws_df(fit_mv)
int.cor <- quantile(
     draws$cor_id__eggnumber_Intercept__meaneggsize_Intercept,
     c(0.025, 0.5, 0.975)
    )

# egg number repeatability (negative binomial; observation scale)

# Extract posterior draws
draws <- as_draws_df(fit_mv)
names(draws)[grep("egg_number", names(draws))]

mu_draws <- posterior_epred(fit_mv, resp = "eggnumber")
mu_bar <- apply(mu_draws, 1, mean)
# Among-female variance
VF <- draws$sd_id__eggnumber_Intercept^2

# Shape parameter
theta <- draws$shape_eggnumber
Vdist <- mu_bar + (mu_bar^2 / theta)
R_obs <- VF / (VF + Vdist)
mean_R <- mean(R_obs)
CI_R <- quantile(R_obs, probs = c(0.025, 0.975))

mean_R
CI_R

# =========================================
# Egg size variation
# =========================================

# earwig_egg_3 <- earwig_egg_2 %>%
#   mutate(
#    brood_factor = factor(brood,
#                           levels = c(1, 2),        # Original numeric values
#                           labels = c("one", "two"))) # New descriptive labels
# 
# egg.size.bf = bf(scale(mean.perim) ~ brood_factor + (1|a|id) + (0+brood_factor||id), sigma ~ brood_factor + (1|a|id)) # need to have (0+brood||id) in order to get sd_id_broodone and sd_id_broodtwo

# egg.size.model <- brm(egg.size.bf,
#                       data = earwig_egg_3,
#                       family = gaussian,
#                       cores = 6,
#                       chains = 6,
#                       warmup = 1000,
#                       iter = 4000,
#                       seed = 34, #make sure to set the seed to make results reproducible
#                       file = "data/processed/egg.size.model.rds",
#                       control = list(adapt_delta = 0.999))
# summary(egg.size.model)
# get_variables(egg.size.model)
# pp_check(egg.size.model,ndraws = 100)
# 
# post <- as_draws_df(egg.size.model)
# 
# ### between
# # Between-female SDs
# sd_int  <- post$sd_id__Intercept
# sd_b2   <- post$sd_id__brood_factortwo
# 
# # Convert to variance
# var_B1 <- sd_int^2
# var_B2 <- sd_int^2 + sd_b2^2
# 
# # Summaries
# summary_B1 <- quantile(var_B1, probs = c(.025, .5, .975))
# summary_B2 <- quantile(var_B2, probs = c(.025, .5, .975))
# 
# summary_B1
# summary_B2
# 
# diff_between <- var_B2 - var_B1
# 
# quantile(diff_between, probs = c(.025, .5, .975))
# mean(diff_between > 0)  # posterior probability brood2 > brood1

# # brood 1
# Median = 0.15
# 95% CI = [0.0006, 0.46]
# 
# # brood 2
# Median = 0.63
# 95% CI = [0.42, 1.02]
# Difference (B2 − B1):
# Median = 0.48
# 95% CI = [0.20, 0.87]
# Posterior probability > 0 = 1.00

# alternative model
 
# egg.size.bf.alt <- bf(
#   scale(mean.perim) ~ brood_factor + (0 + brood_factor | id),
#   sigma ~ brood_factor
# )
# 
# egg.size.model.alt <- brm(
#   egg.size.bf.alt,
#   data = earwig_egg_3,
#   family = gaussian(),
#   chains = 6,
#   cores = 6,
#   iter = 4000,
#   warmup = 1000,
#   file = "data/processed/egg.size.model2.rds",
#   control = list(adapt_delta = 0.999),
#   seed = 34
# )

egg.size.model.alt <- readRDS(file = "data/processed/egg.size.model2.rds")

# between-female
post_alt <- as_draws_df(egg.size.model.alt)

sd_B1_alt <- post_alt$sd_id__brood_factorone
sd_B2_alt <- post_alt$sd_id__brood_factortwo

var_B1_alt <- sd_B1_alt^2
var_B2_alt <- sd_B2_alt^2

quantile(var_B1_alt, c(.025,.5,.975))
quantile(var_B2_alt, c(.025,.5,.975))

diff_alt <- var_B2_alt - var_B1_alt
quantile(diff_alt, c(.025,.5,.975))
mean(diff_alt > 0)

# #Brood 1
# Median ≈ 0.77
# 95% CrI [0.52, 1.23]
# 
# #Brood 2
# Median ≈ 0.64
# 95% CrI [0.42, 1.00]
# 
# #Difference (B2 − B1)
# Median ≈ −0.14
# 95% CrI [−0.60, 0.29]
# Posterior P(B2 > B1) ≈ 0.26
# #Now brood 1 variance is slightly larger, and the difference is uncertain.

### within-clutch
# Linear predictor for sigma
sigma_int  <- post_alt$b_sigma_Intercept
sigma_b2   <- post_alt$b_sigma_brood_factortwo

# Convert from log scale
sigma_B1 <- exp(sigma_int)
sigma_B2 <- exp(sigma_int + sigma_b2)

# Convert to variance
var_within_B1 <- sigma_B1^2
var_within_B2 <- sigma_B2^2

# Summaries
quantile(var_within_B1, c(.025,.5,.975))
quantile(var_within_B2, c(.025,.5,.975))

diff_within <- var_within_B2 - var_within_B1

quantile(diff_within, c(.025,.5,.975))
mean(diff_within > 0)

# # brood 1
# Median = 0.19
# 95% CrI = [0.17, 0.22]
# 
# # brood 2
# Median = 0.32
# 95% CrI = [0.28, 0.37]
# 
# # Difference:
# Median = 0.13
# 95% CrI = [0.096, 0.164]
# Posterior probability > 0 = 1.00

# repeatability per brood
R_B1 <- var_B1_alt / (var_B1_alt + var_within_B1)
quantile(R_B1, c(.025,.5,.975))

R_B2 <- var_B2_alt / (var_B2_alt + var_within_B2)
quantile(R_B2, c(.025,.5,.975))

post_alt$cor_id__brood_factorone__brood_factortwo
quantile(post_alt$cor_id__brood_factorone__brood_factortwo, c(.025,.5,.975))

#loo(egg.size.model, egg.size.model.alt)

#==========================================


# =========================================
# Plots
# =========================================
# between-female correlation
# 
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
#   geom_smooth(method = "lm", se = TRUE, colour="black") +
#   labs(
#     x = "Mean egg number per female",
#     y = "Mean egg size per female"
#   ) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill="white"),
#         axis.title.x = element_text(size=14,face="bold"),
#         axis.title.y = element_text(size=14,face="bold"),
#         axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12))
# #ggsave(p_between, filename="Figure_1.jpg", width=10.83, height=10.83, dpi=300,antialias="default")
# 
# # # within-female correlation with reaction norms
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
# df_within <- df %>%
#   group_by(id) %>%
#   mutate(
#     egg_n_within = egg.number - mean(egg.number)
#   ) %>%
#   ungroup()
# 
# p_within <- ggplot(df_within,
#                    aes(x = egg_n_within,
#                        y = mean.egg.size,
#                        group = id)) +
#   geom_line(alpha = 0.3) +
#   geom_point(aes(colour=brood),alpha = 0.7, size = 2) +
# 
#   # Posterior ribbon (population-level)
#   geom_ribbon(
#     data = pred_within,
#     aes(x = egg_n_within, ymin = lwr, ymax = upr),
#     inherit.aes = FALSE,
#     alpha = 0.25
#   ) +
#   scale_color_discrete(name="Clutch order",labels = c("First", "Second")) +
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
#     y = "Mean egg size"
#   ) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill="white"),
#         axis.title.x = element_text(size=14,face="bold"),
#         axis.title.y = element_text(size=14,face="bold"),
#         axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12))
# #ggsave(p_within, filename="Figure_3.jpg", width=10.83, height=10.83, dpi=300,antialias="default")
# 
# # # residual correlation plot
# # 
# # # m_mv = your multivariate brms model with responses: mean.egg.size and egg.number
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
# p_resid <- ggplot(df_resid, aes(x = resid_egg_number, y = resid_egg_size)) +
#   geom_point(alpha = 0.5, size = 2) +  # residual points
#   geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +  # trend line
#   labs(
#     x = "Residual egg number",
#     y = "Residual egg size"
#   ) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill="white"),
#         axis.title.x = element_text(size=14,face="bold"),
#         axis.title.y = element_text(size=14,face="bold"),
#         axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12))
# #ggsave(p_resid, filename="figure_2.jpg", width=10.83, height=10.83, dpi=300,antialias="default")
# 
# post_alt <- as_draws_df(egg.size.model.alt)
# 
# df_plot <- tibble(
#   var_between_B1 = post_alt$sd_id__brood_factorone^2,
#   var_between_B2 = post_alt$sd_id__brood_factortwo^2,
#   var_within_B1  = var_within_B1,
#   var_within_B2  = var_within_B2
# ) %>%
#   pivot_longer(everything(),
#                names_to = "component",
#                values_to = "variance")
# 
# # Combine posterior draws
# df_between <- tibble(
#   brood = rep(c("first", "second"),
#               each = length(var_B1_alt)),
#   variance = c(var_B1_alt, var_B2_alt),
#   type = "Between-female egg size variation"
# )
# 
# df_within <- tibble(
#   brood = rep(c("first", "second"),
#               each = length(var_within_B1)),
#   variance = c(var_within_B1, var_within_B2),
#   type = "Within-clutch egg size variation"
# )
# 
# df_plot <- bind_rows(df_between, df_within)
# 
# p_variance <- ggplot(df_plot, aes(x = brood, y = variance)) +
#   
#   # posterior draws as semi-transparent dots
#   stat_dotsinterval(slab_fill="darkgrey", slab_color="darkgrey", point_interval = median_qi,
#                     .width = 0.95, quantiles=100) +
#   
#   # median + 95% credible interval
#   stat_summary(fun.data = median_qi,
#                fun.args = list(.width = 0.95),
#                geom = "pointrange",
#                color = "black",
#                size = 0.8) +
#   
#   facet_wrap(~ type, scales = "free_y") +
#   
#   labs(x = "Clutch laying order",
#        y = NULL) +
#   
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill="white"),
#         axis.title.x = element_text(size=14,face="bold"),
#         axis.title.y = element_text(size=14,face="bold"),
#         axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12))
# ggsave(p_variance, filename="figure_4.jpg", width=10.83, height=10.83, dpi=300,antialias="default")
# 


## ---- end



