
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
library(cowplot)

sci_latex <- function(x, digits = 2) {
  if (length(x) != 1 || is.na(x)) return(NA_character_)
  
  s <- as.character(formatC(x, format = "e", digits = digits))
  parts <- strsplit(s, "e")[[1]]
  
  mantissa <- parts[1]
  exponent <- as.integer(parts[2])
  
  paste0(mantissa, " \\times 10^{", exponent, "}")
}

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
#write.csv(df, file = "earwig_data.csv", row.names = FALSE)

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

fit_mv <- brm(
  mv_model,
  data = df,
  prior = priors,
  chains = 6,
  cores = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.95),save_pars = save_pars(all = TRUE),
  file = "data/processed/fit_mv1",
  backend = "cmdstanr",
  seed = 123
)

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

# EGG NUMBER repeatability (negative binomial; observation scale)

draws <- as_draws_df(fit_mv)

sigma2 <- draws$sd_id__eggnumber_Intercept^2
theta  <- draws$shape_eggnumber

mu_draws <- posterior_epred(fit_mv, resp = "eggnumber")

# Preallocate
VF_obs <- numeric(nrow(mu_draws))
Vdist  <- numeric(nrow(mu_draws))

for (i in seq_len(nrow(mu_draws))) {
  mu_i <- mu_draws[i, ]
  
  # Among-female variance (observation scale)
  VF_obs[i] <- mean(mu_i^2) * (exp(sigma2[i]) - 1)
  
  # Distribution variance
  Vdist[i] <- mean(mu_i + (mu_i^2 / theta[i]))
}

R_obs <- VF_obs / (VF_obs + Vdist)

mean_R <- mean(R_obs)
CI_R <- quantile(R_obs, c(0.025, 0.975))

mean_R
CI_R

# =========================================
# Egg size variation
# =========================================

# alternative model

priors <- c(
  # Mean model
  prior(normal(0, 0.5), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  
  # Random slopes
  prior(exponential(2), class = "sd"),
  prior(lkj(2), class = "cor"),
  
  # Residual model
  prior(normal(0, 0.5), class = "Intercept", dpar = "sigma"),
  prior(normal(0, 0.5), class = "b", dpar = "sigma")
)
 
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
#   prior = priors,
#   iter = 4000,
#   warmup = 1000,
#   file = "data/processed/egg.size.model.rds",
#   control = list(adapt_delta = 0.999),
#   seed = 34
# )

egg.size.model.alt <- readRDS(file = "data/processed/egg.size.model.rds")

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

df <- df %>%
  group_by(id) %>%
  mutate(
    mean_eggnumber = mean(egg.number, na.rm = TRUE),
    egg_number_within = egg.number - mean_eggnumber
  ) %>%
  ungroup() %>%
  mutate(
    egg_number_between = mean_eggnumber,
    egg_number_within_z  = as.numeric(scale(egg_number_within)),
    egg_number_between_z = as.numeric(scale(egg_number_between))
  )

new_within <- tibble(
  egg_number_within_z = seq(min(df$egg_number_within_z),
                            max(df$egg_number_within_z),
                            length.out = 100),
  egg_number_between_z = 0,   # hold between constant
  brood = "one"             # or reference level
)

new_between <- tibble(
  egg_number_between_z = seq(min(df$egg_number_between_z),
                             max(df$egg_number_between_z),
                             length.out = 100),
  egg_number_within_z = 0,    # hold within constant
  brood = "one"
)

female_means <- df %>%
  group_by(id) %>%
  summarise(
    mean_eggnumber = mean(egg.number, na.rm = TRUE),
    mean_eggsize   = mean(mean.egg.size, na.rm = TRUE),
    .groups = "drop"
  )

# Within predictions
pred_within <- add_epred_draws(
  fit_mv,
  newdata = new_within,
  resp = "meaneggsize",
  re_formula = NA
) %>%
  mutate(type = "Within-female")

# Between predictions
pred_between <- add_epred_draws(
  fit_mv,
  newdata = new_between,
  resp = "meaneggsize",
  re_formula = NA
) %>%
  mutate(type = "Between-female")

# Combine
pred_all <- bind_rows(pred_within, pred_between)

# PANEL A — Within-female
p_within <- ggplot() +

  # Raw reaction norms
  geom_line(data = df,
            aes(x = egg_number_within_z,
                y = mean.egg.size,
                group = id),
            alpha = 0.3) +

  geom_point(data = df,
             aes(x = egg_number_within_z,
                 y = mean.egg.size, shape=brood),
             size = 3, colour="blue") +

    # Bayesian slope only
  stat_summary(
    data =
      pred_within,
    aes(x = egg_number_within_z,
        y = .epred,
        group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1.3,
    color = "black"
  ) +
  scale_shape_manual(name="Clutch order",labels = c("First", "Second"),values = c(1, 16)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.2, 0.8),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  labs(
    x = "Within-female deviation in egg number",
    y = "Mean egg size (mm)"
  )

# PANEL B — Between-female
p_between <- ggplot() +

  # Female means (raw intuition)
  geom_point(data = female_means,
             aes(x = scale(mean_eggnumber),
                 y = mean_eggsize), colour = "red",
             size = 3) +

  # Bayesian slope only
  stat_summary(
    data = pred_between,
    aes(x = egg_number_between_z,
        y = .epred,
        group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1.3,
    color = "black"
  ) +

  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  labs(
    x = "Between-female deviation in egg number",
    y = "Mean egg size (mm)"
  )


# Combine
figure_1 <- plot_grid(p_within, p_between, ncol=2, nrow=1,
                    labels = c('(a)', '(b)'))
ggsave(figure_1, filename="figure_1.jpg", width=10.83, height=10.83, dpi=300,antialias="default")


# figure 2
egg_seq <- seq(-2, 2, length.out = 100)

# # Within-female grid
new_within <- data.frame(
  brood = "one",
  egg_number_within_z = egg_seq,
  egg_number_between_z = 0
)

# Between-female grid
new_between <- data.frame(
  brood = "one",
  egg_number_within_z = 0,
  egg_number_between_z = egg_seq
)
pred_within <- add_epred_draws(
  fit_mv,
  newdata = new_within,
  resp = "meaneggsize",
  re_formula = NA
)

post_alt <- as_draws_df(egg.size.model.alt)

df_plot <- tibble(
  var_between_B1 = post_alt$sd_id__brood_factorone^2,
  var_between_B2 = post_alt$sd_id__brood_factortwo^2,
  var_within_B1  = var_within_B1,
  var_within_B2  = var_within_B2
) %>%
  pivot_longer(everything(),
               names_to = "component",
               values_to = "variance")

# Combine posterior draws
df_between <- tibble(
  brood = rep(c("first", "second"),
              each = length(var_B1_alt)),
  variance = c(var_B1_alt, var_B2_alt),
  type = "Between-female egg size variation"
)

df_within <- tibble(
  brood = rep(c("first", "second"),
              each = length(var_within_B1)),
  variance = c(var_within_B1, var_within_B2),
  type = "Within-clutch egg size variation"
)

df_plot <- bind_rows(df_between, df_within)

p_variance <- ggplot(df_plot, aes(x = brood, y = variance)) +

  # posterior draws as semi-transparent dots
  stat_dotsinterval(slab_fill="orange", slab_color="orange", point_interval = median_qi,
                    .width = 0.95, quantiles=100) +

  # median + 95% credible interval
  stat_summary(fun.data = median_qi,
               fun.args = list(.width = 0.95),
               geom = "pointrange",
               color = "black",
               size = 0.8) +

  facet_wrap(~ type, scales = "free_y") +

  labs(x = "Clutch laying order",
       y = "Variance in egg size (mm)") +

  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 14),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))
ggsave(p_variance, filename="figure_2.jpg", width=10.83, height=10.83, dpi=300,antialias="default")



## ---- end

