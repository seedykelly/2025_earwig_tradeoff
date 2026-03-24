











tapply(earwig_egg_summary.2$num, earwig_egg_summary.2$id, sd)
sd(earwig_egg_summary.2$egg.perim)
sd(earwig_egg_summary.2$egg_number_within)
sd(earwig_egg_summary.2$egg_number_between)


earwig_egg_summary.2 %>%
  group_by(id) %>%
  summarise(
    mean_eggs = mean(num),
    sd_eggs = sd(num)
  )

# Extract posterior draws for slope and intercept
draws <- as_draws_df(fit.earwig)

# Create grid of x values
xseq <- seq(min(earwig_egg_summary.2$egg_number_within),
            max(earwig_egg_summary.2$egg_number_within),
            length.out = 100)

# Posterior predictions
pred_grid <- expand.grid(
  egg_number_within = xseq,
  egg_number_between = mean(earwig_egg_summary.2$egg_number_between),
  id = NA
)

fitted_vals <- fitted(fit.earwig,
                      newdata = pred_grid,
                      re_formula = NA,
                      summary = TRUE)

pred_grid$mean <- fitted_vals[, "Estimate"]
pred_grid$lwr  <- fitted_vals[, "Q2.5"]
pred_grid$upr  <- fitted_vals[, "Q97.5"]

# Plot
ggplot(earwig_egg_summary.2, aes(x = egg_number_within, y = egg.perim)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = pred_grid, aes(y=mean,ymin = lwr, ymax = upr),
              fill = "steelblue", alpha = 0.3) +
  geom_line(data = pred_grid, aes(y = mean), size = 1.2, color = "black") +
  theme_classic() +
  labs(
    x = "Clutch size deviation (within female)",
    y = "Mean egg size",
    title = "Within-Female Egg Size–Clutch Size Trade-Off"
  )

posterior_summary(fit, probs = c(0.025, 0.975))

#Relationship between within-female clutch size deviation and mean egg size. Lines represent 
#posterior mean (solid) and 95% credible interval (shaded). The near-zero slope indicates negligible within-female trade-off.




# fit <- brm(data = dat, 
#           family = list(gaussian(), gaussian()),
#           bf(mvbind(mean_egg_size, egg_number) ~ 1 + (1|female)) + set_rescor(TRUE),
#           backend = "cmdstanr",
#           chains=4,cores=4,iter = 4000, warmup= 1000)
# summary(fit)$rescor_pars

###
#####
library(brms)
library(dplyr)
library(ggplot2)

set.seed(123)

# Design: 42 females, 2 clutches each
n_females <- 42
n_clutches <- 2
N <- n_females * n_clutches

female <- factor(rep(1:n_females, each = n_clutches))

# Between-female variation in egg size
female_intercept <- rnorm(n_females, 0, 0.1)
female_intercept_rep <- rep(female_intercept, each = n_clutches)

# Clutch-level egg number (Poisson-like)
egg_number <- rpois(N, lambda = 8 + female_intercept_rep*2)

# Within-female small negative trade-off
within_effect <- -0.01  # very small effect
mean_egg_size <- 3.67 +
  female_intercept_rep +
  within_effect * (egg_number - ave(egg_number, female)) +
  rnorm(N, 0, 0.05)  # residual variation

dat <- data.frame(
  female = female,
  clutch = rep(1:2, times = n_females),
  egg_number = egg_number,
  mean_egg_size = mean_egg_size
)

# Create within- and between-female components
dat <- dat %>%
  group_by(female) %>%
  mutate(
    egg_number_between = mean(egg_number),
    egg_number_within  = egg_number - egg_number_between
  ) %>%
  ungroup()

fit <- brm(
  mean_egg_size ~ egg_number_within + egg_number_between + (1|female),
  data = dat,
  family = gaussian(),
  backend = "cmdstanr",
  chains = 4, cores = 4,
  iter = 4000
)

summary(fit)

library(tidybayes)

# Create grid for plotting
xseq <- seq(min(dat$egg_number_within), max(dat$egg_number_within), length.out = 100)

# Posterior predictions
pred_grid <- expand.grid(
  egg_number_within = xseq,
  egg_number_between = mean(dat$egg_number_between),
  female = NA
)

fitted_vals <- fitted(fit,
                      newdata = pred_grid,
                      re_formula = NA,
                      summary = TRUE)

pred_grid$mean <- fitted_vals[, "Estimate"]
pred_grid$lwr  <- fitted_vals[, "Q2.5"]
pred_grid$upr  <- fitted_vals[, "Q97.5"]

# Plot
ggplot(dat, aes(x = egg_number_within, y = mean_egg_size)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = pred_grid, aes(y=mean,ymin = lwr, ymax = upr),
              fill = "steelblue", alpha = 0.3) +
  geom_line(data = pred_grid, aes(y = mean), size = 1.2, color = "black") +
  theme_classic() +
  labs(
    x = "Clutch size deviation (within female)",
    y = "Mean egg size",
    title = "Within-Female Egg Size–Clutch Size Trade-Off"
  )

posterior_summary(fit, probs = c(0.025, 0.975))

######
library(dplyr)
library(brms)

df <- dat %>%
  group_by(id) %>%
  mutate(
    egg_n_mean = mean(egg.number),
    egg_n_within = egg.number - egg_n_mean
  ) %>%
  ungroup()


m_tradeoff <- brm(
  mean.egg.size ~ 
    egg_n_within +        # within-female (trade-off)
    egg_n_mean +          # between-female differences
    (1 | id),             # female identity
  data = df,
  backend = "cmdstanr",
  family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class = "b"),
    prior(normal(0, 1), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  chains = 4,
  cores = 4,
  iter = 4000,
  control = list(adapt_delta = 0.95)
)
print(summary(m_tradeoff), digits = 4)

library(brms)

# random slopes model


library(dplyr)
library(knitr)


# Raw values
effects_raw <- tibble(
  Effect = c("Within-female slope", "Between-female slope", "Random-slope SD"),
  Estimate_raw = c(0.0035, -0.0041, 0.00463),
  SD_X = c(13.4981, 12.984, 13.4981),
  SD_Y = 0.2177  # observed SD of mean.egg.size
)

# Standardize
effects_std <- effects_raw %>%
  mutate(
    Estimate_std = Estimate_raw * (SD_X / SD_Y)
  )

# Format table
kable(effects_std,
      col.names = c("Effect", "Raw Estimate", "SD Predictor", "SD Outcome", "Standardized Estimate"),
      digits = 3)

# Approximate credible intervals for illustration
effects_plot <- tibble(
  Effect = c("Within-female slope", "Between-female slope", "Random-slope SD"),
  Raw = c(0.0035, -0.0041, 0.00463),
  Raw_L95 = c(-0.0005, -0.0081, 0.00027),
  Raw_U95 = c(0.0077, -0.0002, 0.0111),
  Std = c(0.254, -0.284, 0.334),
  Std_L95 = c(-0.036, -0.486, 0.001),
  Std_U95 = c(0.526, -0.012, 0.476)
)

# Plot
ggplot(effects_plot, aes(x=Effect)) +
  geom_point(aes(y=Raw, color="Raw"), size=3) +
  geom_errorbar(aes(ymin=Raw_L95, ymax=Raw_U95, color="Raw"), width=0.2) +
  geom_point(aes(y=Std, color="Standardized"), size=3, shape=17) +
  geom_errorbar(aes(ymin=Std_L95, ymax=Std_U95, color="Standardized"), width=0.2) +
  scale_color_manual(values=c("Raw"="blue", "Standardized"="red")) +
  labs(y="Effect size", color="Scale") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, hjust=1))

# Female means (between-female)
df_between <- df %>%
  group_by(id) %>%
  summarise(
    egg_n_mean = mean(egg.number),
    egg_size_mean = mean(mean.egg.size),
    .groups = "drop"
  )

# Within-female centered data
df_within <- df %>%
  group_by(id) %>%
  mutate(
    egg_n_within = egg.number - mean(egg.number)
  ) %>%
  ungroup()

p_between <- ggplot(df_between,
                    aes(x = egg_n_mean, y = egg_size_mean)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(
    x = "Mean egg number per female",
    y = "Mean egg size per female",
    title = "A) Between-female relationship"
  ) +
  theme_classic(base_size = 12)

p_within <- ggplot(df_within,
                   aes(x = egg_n_within,
                       y = mean.egg.size,
                       group = id)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(aes(group = 1),
              method = "lm",
              se = TRUE,
              color = "black",
              linewidth = 1) +
  labs(
    x = "Within-female deviation in egg number",
    y = "Mean egg size",
    title = "B) Within-female plasticity"
  ) +
  theme_classic(base_size = 12)

new_between <- tibble(
  egg_n_mean   = seq(min(df$egg_n_mean),
                     max(df$egg_n_mean),
                     length.out = 100),
  egg_n_within = 0
)

epred_between <- posterior_epred(
  m_rs,
  newdata = new_between,
  re_formula = NA
)

pred_between <- new_between %>%
  mutate(
    mean = colMeans(epred_between),
    lwr  = apply(epred_between, 2, quantile, 0.025),
    upr  = apply(epred_between, 2, quantile, 0.975)
  )

library(ggplot2)

p_between <- ggplot(df_between,
                    aes(x = egg_n_mean, y = egg_size_mean)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_ribbon(data = pred_between,
              aes(y= mean,x=egg_n_mean,ymin = lwr, ymax = upr),
              inherit.aes = FALSE,
              alpha = 0.25) +
  geom_line(data = pred_between,
            aes(x = egg_n_mean, y = mean),
            linewidth = 1) +
  labs(
    x = "Mean egg number per female",
    y = "Mean egg size per female",
    title = "A) Between-female relationship"
  ) +
  theme_classic(base_size = 12)




## residual correlation

m_mv <- brm(data = dat, 
            family = list(gaussian(), gaussian()),
            bf(mvbind(mean.egg.size, egg.number) ~ 1 + (1|id)) + set_rescor(TRUE),
            backend = "cmdstanr",
            chains=4,cores=4,iter = 4000, warmup= 1000)
print(summary(m_mv), digits = 4)
get_variables(m_mv)



# Load required packages
library(dplyr)

# fixed effect of brood number

