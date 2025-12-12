#rm(list=ls())

# ---- load_libraries ----
library(tidyverse); library(readr)
library(dplyr); library(lme4)  # to run mixed models
library (lmerTest); library(DataCombine)
library(broman); library(broom.mixed)
library(stringr); library(ggplot2)
library(forcats); library(flextable); library(factoextra)
library(officer); library(purrr); library(sjPlot)
library(ftExtra); library(officedown); library(scales)
library(cowplot); library(performance)
library(Hmisc); library(pbkrtest); library(rstatix)
library(scales); library(rethinking)
library(brms); library(bestNormalize)
library(MASS); library(ggpubr); library(tidybayes)
library(rptR); library(ppcor); library(report)
## ---- end

# install.packages("remotes")
# remotes::install_github("stan-dev/cmdstanr")
# install.packages(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
# devtools::install_github("rmcelreath/rethinking")
# ---- load_data ----
rep1 <- readRDS(file = "data/processed/rep1.rds")
rep2 <- readRDS(file = "data/processed/rep2.rds")
rep3 <- readRDS(file = "data/processed/rep3.rds")
rep4 <- readRDS(file = "data/processed/rep4.rds")
total.data <- readRDS(file = "data/processed/total.data.rds")
earwig_egg_summary.2 <- readRDS(file = "data/processed/earwig_egg_summary.2.rds")
earwig_egg_summary.2 <- readRDS(file = "data/processed/earwig_egg_summary.2.rds")
earwig_egg_summary.3 <- readRDS(file = "data/processed/earwig_egg_summary.3.rds")
trade_off.brms_cov.mass <- readRDS(file = "data/processed/trade_off.brms_cov.mass.rds")
trade_off.brms_cov.pro <- readRDS(file = "data/processed/trade_off.brms_cov.pro.rds")
pheno.corr.r <- readRDS(file = "data/processed/pheno.corr.r.rds")
res_cor <- readRDS(file = "data/processed/res_cor.rds")
brood.va.delta <- readRDS(file = "data/processed/brood.va.delta.rds")
brood.vw.delta <- readRDS(file = "data/processed/brood.vw.delta.rds")
size.comparison <- readRDS(file = "data/processed/size.comparison.rds")
rep.num <- readRDS(file = "data/processed/rep.num.rds")
rep.size <- readRDS(file = "data/processed/rep.size.rds")

## ---- end

#### DATA READ and ANALYSIS ####
earwig_body<-read.csv("data/raw/earwig_body_measurements.csv")

# egg measurement file list
#earwig_egg_list<-read.csv("data/raw/female_pronotum_photos/pronotum_measurements.csv")
list_of_files <- list.files(path = "data/raw/egg_data_files",
                            recursive = TRUE,
                            pattern = "\\.csv$",
                            full.names = TRUE)

earwig_egg <- readr::read_csv(list_of_files)

#saveRDS(earwig_egg, file = "data/processed/earwig_egg.rds")
#saveRDS(earwig_body, file = "data/processed/earwig_body.rds")

earwig_body<-readRDS(file ="data/processed/earwig_body.rds")
earwig_egg<-readRDS(file ="data/processed/earwig_egg.rds")
#earwig_egg$brood <- as.numeric(earwig_egg$brood)
earwig_egg$num <- as.factor(earwig_egg$num)
earwig_egg$id <- as.factor(earwig_egg$id)

#### REPEATABILITY ####
# make wide body data long
data_long.body <- gather(earwig_body, trait, measurement, pronotum_length_1:body_mass_3, factor_key=TRUE)
data_long.body
data_long.body$brood <- as.numeric(data_long.body$brood)
data_long.body$num <- as.factor(data_long.body$num)
data_long.body$id <- as.factor(data_long.body$id)

# repeatability of pronotum length
data_long.body.pronotum <- data_long.body %>%
  filter(str_detect(trait, "pronotum_")) %>%
  separate(trait, into = c("trait", "part2"), sep = "_(?=[^_]*_)", extra = "merge")

rep1 <- rpt(measurement ~ 1 +  (1| id), grname = c("id"),
            data = data_long.body.pronotum, datatype = "Gaussian", nboot = 1000, npermut = 0)
saveRDS(rep1, file = "data/processed/rep1.rds")

# repeatability of body mass
data_long.body.mass <- data_long.body %>%
  filter(str_detect(trait, "body_")) %>%
  separate(trait, into = c("trait", "part2"), sep = "_(?=[^_]*_)", extra = "merge")

rep2 <- rpt(measurement ~ 1 +  (1| id), grname = c("id"),
            data = data_long.body.mass, datatype = "Gaussian", nboot = 1000, npermut = 0)
saveRDS(rep2, file = "data/processed/rep2.rds")

# calculate average pronotum length and body mass
earwig_body_2 <- earwig_body %>%
  rowwise() %>%
  mutate(mean.pronotum = mean(c_across(pronotum_length_1:pronotum_length_3)), mean.mass = mean(c_across(body_mass_1:body_mass_3))) %>%
  ungroup()

# make wide egg data long
data_long <- gather(earwig_egg, perimeter, measurement, perim:perim3, factor_key=TRUE)
data_long
data_long$brood <- as.numeric(data_long$brood)
data_long$num <- as.factor(data_long$num)
data_long$id <- as.factor(data_long$id)

# analyze Rpt for each brood separately
data_long_2 <- data_long %>%
  unite(peregg, c("num", "id"))

data_long_3 <- data_long_2 %>%
  filter(brood==1)

rep3 <- rpt(measurement ~ 1 +  (1| peregg), grname = c("peregg"),
            data = data_long_3, datatype = "Gaussian", nboot = 1000, npermut = 0)
saveRDS(rep3, file = "data/processed/rep3.rds")

data_long_4 <- data_long_2 %>%
  filter(brood==2)

rep4 <- rpt(measurement ~ 1 +  (1| peregg), grname = c("peregg"),
            data = data_long_4, datatype = "Gaussian", nboot = 1000, npermut = 0)
saveRDS(rep4, file = "data/processed/rep4.rds")

# calculate average egg perimeter
earwig_egg_2 <- earwig_egg %>%
  rowwise() %>%
  mutate(mean.perim = mean(c_across(perim:perim3))) %>%
  ungroup()
 
#### PRINCIPAL COMPONENTS ANALYSIS ####

pc <- prcomp(earwig_body_2[,c(-1, -2, -3,-4,-5, -6, -7,-8)],
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

saveRDS(earwig_body, file = "data/processed/earwig_body.rds")

###

egg.plot <- ggplot(brood_1_data, aes(x=num, y=measurement, colour=num)) +
  geom_point() +
  facet_wrap(~id, ncol=10)

plot(earwig_body_2$mean.pronotum, earwig_body_2$mean.mass)

# time between broods
earwig_egg.date <- earwig_egg %>%
  group_by(id, brood) %>%
  summarise(date=mean(date)) %>%
  mutate(brood = recode(brood, "1" = "one", "2" = "two")) %>%
  pivot_wider(names_from = brood, values_from = date, names_glue = "{brood}_{.value}")

data_with_diff <- earwig_egg.date %>%
  mutate(duration_days = as.numeric(difftime(two_date, one_date, units = "days"))) %>%
  ungroup() %>%
  summarise(mean = mean(duration_days), sd = sd(duration_days)) #77.8 d

# attach pc1 to egg data
total.data <- earwig_egg_2 %>%
  inner_join(earwig_body_2, by="id") %>%
  mutate(brood = recode(brood, "1" = "one", "2" = "two")) %>%
  ungroup() %>%
  dplyr::select(-c(num.x, num.y)) %>%
  group_by(id,brood) %>%
  mutate(num=n())

#total.data$brood.z <- (total.data$brood - 1)
total.data$mean.perim.z <- scale(total.data$mean.perim, center = TRUE, scale = TRUE)
total.data$num.z <- scale(total.data$num, center = TRUE, scale = TRUE)

saveRDS(total.data, file = "data/processed/total.data.rds")

earwig_egg_summary.2 <- total.data %>%
  mutate(brood = recode(brood, "1" = "one", "2" = "two")) %>%
  group_by(id, brood) %>%
  summarise(num=n(), egg.perim = mean(mean.perim),pc1=mean(pc1), mean.pro=mean(mean.pronotum), mean.massy=mean(mean.mass), sd.egg=sd(mean.perim), CV=sd.egg/egg.perim) 
saveRDS(earwig_egg_summary.2, file = "data/processed/earwig_egg_summary.2.rds")

earwig_egg_summary.3 <- total.data %>%
  dplyr::select(-perim,-perim2, -perim3) %>%
  mutate(brood = recode(brood, "1" = "one", "2" = "two")) %>%
  group_by(id, brood) %>%
  summarise(num=mean(num), mean.egg.perim = mean(mean.perim)) %>%
  group_by(brood) %>%
  dplyr::summarise(n=n(), mean.num=mean(num), sd.num=sd(num, na.rm=FALSE), min.num=min(num),
          max.num=max(num), egg.perim = mean(mean.egg.perim), sd.egg=sd(mean.egg.perim), min.size=min(mean.egg.perim),
          max.size=max(mean.egg.perim)) 
saveRDS(earwig_egg_summary.3, file = "data/processed/earwig_egg_summary.3.rds")

size.plot <- ggplot(earwig_egg_summary.2, aes(x=mean.massy, y=num, colour=brood)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x)

size.plot <- ggplot(earwig_egg_summary.2, aes(x=mean.pro, y=num, colour=brood)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x)

size.plot <- ggplot(earwig_egg_summary.2, aes(x=mean.pro, y=egg.perim, colour=brood)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x)
  
#### SCALED MASS ####

library(smatr)
size.data <- earwig_egg_summary.2 %>%
  filter(brood=="one")
fit<-sma(log(size.data$mean.massy)~log(size.data$mean.pro))
fit#no sex difference in slope p=0.821
summary(fit)
plot(fit)
plot(fit, which="qq")#assumptions are met

#reference population
L0<-mean(size.data$mean.pro)
L0
fit<-sma(log(size.data$mean.massy)~log(size.data$mean.pro))
fit
b.mass<-1.833134
Mi_hat<-(size.data$mean.massy)*(L0/size.data$mean.pro)^b.mass
size.data$Mi_hat <- Mi_hat
size.data.2 <- size.data %>%
  dplyr::select(id, Mi_hat)

total.data.1 <- total.data %>%
  left_join(size.data.2, by="id")





#### BRMS analysis ####

hist(earwig_egg_summary.2$num)
hist(total.data$mean.perim)

# PLOT OF RAW DATA

tradeoff.plot <- ggplot(earwig_egg_summary.2, aes(x=num, y=egg.perim, colour=brood, group=id)) +
  geom_line(colour="black") +
  geom_point(size=5) +
  ylab("Mean egg size (mm)") +
  xlab("Clutch size") +
  labs(color = "Clutch order") +
  scale_color_discrete(labels = c("First", "Second")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))
ggsave(tradeoff.plot, filename="figure_1.jpg", width=10.83, height=10.83, dpi=300,antialias="default")    

priors <- c(
  set_prior("normal(-1, 1)", class = "b", resp = c("scalemeanperim","scalenum")),
  set_prior("normal(-1, 1)", class = "Intercept", resp = c("scalemeanperim","scalenum")),
  set_prior("normal(0, 1)", class = "sd", resp = c("scalemeanperim","scalenum")),
  set_prior('lkj(2)', class = 'rescor'),
  set_prior("lkj(2)", class = "cor"))

# number_perim_cov.1 <- bf(mvbind(scale(mean.perim),scale(num)) ~ 1 + brood + (1|p|id)) + set_rescor(TRUE) #
# trade_off.brms_cov.1 <- brm(number_perim_cov.1, data = total.data,
#                             family = gaussian(),
#                               cores = 6,
#                               chains = 6,
#                               warmup = 3000,
#                               iter = 10000,
#                               thin = 2,
#                               prior=priors,
#                               file = "data/processed/trade_off.brms_cov.1",
#                               control = list(adapt_delta = 0.999),
#                               seed = 12345)
# summary(trade_off.brms_cov.1)
# get_variables(trade_off.brms_cov.1)
# trade_off.brms_cov.1 <- add_criterion(trade_off.brms_cov.1, "loo")

# need to see if model controlling for body size is better than null model and determine which body size proxy is best
# number_perim_cov.2.mass <- bf(mvbind(scale(mean.perim), scale(num)) ~ 1 + brood + scale(mean.mass) + (1|p|id)) + set_rescor(TRUE) # this model is the best
# trade_off.brms_cov.mass <- brm(number_perim_cov.2.mass, data = total.data,
#                             family = gaussian(),
#                             cores = 6,
#                             chains = 6,
#                             thin=2,
#                             warmup = 3000,
#                             iter = 10000,
#                             prior=priors,
#                             control = list(adapt_delta = 0.999),
#                             file = "data/processed/trade_off.brms_cov.mass",
#                             seed = 12345)
# summary(trade_off.brms_cov.mass)
# trade_off.brms_cov.mass <- add_criterion(trade_off.brms_cov.mass, "loo")
# get_variables(trade_off.brms_cov.mass)
# tidy(trade_off.brms_cov.mass)

# random regression (brood|id) says that brood slope varies for each level of id
number_perim_cov.2.pro <- bf(mvbind(scale(mean.perim), scale(num)) ~ 1 + brood + scale(mean.pronotum) + (1|a|id)) + set_rescor(TRUE)
trade_off.brms_cov.pro <- brm(number_perim_cov.2.pro, data = total.data,
                               family = gaussian(),
                               cores = 6,
                               chains = 6,
                               thin=2,
                               warmup = 3000,
                               iter = 10000,
                               prior=priors,
                               sample_prior = TRUE,
                               control = list(adapt_delta = 0.999),
                               file = "data/processed/trade_off.brms_cov.pro",
                               seed = 12345)
summary(trade_off.brms_cov.pro)
trade_off.brms_cov.pro <- add_criterion(trade_off.brms_cov.pro, "loo")

pp_check(trade_off.brms_cov.pro, resp = "scalenum", ndraws = 1000)
pp_check(trade_off.brms_cov.pro, resp = "scalemeanperim", ndraws = 1000)


# number_perim_cov.2.pc1 <- bf(mvbind(scale(mean.perim), scale(num)) ~ 1 + brood + pc1 + (1|p|id)) + set_rescor(TRUE) # mean.mass is best
# trade_off.brms_cov.pc1 <- brm(number_perim_cov.2.pc1, data = total.data,
#                               family = gaussian(),
#                               cores = 6,
#                               chains = 6,
#                               thin=2,
#                               warmup = 3000,
#                               iter = 10000,
#                               prior=priors,
#                               sample_prior = TRUE,
#                               control = list(adapt_delta = 0.999),
#                               file = "data/processed/trade_off.brms_cov.pc1",
#                               seed = 12345)
# summary(trade_off.brms_cov.pc1)
# trade_off.brms_cov.pc1 <- add_criterion(trade_off.brms_cov.pc1, "loo")

# number_perim_cov.2.smi <- bf(mvbind(scale(mean.perim), scale(num)) ~ 1 + brood + scale(Mi_hat) + (1|p|id)) + set_rescor(TRUE) # mean.mass is best
# trade_off.brms_cov.smi <- brm(number_perim_cov.2.smi, data = total.data.1,
#                               family = gaussian(),
#                               cores = 6,
#                               chains = 6,
#                               thin=2,
#                               warmup = 3000,
#                               iter = 10000,
#                               prior=priors,
#                               sample_prior = TRUE,
#                               control = list(adapt_delta = 0.999),
#                               file = "data/processed/trade_off.brms_cov.smi",
#                               seed = 12345)
# summary(trade_off.brms_cov.smi)
# trade_off.brms_cov.smi <- add_criterion(trade_off.brms_cov.smi, "loo")

size.comparison <- loo_compare(trade_off.brms_cov.mass, trade_off.brms_cov.pro, trade_off.brms_cov.pc1, trade_off.brms_cov.smi,trade_off.brms_cov.1, criterion = "loo")
size.comparison.report <-report(size.comparison)

saveRDS(size.comparison, file = "data/processed/size.comparison.rds")

pp_check(trade_off.brms_cov.pro, resp="scalenum",ndraws = 100)
pp_check(trade_off.brms_cov.pro, resp="scalemeanperim",ndraws = 100)

VarCorr(trade_off.brms_cov.smi)
# r <- -0.1298779 + 0.05960063 # = -0.07
# cor.test(~num + egg.perim, method = "pearson", data=earwig_egg_summary.2) # no trade-off = -0.04
# these two are very similar

pheno.corr.r <- trade_off.brms_cov.smi %>%
  spread_draws(rescor__scalemeanperim__scalenum, cor_id__scalemeanperim_Intercept__scalenum_Intercept) %>%
  mean_qi(r=rescor__scalemeanperim__scalenum+cor_id__scalemeanperim_Intercept__scalenum_Intercept)

saveRDS(pheno.corr.r, file = "data/processed/pheno.corr.r.rds")

res_cor <- trade_off.brms_cov.smi %>%
  spread_draws(rescor__scalemeanperim__scalenum) %>%
  mean_qi(r=rescor__scalemeanperim__scalenum)

saveRDS(res_cor, file = "data/processed/res_cor.rds")

# correlation between clutch and egg size per brood
brood.one.data <- total.data %>%
  filter(brood=="one") %>%
  ungroup()

number_perim_cov.2.pro.1 <- bf(mvbind(scale(mean.perim)) ~ 1 + num + scale(mean.pronotum), family=gaussian) # 
trade_off.brms_cov.pro.1 <- brm(number_perim_cov.2.pro.1, data = brood.one.data,
                              cores = 6,
                              chains = 6,
                              thin=2,
                              warmup = 300,
                              iter = 1000,
                              #prior=priors,
                              sample_prior = TRUE,
                              control = list(adapt_delta = 0.999),
                              #file = "data/processed/trade_off.brms_cov.pro",
                              seed = 12345)
summary(trade_off.brms_cov.pro.1)

brood.two.data <- total.data %>%
  filter(brood=="two") %>%
  ungroup()

number_perim_cov.2.pro.2 <- bf(mvbind(scale(mean.perim)) ~ 1 + num + scale(mean.pronotum), family=gaussian) # 
trade_off.brms_cov.pro.2 <- brm(number_perim_cov.2.pro.2, data = brood.two.data,
                                cores = 6,
                                chains = 6,
                                thin=2,
                                warmup = 300,
                                iter = 1000,
                                #prior=priors,
                                sample_prior = TRUE,
                                control = list(adapt_delta = 0.999),
                                #file = "data/processed/trade_off.brms_cov.pro",
                                seed = 12345)
summary(trade_off.brms_cov.pro.2)

# intercept.plot <- trade_off.brms_cov.mass %>%
#       spread_draws(r_id__scalemeanperim[id,term], r_id__scalenum[id, term]) %>%
#       group_by(id, term) %>%   # this line not necessary (done by spread_draws)
#       mean_qi(r_id__scalenum,r_id__scalemeanperim) %>%
#       ggplot(aes(x=r_id__scalenum, y=r_id__scalemeanperim)) +
#       geom_point(size=4) +
#   ylab("Egg size (Intercept)") +
#   xlab("Clutch size (Intercept)") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill="white"),
#         axis.title.x = element_text(size=14,face="bold"),
#         axis.title.y = element_text(size=14,face="bold"),
#         axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12))
# ggsave(intercept.plot, filename="figure_2.jpg", width=10.83, height=10.83, dpi=300,antialias="default")    

#### AMONG AND WITHIN EGG SIZE VARIATION ####
egg.size.bf = bf(scale(mean.perim) ~ mean.pronotum + brood + (1|a|id) + (0+brood||id), sigma ~ brood + (1|a|id)) # need to have (0+brood||id) in order to get sd_id_broodone and 
# sd_id_broodtwo

#run the model
egg.size.model <- brm(egg.size.bf,
                                data = total.data.1,
                                family = gaussian,
                                cores = 6, 
                                chains = 6, 
                                warmup = 1000,
                                iter = 4000,
                                seed = 34, #make sure to set the seed to make results reproducible
                                #file = "data/processed/egg.size.model.rds",
                                control = list(adapt_delta = 0.999))
summary(egg.size.model)
get_variables(egg.size.model)
pp_check(egg.size.model,ndraws = 100)

#among
Va.brood.one <- as_draws_df(egg.size.model)$"sd_id__broodone"^2
Va.brood.two <- as_draws_df(egg.size.model)$"sd_id__broodtwo"^2

#within
Vw.brood.one <- exp(as_draws_df(egg.size.model)$"b_sigma_broodone")^2
Vw.brood.two <- exp(as_draws_df(egg.size.model)$"b_sigma_broodtwo")^2

post.data.brood = data.frame(Va.brood.one, Va.brood.two, Vw.brood.one, Vw.brood.two) 

# among-female egg size variation
brood.Va = post.data.brood %>%
  dplyr::select(starts_with("Va.brood.")) %>%
  pivot_longer(cols = starts_with("Va.brood."),
               names_to = 'brood',
               names_prefix = "Va.brrod.",
               values_to ="Estimate")

Va_brood = brood.Va %>%
  dplyr::group_by(brood) %>%
  dplyr::summarise(Va_brood = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3))%>%
  as.data.frame()

post.data.brood$delta.va.brood=
  with(post.data.brood, Va.brood.one-Va.brood.two)

va.delta.brood = post.data.brood %>%
  dplyr::select(starts_with("delta.va.brood")) %>%
  pivot_longer(cols = starts_with("delta.va.brood"),
               names_to = 'Contrast',
               names_prefix = "delta.va.brood",
               values_to ="Estimate")

brood.va.delta = va.delta.brood %>%
  dplyr::group_by(Contrast)%>%
  dplyr::summarise(va.delta.brood = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3)) %>%
  as.data.frame()

saveRDS(brood.va.delta, file = "data/processed/brood.va.delta.rds")

# within-female egg size variation
brood.Vw = post.data.brood %>%
  dplyr::select(starts_with("Vw.brood.")) %>%
  pivot_longer(cols = starts_with("Vw.brood."),
               names_to = 'brood',
               names_prefix = "Vw.brrod.",
               values_to ="Estimate")

Vw_brood = brood.Vw %>%
  dplyr::group_by(brood) %>%
  dplyr::summarise(Vw_brood = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3))%>%
  as.data.frame()

post.data.brood$delta.vw.brood=
  with(post.data.brood, Vw.brood.one-Vw.brood.two)

vw.delta.brood = post.data.brood %>%
  dplyr::select(starts_with("delta.vw.brood")) %>%
  pivot_longer(cols = starts_with("delta.vw.brood"),
               names_to = 'Contrast',
               names_prefix = "delta.vw.brood",
               values_to ="Estimate")

brood.vw.delta = vw.delta.brood %>%
  dplyr::group_by(Contrast)%>%
  dplyr::summarise(vw.delta.brood = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3)) %>%
  as.data.frame()

saveRDS(brood.vw.delta, file = "data/processed/brood.vw.delta.rds")

# repeatability of egg size
rep.size <- rpt(egg.perim ~ 1 + + (1|id), grname = c("id"),
            data = earwig_egg_summary.2, datatype = "Gaussian", nboot = 1000, npermut = 0)
saveRDS(rep.size, file = "data/processed/rep.size.rds")

total.data.number <- total.data %>%
  group_by(id, brood) %>%
  filter(row_number()==1)

rep.num <- rptPoisson(num ~ 1 + (1|id), grname = c("id"),
                data = earwig_egg_summary.2, link="log", nboot = 1000, npermut = 0)
saveRDS(rep.num, file = "data/processed/rep.num.rds")

# egg size variation (CV) vs clutch size
brood.one <- earwig_egg_summary.2 %>%
  filter(brood=="one")
cor.test(brood.one$CV, brood.one$num, method="spearman")

brood.two <- earwig_egg_summary.2 %>%
  filter(brood=="two")
cor.test(brood.two$sd.egg, brood.two$num, method="spearman")

##PLOTS

var.among.brood <- post.data.brood %>%
  dplyr::select(Va.brood.one, Va.brood.two) %>%
  gather(brood, value,  factor_key=TRUE) %>%
  mutate(order=if_else(str_detect(brood, '.one'), 'first', 'second')) %>%
  mutate(brood = str_remove(brood, "Va.brood."))

plot.Va.brood <- var.among.brood %>%
  mutate(clutch = factor(order, levels=c("first","second"))) %>%
  ggplot(aes(x = order, y = value)) +
  stat_dotsinterval(slab_fill="darkgrey", slab_color="darkgrey", point_interval = median_qi,
                    .width = 0.95, quantiles=100) +
  ylab("Among-female egg size variation") +
  xlab("Clutch laying order") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))
#ggsave(filename="figure_1.tiff", width=6.83, height=6.83, dpi=300,antialias="default")

var.within.brood <- post.data.brood %>%
  dplyr::select(Vw.brood.one, Vw.brood.two) %>%
  gather(brood, value,  factor_key=TRUE) %>%
  mutate(order=if_else(str_detect(brood, '.one'), 'first', 'second')) %>%
  mutate(brood = str_remove(brood, "Vw.brood."))

plot.Vw.brood <- var.within.brood %>%
  mutate(clutch = factor(order, levels=c("first","second"))) %>%
  ggplot(aes(x = order, y = value)) +
  stat_dotsinterval(slab_fill="darkgrey", slab_color="darkgrey", point_interval = median_qi,
                    .width = 0.95, quantiles=100) +
  ylab("Within-clutch egg size variation") +
  xlab("Clutch laying order") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))
#ggsave(filename="figure_1.tiff", width=6.83, height=6.83, dpi=300,antialias="default")

figure_2<-plot_grid(plot.Va.brood, plot.Vw.brood, ncol=2, nrow=1,
                    labels = c('(a)', '(b)'))
ggsave(filename="figure_2.jpg", width=10.83, height=10.83, dpi=300,antialias="default")    

### KOCH
koch<-read.csv("data/raw/koch_data.csv")

koch <- koch %>%
  filter(!is.na(size))

koch.one <- koch %>%
  filter(clutch==1)
koch.two <- koch %>%
  filter(clutch==2)
cor.test(koch$number, koch$size, method="spearman")
cor.test(koch.one$number, koch.one$size, method="pearson")
cor.test(koch.two$number, koch.two$size, method="spearman")

#### junk ####
# tests


egg.size.bf.1 = bf(scale(mean.perim) ~ mean.mass + brood + (1|id) + (0+brood||id), sigma ~ 0 + brood + (1|id))
#run the model
egg.size.model.1 <- brm(egg.size.bf.1,
                      data = total.data,
                      family = gaussian,
                      cores = 6, 
                      chains = 6, 
                      warmup = 1000,
                      iter = 4000,
                      seed = 34, #make sure to set the seed to make results reproducible
                      #file = "data/processed/egg.size.model.1.rds",
                      control = list(adapt_delta = 0.95))
summary(egg.size.model.1)
get_variables(egg.size.model.1)

egg.size.bf.2 = bf(scale(mean.perim) ~ mean.mass + brood + (1|id) + (0+brood|id), sigma ~ 0 + brood + (1|id)) ## very interested in this one; random intercept for id in both part sof model and egg size in clutch one is correlated with clutch two
#run the model
egg.size.model.2 <- brm(egg.size.bf.2,
                        data = total.data,
                        family = gaussian,
                        cores = 6, 
                        chains = 6, 
                        warmup = 1000,
                        iter = 4000,
                        seed = 34, #make sure to set the seed to make results reproducible
                        #file = "data/processed/egg.size.model.2.rds",
                        control = list(adapt_delta = 0.95))
summary(egg.size.model.2)
get_variables(egg.size.model.2)

ranef(egg.size.model.2)

egg.size.bf.3 = bf(scale(mean.perim) ~ mean.mass + brood + (0+brood|id), sigma ~ 0 + brood + (1|id))
#run the model
egg.size.model.3 <- brm(egg.size.bf.3,
                        data = total.data,
                        family = gaussian,
                        cores = 6, 
                        chains = 6, 
                        warmup = 1000,
                        iter = 4000,
                        seed = 34, #make sure to set the seed to make results reproducible
                        #file = "data/processed/egg.size.model.3.rds",
                        control = list(adapt_delta = 0.95))
summary(egg.size.model.3)
get_variables(egg.size.model.3)


egg.size.bf.4 = bf(scale(mean.perim) ~ mean.mass + brood + (0+brood|id), sigma ~ 0 + brood)
#run the model
egg.size.model.4 <- brm(egg.size.bf.4,
                        data = total.data,
                        family = gaussian,
                        cores = 6, 
                        chains = 6, 
                        warmup = 1000,
                        iter = 4000,
                        seed = 34, #make sure to set the seed to make results reproducible
                        #file = "data/processed/egg.size.model.4.rds",
                        control = list(adapt_delta = 0.95))
summary(egg.size.model.4)
get_variables(egg.size.model.4)

egg.size.bf.6 = bf(scale(mean.perim) ~ mean.mass + brood + (1+brood|id), sigma ~ 0 + brood)
#run the model
egg.size.model.6 <- brm(egg.size.bf.6,
                        data = total.data,
                        family = gaussian,
                        cores = 6, 
                        chains = 6, 
                        warmup = 1000,
                        iter = 4000,
                        seed = 34, #make sure to set the seed to make results reproducible
                        #file = "data/processed/egg.size.model.4.rds",
                        control = list(adapt_delta = 0.95))
summary(egg.size.model.6)
get_variables(egg.size.model.6)


egg.size.bf.5 = bf(scale(mean.perim) ~ mean.mass + brood + (1|a|id) + (0+brood|id), sigma ~ 0 + brood + (1|a|id)) ## very interested in this one; random intercept for id in both part sof model and egg size in clutch one is correlated with clutch two
#run the model
egg.size.model.5 <- brm(egg.size.bf.5,
                        data = total.data,
                        family = gaussian,
                        cores = 6, 
                        chains = 6, 
                        warmup = 1000,
                        iter = 4000,
                        seed = 34, #make sure to set the seed to make results reproducible
                        #file = "data/processed/egg.size.model.2.rds",
                        control = list(adapt_delta = 0.95))
summary(egg.size.model.5)
get_variables(egg.size.model.5)

egg.size.bf.7 = bf(scale(mean.perim) ~ mean.mass + brood + (1+brood||id), sigma ~ 0 + brood) ## very interested in this one; random intercept for id in both part sof model and egg size in clutch one is correlated with clutch two
#run the model
egg.size.model.7 <- brm(egg.size.bf.7,
                        data = total.data,
                        family = gaussian,
                        cores = 6, 
                        chains = 6, 
                        warmup = 1000,
                        iter = 4000,
                        seed = 34, #make sure to set the seed to make results reproducible
                        #file = "data/processed/egg.size.model.2.rds",
                        control = list(adapt_delta = 0.95))
summary(egg.size.model.7)
get_variables(egg.size.model.7)








nd <- tibble(x_s = seq(from = -3, to = 3, length.out = d %>% nrow()))

fitted(trade_off.brms_cov.2a) %>% 
  data.frame() %>% 
  bind_cols(nd) %>% 

number_perim_2 <- bf(mvbind(scale(mean.perim),scale(num)) ~ 1 + brood + (1 + brood|p|id)) + set_rescor(TRUE) # gives two values for each id
trade_off.brms_2 <- brm(number_perim_cov.2, data = total.data.2,
                            cores = 4,
                            chains = 4,
                            warmup = 300,
                            iter = 600,
                            thin = 2,
                            seed = 12345)
summary(trade_off.brms_2)
get_variables(trade_off.brms_2, summary=TRUE)


coef(trade_off.brms_cov.2a) 

corr_draws <- trade_off.brms_cov.2 %>%
  spread_draws(Intercept_num) 

posterior_draws <- as_draws_df(trade_off.brms_cov.2)

VarCorr(trade_off.brms_cov.2)

trade_off.brms_cov.2a %>%
  spread_draws(r_id__scalemeanperim[id,term], r_id__scalenum[id, term],Intercept_scalenum, Intercept_scalemeanperim, sigma_scalenum, sigma_scalemeanperim) %>%
  group_by(id, term) %>%   # this line not necessary (done by spread_draws)
  mean_qi(Intercept_scalenum,sigma_scalemeanperim,sigma_scalenum) %>%
  #filter(id=="8WOPO") %>%
  ggplot(aes(x=mean.num, y=mean.perim, group=id)) +
  geom_point() +
  geom_line()

# egg size
df_1 <-
  # with this line we select each of the 20 cafe's posterior mean (i.e., Estimate)
  # for both `Intercept` and `afternoon`
  coef(trade_off.brms_cov.2a)$id[ , 1, 1:2] %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  rename(one = scalemeanperim_Intercept) %>%
  mutate(two = one + scalemeanperim_brood) %>%
  dplyr::select(c(-scalemeanperim_brood)) %>%
  pivot_longer(cols = one:two, 
               names_to = "brood",
               values_to = "size")
# brood size 
df_2 <-
  # with this line we select each of the 20 cafe's posterior mean (i.e., Estimate)
  # for both `Intercept` and `afternoon`
  coef(trade_off.brms_cov.2a)$id[ , 1, 3:4] %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  rename(one = scalenum_Intercept) %>%
  mutate(two = one + scalenum_brood) %>%
  dplyr::select(c(-scalenum_brood)) %>%
  pivot_longer(cols = one:two, 
               names_to = "brood",
               values_to = "number")

egg.corr <- df_1 %>%
  left_join(df_2, by= c("id", "brood"))

ggplot(egg.corr, aes(x=number, y=size, group=id)) +
  geom_point() +
  geom_line()

ggplot(earwig_egg_summary.2, aes(x=num, y=egg.perim, colour=brood, group=id)) +
  geom_point() +
  geom_line()

# brood one
cor.test(~num + egg.perim, method = "spearman", data=earwig_egg_summary.2,  subset = (brood == "one")) # no trade-off
cor.test(~number + size, method = "spearman", data=egg.corr,  subset = (brood == "one")) # positive
corr.data1 <- earwig_egg_summary.2 %>%
  filter(brood=="one")
  pcor.test(corr.data1$num, corr.data1$egg.perim, corr.data1$pc1, method = "spearman") # no trade-off

# brood two
cor.test(~num + egg.perim, method = "spearman", data=earwig_egg_summary.2,  subset = (brood == "two")) # trade-off
cor.test(~number + size, method = "spearman", data=egg.corr,  subset = (brood == "two")) # no trade-off
corr.data2 <- earwig_egg_summary.2 %>%
  filter(brood=="two")
pcor.test(corr.data2$num, corr.data2$egg.perim, corr.data2$pc1, method = "spearman") # trade-off

# pooled
cor.test(~num + egg.perim, method = "pearson", data=earwig_egg_summary.2) # no trade-off
cor.test(~number + size, method = "spearman", data=egg.corr) # positve
pcor.test(earwig_egg_summary.2$num, earwig_egg_summary.2$egg.perim, earwig_egg_summary.2$pc1, method = "spearman") # no trade-off

# among-individual correlation
draws.among.corr.brms.cov <- trade_off.brms_cov.2 %>%
  gather_draws(
    cor_id__meanperim_Intercept__num_Intercept
  ) %>%
  mean_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__meanperim_Intercept__num_Intercept"  ~ "size—number"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval) %>%
  mutate(category = rep(c('among'), 15, length.out = n()))

# within-individual correlation
draws.within.corr.brms.cov <- trade_off.brms_cov.2 %>%
  gather_draws(
    rescor__meanperim__num
  ) %>%
  mean_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "rescor__meanperim__num"  ~ "size—number"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval) %>%
  mutate(category = rep(c('within'), 1, length.out = n())) 

number_perim_cov.4 <- bf(mvbind(mean.perim,num) ~ 1 + pc1 + brood + (1 + brood|p|id)) + set_rescor(FALSE) # gives two values for each id
trade_off.brms_cov.4 <- brm(number_perim_cov.4, data = total.data.2,
                            family = c(gaussian, gaussian),
                            cores = 4,
                            chains = 4,
                            warmup = 300,
                            iter = 600,
                            thin = 2,
                            seed = 12345)
summary(trade_off.brms_cov.4)
get_variables(trade_off.brms_cov.4)


VarCorr(trade_off.brms_cov.4)

cov.num.perim <- -0.258466479
#residual + individual
num.var <- 5.305397e-05 + 14.23216253 + 14.09075647 
perim.var <- 2.695111e-01 + 0.05795814 + 0.06313121

cov.num.perim/sqrt(num.var*perim.var) # phenotypic corr

# raw correlation
number_perim.cov.raw <- bf(mvbind(mean.perim,num) ~ 1 + pc1 + brood + (1 + brood|p|id)) + set_rescor(FALSE)
trade_off.brms.cov.raw <- brm(number_perim.cov.raw, data = total.data.2,
                              family=c(gaussian,poisson),
                      cores = 4,
                      chains = 4,
                      warmup = 300,
                      iter = 600,
                      thin = 2,
                      seed = 12345)
summary(trade_off.brms.raw)
get_variables(trade_off.brms.raw)

draws.raw.corr.brms.cov <- trade_off.brms.raw %>%
  gather_draws(
    cor_id__meanperim_Intercept__num_Intercept
  ) %>%
  mean_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__meanperim_Intercept__num_Intercept"  ~ "size—number"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval) %>%
  mutate(category = rep(c('raw'), 15, length.out = n()))
