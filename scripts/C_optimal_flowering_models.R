#####################################################
##                                                 ##
##             Optimal Flowering Time              ##
##                                                 ##
##       Linear and Inverse model comparison       ##
##                                                 ##
##            JJ & JP - 4th Nov 2024               ##
##                                                 ##
#####################################################

# Solidago and Lythrum
# Using actual flowering date rather than relative day
# Capped growing season length (5 consecutive days > 5C as start of growing season - converse for end)
# brms model comparing constant f to Iwasa-Cohen model of varying f

rm(list = ls())
options(width = 100)

library(tidyverse)
library(brms)
library(latex2exp) 
library(patchwork)
library(flextable)
library(cowplot)
library(gridExtra)

#_______________________________________________________________________________
#### 1. load data and wrangle ####

## Lythrum
load("data/lythrum_99.RData")

# first flowering date to adjust
lythrum_first_flower <- 180 - 137 # first flower = 29 June (day 180); umea sos in 1999 = 137

# data cleaning
lythrum_mdat <- filter(lythrum99, capped_season_length > 0) %>% 
  mutate(GS_f = paste0(round(capped_season_length, 1), "_days"),
         population = paste0("pop", pop)) %>% 
  mutate(population = factor(population, levels = paste0("pop", unique(lythrum99$pop)[order(unique(lythrum99$pop))])),
         flowering_start = flowering_start + (lythrum_first_flower - 1))

## Solidago
load("data/solidago_05.RData")

# first flowering date to adjust data
solidago_first_flower <- 160 - 91 # first flower = 9 June (day 160); uppsala sos in 2005 = 91

# data cleaning
solidago_mdat <- solidago05 %>% 
  mutate(GS_f = paste0(round(capped_season_length, 1), "_days"),
         population = paste0("pop", pop)) %>% 
  mutate(population = factor(population, levels = paste0("pop", unique(solidago05$pop)[order(unique(solidago05$pop))])),
         flowering_start = flowering_start + (solidago_first_flower - 1))

## Colours
lythrum_palette <- c("#8C4681", "#A6449F", "#172601", "#3E5902", "#567639")
lythrum_colour <- lythrum_palette[2]

solidago_palette <- c("#0e2d10", "#214d08", "#90980c", "#d2c200", "#fce708")
solidago_colour <- solidago_palette[4]

#_______________________________________________________________________________
#### 2. Data checks ####

## Distribution checks
ggplot(lythrum_mdat, aes(x = capped_season_length)) +
  geom_histogram(bins = 10)

ggplot(lythrum_mdat, aes(x = flowering_start)) +
  geom_histogram(bins = 35) +
  geom_vline(aes(xintercept = mean(flowering_start))) +
  facet_wrap(~ GS_f, nrow = 7) +
  theme_bw() 

ggplot(solidago_mdat, aes(x = capped_season_length)) +
  geom_histogram(bins = 10)

ggplot(solidago_mdat, aes(x = flowering_start)) +
  geom_histogram(bins = 35) +
  geom_vline(aes(xintercept = mean(flowering_start))) +
  facet_wrap(~ GS_f, nrow = 7) +
  theme_bw() 

## Raw data plots
rawplot_lythrum <- ggplot(lythrum_mdat, 
       aes(x = capped_season_length, 
           y = flowering_start)) +
  geom_jitter(size = 3, alpha = 0.7, colour = lythrum_colour) +
  stat_summary(geom = "errorbar", 
               aes(group = population, colour = NULL), 
               size = 1, 
               width = 0.4) +  
  stat_summary(geom = "point", 
               aes(group = population, colour = NULL), 
               size = 4) +
  labs(x = expression(paste("Growing Season Length (", italic("T"), ") [days]")), 
       y = expression(paste("Flowering time from start of season (", italic(t[italic("2")]), ")")),
       title = "Lythrum") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "italic"))

rawplot_solidago <- ggplot(solidago_mdat, 
                          aes(x = capped_season_length, 
                              y = flowering_start)) +
  geom_jitter(size = 3, alpha = 0.7, colour = solidago_colour) +
  stat_summary(geom = "errorbar", 
               aes(group = population, colour = NULL), 
               size = 1, 
               width = 0.4) +  
  stat_summary(geom = "point", 
               aes(group = population, colour = NULL), 
               size = 4) +
  labs(x = expression(paste("Growing Season Length (", italic("T"), ") [days]")),
       y = expression(paste("Flowering time from start of season (", italic(t[italic("2")]), ")")),
       title = "Solidago") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "italic"))

rawplot_combined <- grid.arrange(rawplot_lythrum, rawplot_solidago, ncol = 1)
ggsave(rawplot_combined, 
       filename = "output/rawplot_optimal_flowering20241010.jpeg",
       width = 12, height = 20, units = "cm", dpi = 600)

#_______________________________________________________________________________
#### 3. Constant f model ####

set.seed(999)
lythrum_constant_f <- brm(bf(flowering_start ~ capped_season_length - (1/f),
                           f ~ 1, nl = TRUE),
                        data = lythrum_mdat, 
                        family = gaussian,
                        prior = c(prior(normal(0, 0.1), nlpar =  f, lb = 0),
                                  prior(normal(0, 5), class = sigma, lb = 0)), 
                        control = list(adapt_delta = 0.97, max_treedepth = 12),
                        chains = 4, cores = 4, 
                        iter = 4000, warmup = 2000)

set.seed(999)
solidago_constant_f <- brm(bf(flowering_start ~ capped_season_length - (1/f),
                             f ~ 1, nl = TRUE),
                          data = solidago_mdat, 
                          family = gaussian,
                          prior = c(prior(normal(0, 0.1), nlpar =  f, lb = 0),
                                    prior(normal(0, 5), class = sigma, lb = 0)), 
                          control = list(adapt_delta = 0.97, max_treedepth = 12),
                          chains = 4, cores = 4, 
                          iter = 4000, warmup = 2000)

#_______________________________________________________________________________
#### 4. Iwasa-Cohen model ####

# example priors (same for both species)
get_prior(bf(flowering_start ~ capped_season_length - (1/(a/(1 + exp(-b*(capped_season_length + c))))),
             a + b + c ~ 1, nl = TRUE), data = lythrum_mdat, family = gaussian)

## Lythrum
set.seed(999)
lythrum_iwasa_cohen <- 
  brm(bf(flowering_start ~ capped_season_length - (1/(a/(1 + exp(-b*(capped_season_length + c))))),
         a + b + c ~ 1, nl = TRUE),
      data = lythrum_mdat, family = gaussian,
      prior = c(prior(normal(0.05, 1), nlpar =  a, lb = 0.0001, ub = 1),
                prior(normal(-0.05, 1), nlpar =  b, lb = -1, ub = -0.0001),
                prior(normal(-150, 10), nlpar = c, lb = -200, ub = -100),
                prior(normal(0, 5), class = sigma, lb = 0)), 
      control = list(adapt_delta = 0.97, max_treedepth = 12),
      chains = 4, cores = 4, 
      iter = 4000, warmup = 2000)

## Solidago
set.seed(999)
solidago_iwasa_cohen <- 
  brm(bf(flowering_start ~ capped_season_length - (1/(a/(1 + exp(-b*(capped_season_length + c))))),
         a + b + c ~ 1, nl = TRUE),
      data = solidago_mdat, 
      family = gaussian,
      prior = c(prior(normal(0.5, 0.75), nlpar =  a, lb = 0.0001, ub = 1),
                prior(normal(-0.5, 0.75), nlpar =  b, lb = -1, ub = -0.0001),
                prior(normal(-150, 10), nlpar = c, lb = -200, ub = -100),
                prior(normal(0, 5), class = sigma, lb = 0)), 
      control = list(adapt_delta = 0.97, max_treedepth = 12),
      chains = 4, cores = 4, 
      iter = 4000, warmup = 2000)

#_______________________________________________________________________________
#### 5. Model comparisons ####

## Lythrum
lythrum_constant_f <- add_criterion(lythrum_constant_f, criterion = c("loo","waic"))
lythrum_iwasa_cohen <- add_criterion(lythrum_iwasa_cohen, criterion = c("loo","waic"))

as.data.frame(loo_compare(lythrum_constant_f, lythrum_iwasa_cohen, criterion = "loo")) 

as.data.frame(loo_compare(lythrum_constant_f, lythrum_iwasa_cohen, criterion = "loo")) %>% 
  mutate(Model = c("Iwasa-Cohen model [f(T)]", "Constant f")) %>%
  dplyr::select(Model, elpd = elpd_loo, elpd_diff,
                `elpd se` = se_elpd_loo, `LOOIC` = looic) %>%
  mutate(across(elpd:LOOIC, \(x) round(x, 1))) %>%
  flextable(cwidth = 1) %>%
  width(j = 1, width = 2.5) %>%
  compose(part = "header", j = "Model", value = as_paragraph("Lythrum model")) %>%
  compose(part = "header", j = "elpd_diff", value = as_paragraph("\U0394", "elpd")) %>%
  theme_zebra(odd_header = alpha(lythrum_colour, alpha = 0.3)) %>%
  save_as_image("output/model_comparison_lythrum_20241104.png")

## Solidago
solidago_constant_f <- add_criterion(solidago_constant_f, criterion = c("loo","waic"))
solidago_iwasa_cohen <- add_criterion(solidago_iwasa_cohen, criterion = c("loo","waic"))

as.data.frame(loo_compare(solidago_constant_f, solidago_iwasa_cohen, criterion = "loo")) 

as.data.frame(loo_compare(solidago_constant_f, solidago_iwasa_cohen, criterion = "loo")) %>% 
  mutate(Model = c("Iwasa-Cohen model [f(T)]", "Constant f")) %>%
  dplyr::select(Model, elpd = elpd_loo, elpd_diff,
                `elpd se` = se_elpd_loo, `LOOIC` = looic) %>%
  mutate(across(elpd:LOOIC, \(x) round(x, 1))) %>%
  flextable(cwidth = 1) %>%
  width(j = 1, width = 2.5) %>%
  compose(part = "header", j = "Model", value = as_paragraph("Solidago model")) %>%
  compose(part = "header", j = "elpd_diff", value = as_paragraph("\U0394", "elpd")) %>%
  theme_zebra(odd_header = alpha(solidago_colour, alpha = 0.3)) %>%
  save_as_image("output/model_comparison_solidago_20241104.png")

#_______________________________________________________________________________
#### 6. Posterior draws for f ####

#### 6a. Lythrum
post_draws_f_lythrum <- 
  as_draws_df(lythrum_iwasa_cohen) %>% 
  as_tibble() %>% 
  dplyr::select(1:3) %>% 
  slice(rep(1:n(), each = 100)) %>% 
  mutate(capped_season_length = rep(seq(134,206, length.out = 100), 8000),
         f = b_a_Intercept/(1 + exp(-b_b_Intercept*(capped_season_length + b_c_Intercept))))

f_lythrum <- post_draws_f_lythrum %>% 
  group_by(capped_season_length) %>% 
  summarise(f_mn = mean(f),
            f_50l = quantile(f, prob = 0.25),
            f_50u = quantile(f, prob = 0.75),
            f_min = quantile(f, prob = 0.025),
            f_max = quantile(f, prob = 0.975))

f_T_lythrum <- 
  ggplot(f_lythrum, aes(x = capped_season_length, y = f_mn)) +
  geom_ribbon(aes(ymax = f_max, ymin = f_min), alpha = 0.2, fill = lythrum_colour) +
  geom_line(linewidth = 0.7, colour = lythrum_colour) +
  annotate("text", label = TeX(r"($f(T) = \frac{a}{1 + e^{-b(T + c)}})"), 
           x= 180, y = 0.0103, size = 2.5) +
  annotate("text", label = "a = 0.02", x = 140, y = 0.0085, size = 2, hjust = 0) +
  annotate("text", label = "b = -0.01", x = 140, y = 0.0080, size = 2, hjust = 0) +
  annotate("text", label = "c = -144", x = 140, y = 0.0075, size = 2, hjust = 0) +
  scale_x_continuous(breaks = seq(130,210, 20)) +
  labs(x = expression(italic(T)), y = expression(italic(f))) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), 
        plot.background = element_blank())

#### 6b. Solidago
post_draws_f_solidago <- 
  as_draws_df(solidago_iwasa_cohen) %>% 
  as_tibble() %>% 
  dplyr::select(1:3) %>% 
  slice(rep(1:n(), each = 100)) %>% 
  mutate(capped_season_length = rep(seq(116,160, length.out = 100), 8000),
         f = b_a_Intercept/(1 + exp(-b_b_Intercept*(capped_season_length + b_c_Intercept))))

f_solidago <- post_draws_f_solidago %>% 
  group_by(capped_season_length) %>% 
  summarise(f_mn = mean(f),
            f_50l = quantile(f, prob = 0.25),
            f_50u = quantile(f, prob = 0.75),
            f_min = quantile(f, prob = 0.005),
            f_max = quantile(f, prob = 0.995))

f_T_solidago <- 
  ggplot(f_solidago, aes(x = capped_season_length, y = f_mn)) +
  geom_ribbon(aes(ymax = f_max, ymin = f_min), alpha = 0.2, fill = solidago_colour) +
  geom_line(linewidth = 0.7, colour = solidago_colour) +
  annotate("text", label = TeX(r"($f(T) = \frac{a}{1 + e^{-b(T + c)}})"), 
           x= 144, y = 0.0214, size = 2.5) +
  annotate("text", label = "a = 0.03", x = 120, y = 0.017, size = 2, hjust = 0) +
  annotate("text", label = "b = -0.04", x = 120, y = 0.0155, size = 2, hjust = 0) +
  annotate("text", label = "c = -162", x = 120, y = 0.014, size = 2, hjust = 0) +
  scale_x_continuous(breaks = seq(100,160, 20)) +
  scale_y_continuous(breaks = seq(0.0125, 0.025, 0.005)) +
  labs(x = expression(italic(T)), y = expression(italic(f))) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        plot.background = element_blank())

#__________________________________________________________________________________________
#### 7. Posterior smoothed predictions ####

## 99% parameter uncertainty

### 7a. Lythrum
## sigma 50% quantile
lythrum_sigma_50 <- as_draws_df(lythrum_iwasa_cohen) %>% 
  pull(sigma) %>% 
  quantile(., 0.5)

## getting posterior parameter values
pred_smooth_lythrum <- as_draws_df(lythrum_iwasa_cohen) %>% 
  dplyr::select(b_a_Intercept:b_c_Intercept) %>% 
  mutate(sim = 1:n()) %>% 
  slice(rep(1:n(), each = 100)) %>% 
  mutate(capped_season_length = rep(seq(134,206, length.out = 100), 8000),
         f = b_a_Intercept/(1 + exp(-b_b_Intercept*(capped_season_length + b_c_Intercept))),
         t_2 = capped_season_length - 
           (1/(b_a_Intercept/(1 + exp(-b_b_Intercept*(capped_season_length + b_c_Intercept)))))) %>% 
  group_by(capped_season_length) %>% 
  summarise(fit50 = quantile(t_2,  0.50),
            lwr = quantile(t_2, 0.005), upr = quantile(t_2, 0.995),
            lwr_full = lwr - lythrum_sigma_50,
            upr_full = upr + lythrum_sigma_50)

## mean +/- se of raw data for each population
lythrum_average <- lythrum_mdat %>% 
  group_by(pop) %>% 
  summarise(mn = mean(flowering_start),
            se = sd(flowering_start)/(sqrt(n())),
            capped_season_length = capped_season_length[1])

## fitted plot
fitted_lythrum <- 
  ggplot(lythrum_average, aes(x = capped_season_length, y = mn)) +
  geom_errorbar(aes(ymin = mn - se, ymax = mn + se, y = NULL), width = 0, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_ribbon(data = pred_smooth_lythrum, aes(ymax = upr, ymin = lwr, y = NULL), 
              alpha = 0.3, fill = lythrum_colour) +
  geom_line(data = pred_smooth_lythrum, stat = "identity", aes(y = fit50), 
            linewidth = 1, colour = lythrum_colour) +
  annotate("text", label = TeX(r"($t_2 = T - \frac{1}{f(T)})"), x = 193, y = 75.5) +
  labs(x = expression(paste("Growing Season Length (", italic("T"), ") [days]")), 
       y = expression(paste("Flowering time from start of season (", italic(t[italic("2")]), ")")), 
       tag = "a)") +
  coord_cartesian(ylim = c(40,78)) +
  scale_y_continuous(breaks = seq(0, 120, 10)) +
  scale_x_continuous(breaks = seq(100,220, by = 10)) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

## full plot
lythrum_plot_full <- ggdraw(fitted_lythrum) +
  draw_plot(f_T_lythrum, x = 0.18, y = 0.62, width = 0.42, height = 0.3) 

### 7a. Solidago
## sigma 50% quantile
solidago_sigma_50 <- as_draws_df(solidago_iwasa_cohen) %>% 
  pull(sigma) %>% 
  quantile(., 0.5)

## getting posterior average parameter values
pred_smooth_solidago <- as_draws_df(solidago_iwasa_cohen) %>% 
  dplyr::select(b_a_Intercept:b_c_Intercept) %>% 
  mutate(sim = 1:n()) %>% 
  slice(rep(1:n(), each = 100)) %>% 
  mutate(capped_season_length = rep(seq(116,160, length.out = 100), 8000),
         f = b_a_Intercept/(1 + exp(-b_b_Intercept*(capped_season_length + b_c_Intercept))),
         t_2 = capped_season_length - 
           (1/(b_a_Intercept/(1 + exp(-b_b_Intercept*(capped_season_length + b_c_Intercept)))))) %>% 
  group_by(capped_season_length) %>% 
  summarise(fit50 = quantile(t_2,  0.50),
            lwr = quantile(t_2, 0.005), upr = quantile(t_2, 0.995),
            lwr_full = lwr - solidago_sigma_50,
            upr_full = upr + solidago_sigma_50)

## mean +/- se of raw data for each population
solidago_average <- solidago_mdat %>% 
  group_by(pop) %>% 
  summarise(mn = mean(flowering_start),
            se = sd(flowering_start)/(sqrt(n())),
            capped_season_length = capped_season_length[1])

## fitted plot
fitted_solidago <- 
  ggplot(solidago_average, aes(x = capped_season_length, y = mn)) +
  geom_errorbar(aes(ymin = mn - se, ymax = mn + se, y = NULL), width = 0, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_ribbon(data = pred_smooth_solidago, aes(ymax = upr, ymin = lwr, y = NULL), 
              alpha = 0.3, fill = solidago_colour) +
  geom_line(data = pred_smooth_solidago, stat = "identity", aes(y = fit50), 
            linewidth = 1, colour = solidago_colour) +
  annotate("text", label = TeX(r"($t_2 = T - \frac{1}{f(T)})"), x = 152, y = 98) +
  labs(x = expression(paste("Growing Season Length (", italic("T"), ") [days]")), 
       y = expression(paste("Flowering time from start of season (", italic(t[italic("2")]), ")")), 
       tag = "b)") +
  coord_cartesian(ylim = c(70,100)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  scale_x_continuous(breaks = seq(100,170, by = 10)) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

## full plot
solidago_plot_full <- ggdraw(fitted_solidago) +
  draw_plot(f_T_solidago, x = 0.2, y = 0.62, width = 0.42, height = 0.3) 

#_______________________________________________________________________________
#### 8. Save plot ####

ggsave(lythrum_plot_full + solidago_plot_full, 
       filename = "output/predictions_20241104.jpeg",
       width = 22, height = 11, units = "cm", dpi = 700)

#_______________________________________________________________________________
#### 9. Parameter credible intervals for main text ####

as_draws_df(lythrum_iwasa_cohen) %>% 
  dplyr::select(b_a_Intercept:b_c_Intercept) %>% 
  mutate(sim = 1:n()) %>% 
  pivot_longer(-sim) %>% 
  group_by(name) %>% 
  summarise(mn50 = quantile(value, 0.5),
            mn = mean(value),
            lwr = quantile(value, 0.005),
            upr = quantile(value, 0.995))
  
as_draws_df(solidago_iwasa_cohen) %>% 
  dplyr::select(b_a_Intercept:b_c_Intercept) %>% 
  mutate(sim = 1:n()) %>% 
  pivot_longer(-sim) %>% 
  group_by(name) %>% 
  summarise(mn50 = quantile(value, 0.5),
            mn = mean(value),
            lwr = quantile(value, 0.005),
            upr = quantile(value, 0.995))
