#####################################################
##                                                 ##
##           In situ flowering estimates           ##
##                                                 ##
##     Given in situ start-of-spring and "t_2"     ##
##                                                 ##
##               JP - 12th Sept 2024               ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)
library(patchwork)
library(grid)
library(gridExtra)

#_______________________________________________________________________________
## Load processed data, and get mean SOS (start of spring) and mean gsl for each pop in situ, for the chosen year range

# Lythrum
load("data/lythrum_seasoncaps_predict.list.RData")
lythrum_sos <- lythrum_seasoncaps_predict.list[[1]][,-1]
lythrum_gsl <- lythrum_seasoncaps_predict.list[[3]][,-1]
lythrum_summary <- tibble(pop = as.numeric(colnames(lythrum_sos)),
                           mean_sos =
                             apply(lythrum_sos, 2, function(x){mean(as.numeric(x))}),
                           mean_gsl = 
                             apply(lythrum_gsl, 2, function(x){mean(as.numeric(x))})) 
load("data/umea_predict.list.RData")
umea_cg_sos <- as.numeric(umea_predict.list[[1]][which(umea_predict.list[[1]][,1] == "1999"), 2]) # SOS in experiment year (1999)

# Solidago
load("data/solidago_seasoncaps_predict.list.RData")
solidago_boreal <- c(16, 24, 3, 9, 6, 19, 13, 10) # pop id's that correspond to boreal subsample, consistent with main analysis
solidago_sos <- solidago_seasoncaps_predict.list[[1]][,-1]
solidago_sos <- solidago_sos[,solidago_boreal]
solidago_gsl <- solidago_seasoncaps_predict.list[[3]][,-1]
solidago_gsl <- solidago_gsl[,solidago_boreal]
solidago_summary <- tibble(pop = as.numeric(colnames(solidago_sos)),
                          mean_sos =
                            apply(solidago_sos, 2, function(x){mean(as.numeric(x))}),
                          mean_gsl = 
                            apply(solidago_gsl, 2, function(x){mean(as.numeric(x))})) 
load("data/uppsala_predict.list.RData")
uppsala_cg_sos <- as.numeric(uppsala_predict.list[[1]][which(uppsala_predict.list[[1]][,1] == "2005"), 2]) # SOS in experiment year (2005)


## Colours
lythrum_palette <- c("#8C4681", "#A6449F", "#172601", "#3E5902", "#567639")
lythrum_colour <- lythrum_palette[2]

solidago_palette <- c("#0e2d10", "#214d08", "#90980c", "#d2c200", "#fce708")
solidago_colour <- solidago_palette[4]


#_______________________________________________________________________________
## Fit linear regression on SOS ~ GSL, and add this y=ax+b function to fitted t_2(T) nonlinear function (params obtained from fit)

#---Lythrum---#

lythrum_gsl_range <- seq(min(lythrum_summary$mean_gsl), max(lythrum_summary$mean_gsl), length.out = 100) # range of GSL for Lythrum pops

# IN SITU SoS -- linear regression
lythrum_lm <- lm(mean_sos ~ mean_gsl, data = lythrum_summary) # linear regression of SOS to GSL
lythrum_sos_coeff_a <- coef(lythrum_lm)[1] #intercept
lythrum_sos_coeff_b <- coef(lythrum_lm)[2] #slope

# CG flowering -- add first flower day of year to all t2 values
lythrum_first_flower <- 180 # first flower, to which all data were normalized in main analysis was day 180 (June 29)
lythrum_flowering.f <- function(x) {
  return(x - (1/(0.021/(1 + exp(0.01 * (x - 144))))) + lythrum_first_flower)
  } 

# IN SITU flowering -- compute diff's between CG SoS (Umea) and IN SITU SoS
lythrum_sos_diffs <- function(x) {
  return(lythrum_sos_coeff_a + lythrum_sos_coeff_b*x - umea_cg_sos)
}

# -- then add SoS diffs to fitted flowering function
lythrum_insitu.f <- function(x) {
  return(lythrum_sos_diffs(x) + lythrum_flowering.f(x))
}

# final curve of estimated in situ flowering times
lythrum_insitu_est <- lythrum_insitu.f(lythrum_gsl_range)


#---Solidago---#

solidago_gsl_range <- seq(min(solidago_summary$mean_gsl), max(solidago_summary$mean_gsl), length.out = 100) # range of GSL for Solidago pops

# IN SITU SoS -- linear regression
solidago_lm <- lm(mean_sos ~ mean_gsl, data = solidago_summary) # linear regression of SOS to GSL
solidago_sos_coeff_a <- coef(solidago_lm)[1] #intercept
solidago_sos_coeff_b <- coef(solidago_lm)[2] #slope

# CG flowering -- add first flower day of year to all t2 values
solidago_first_flower <- 158 # first flower, to which data were normalized in main analysis, was day 158
solidago_flowering.f <- function(x) {
  return(x - (1/(0.026/(1 + exp(0.045 * (x - 161))))) + solidago_first_flower)
} 

# IN SITU flowering -- compute diff's between CG SoS (Uppsala) and IN SITU SoS
solidago_sos_diffs <- function(x) {
  return(solidago_sos_coeff_a + solidago_sos_coeff_b*x - uppsala_cg_sos)
}
# -- then add SoS diffs to fitted flowering function
solidago_insitu.f <- function(x) {
  return(solidago_sos_diffs(x) + solidago_flowering.f(x))
}

# final curve of estimated in situ flowering times

solidago_insitu_est <- solidago_insitu.f(solidago_gsl_range)


#_______________________________________________________________________________
## inset plots - SoS variation

# Combined aesthetic for the legend
combined_aesthetic <- function(col, linetype) {
  paste(col, linetype, sep = "_")
}

# Custom labels for the legend
custom_labels_sos <- c(
  "sos_cg" = "Common garden",
  "sos_insitu" = expression(italic("in situ"))
)

lythrum_insitucompar.p <- ggplot() +
  geom_line(aes(x = c(-Inf, Inf), y = umea_cg_sos, 
                group = combined_aesthetic("sos", "cg"), 
                col = combined_aesthetic("sos", "cg"), 
                linetype = combined_aesthetic("sos", "cg")), 
                linewidth = 1.2, alpha = 0.6) + # horizontal line showing common SOS starting line in Umea cg experiment
  
  annotate(geom = "text", label = expression(paste(italic("Ume" * "\u00E5 "), "(63.5"*~degree*"N)")), x = 172, y = umea_cg_sos + 4, size = 5, col = "gray55") +
  
  geom_point(data = lythrum_summary, 
             aes(x = as.numeric(mean_gsl), y = as.numeric(mean_sos), 
                 col = combined_aesthetic("sos", "insitu")), 
             size = 3, shape = 18, alpha = 1) + # SOS against GSL gradient
  geom_segment(aes(x = min(lythrum_gsl_range), 
                   xend = max(lythrum_gsl_range), 
                   y = lythrum_sos_coeff_a + lythrum_sos_coeff_b * min(lythrum_gsl_range), 
                   yend = lythrum_sos_coeff_a + lythrum_sos_coeff_b * max(lythrum_gsl_range),
                   col = combined_aesthetic("sos", "insitu"), 
                   linetype = combined_aesthetic("sos", "insitu")),
               linewidth = 1.2, alpha = 1) + # linear fit of SOS~GSL
  scale_color_manual(values = c(
    "sos_cg" = "black",
    "sos_insitu" = "black"), 
    labels = custom_labels_sos) +
  scale_linetype_manual(values = c(
    "sos_cg" = 11,
    "sos_insitu" = "solid"), 
    labels = custom_labels_sos) +
  guides(col = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) + 

  geom_point(aes(x = min(lythrum_gsl_range), y = lythrum_summary$mean_sos[which.min(lythrum_summary$mean_gsl)]), shape = 1, size = 5, col = "gray55") +
  geom_point(aes(x = max(lythrum_gsl_range), y = lythrum_summary$mean_sos[which.max(lythrum_summary$mean_gsl)]), shape = 1, size = 5, col = "gray55") +
  
  annotate(geom = "text", label = expression("66.1"*~degree*"N"), x = min(lythrum_gsl_range) + 12, y = lythrum_summary$mean_sos[which.min(lythrum_summary$mean_gsl)] - 17, size = 5, col = "gray55") +
  annotate(geom = "text", label = expression("57.4"*~degree*"N"), x = max(lythrum_gsl_range) - 28, y = lythrum_summary$mean_sos[which.max(lythrum_summary$mean_gsl)] + 2, size = 5, col = "gray55") +
  
  geom_curve(aes(x = min(lythrum_gsl_range), y = lythrum_summary$mean_sos[which.min(lythrum_summary$mean_gsl)], 
                 xend = min(lythrum_gsl_range) + 7, yend = lythrum_summary$mean_sos[which.min(lythrum_summary$mean_gsl)] - 13),
             arrow = arrow(length = unit(0.1, "inch")), size = 0.6,
             color = "gray50", curvature = 0.1) +
  geom_curve(aes(x = max(lythrum_gsl_range), y = lythrum_summary$mean_sos[which.max(lythrum_summary$mean_gsl)], 
                 xend = max(lythrum_gsl_range) - 14, yend = lythrum_summary$mean_sos[which.max(lythrum_summary$mean_gsl)] + 1),
             arrow = arrow(length = unit(0.1, "inch")), size = 0.6,
             color = "gray50", curvature = -0.1) +
  
  labs(x = expression(italic("T")),
       y = "Start of spring") +
  
  ylim(lythrum_summary$mean_sos[which.min(lythrum_summary$mean_sos)] - 5, lythrum_summary$mean_sos[which.max(lythrum_summary$mean_sos)] + 7) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        #legend.text = element_text(size = 11),
        legend.position = "none",
        plot.margin= unit(c(0, 0, 0, 0), "pt"))


solidago_insitucompar.p <- ggplot() +
  geom_line(aes(x = c(-Inf, Inf), y = uppsala_cg_sos, 
                group = combined_aesthetic("sos", "cg"), 
                col = combined_aesthetic("sos", "cg"), 
                linetype = combined_aesthetic("sos", "cg")), 
            linewidth = 1.2, alpha = 0.6) + # horizontal line showing common SOS starting line in Uppsala cg experiment
  
  annotate(geom = "text", label = expression(paste(italic("Uppsala "), "(59.8"*~degree*"N)")), x = 135, y = uppsala_cg_sos + 5, size = 5, col = "gray55") +

  geom_point(data = solidago_summary, 
             aes(x = as.numeric(mean_gsl), y = as.numeric(mean_sos), 
                 col = combined_aesthetic("sos", "insitu")), 
             size = 3, shape = 18, alpha = 1) + # SOS against GSL gradient
  geom_segment(aes(x = min(solidago_gsl_range), 
                   xend = max(solidago_gsl_range), 
                   y = solidago_sos_coeff_a + solidago_sos_coeff_b * min(solidago_gsl_range), 
                   yend = solidago_sos_coeff_a + solidago_sos_coeff_b * max(solidago_gsl_range),
                   col = combined_aesthetic("sos", "insitu"), 
                   linetype = combined_aesthetic("sos", "insitu")),
               linewidth = 1.2, alpha = 1) + # linear fit of SOS~GSL
  scale_color_manual(values = c(
    "sos_cg" = "black",
    "sos_insitu" = "black"), 
    labels = custom_labels_sos) +
  scale_linetype_manual(values = c(
    "sos_cg" = 11,
    "sos_insitu" = "solid"), 
    labels = custom_labels_sos) +
  guides(col = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) + 
  
  geom_point(aes(x = min(solidago_gsl_range), y = solidago_summary$mean_sos[which.min(solidago_summary$mean_gsl)]), shape = 1, size = 5, col = "gray55") +
  geom_point(aes(x = max(solidago_gsl_range), y = solidago_summary$mean_sos[which.max(solidago_summary$mean_gsl)]), shape = 1, size = 5, col = "gray55") +
  
  annotate(geom = "text", label = expression("67.9"*~degree*"N"), x = min(solidago_gsl_range) + 13, y = solidago_summary$mean_sos[which.min(solidago_summary$mean_gsl)] - 25, size = 5, col = "gray55") +
  annotate(geom = "text", label = expression("61.3"*~degree*"N"), x = max(solidago_gsl_range) - 10, y = solidago_summary$mean_sos[which.max(solidago_summary$mean_gsl)] - 17, size = 5, col = "gray55") +
  
  geom_curve(aes(x = min(solidago_gsl_range), y = solidago_summary$mean_sos[which.min(solidago_summary$mean_gsl)], 
                 xend = min(solidago_gsl_range) + 5, yend = solidago_summary$mean_sos[which.min(solidago_summary$mean_gsl)] - 23),
             arrow = arrow(length = unit(0.1, "inch")), size = 0.6,
             color = "gray55", curvature = 0.3) +
  geom_curve(aes(x = max(solidago_gsl_range), y = solidago_summary$mean_sos[which.max(solidago_summary$mean_gsl)], 
                 xend = max(solidago_gsl_range) - 2, yend = solidago_summary$mean_sos[which.max(solidago_summary$mean_gsl)] - 14),
             arrow = arrow(length = unit(0.1, "inch")), size = 0.6,
             color = "gray55", curvature = -0.3) +
  labs(x = expression(italic("T")),
       y = "Start of spring") +
  ylim(uppsala_cg_sos - 5, solidago_summary$mean_sos[which.max(solidago_summary$mean_sos)] + 3) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        #legend.text = element_text(size = 11),
        legend.position = "none",
        plot.margin= unit(c(0, 0, 0, 0), "pt"))


#_______________________________________________________________________________
# in situ flowering predictions, and extrapolations with +20 days of GSL expansion

# Custom labels for the legend
custom_labels_flowering <- c(
  "flwr_cg" = "Common garden",
  "flwr_insitu" = expression(italic("in situ")),
  "forecast" = "Forecast"
)

hypothetical_lythrum_gsl_expand <-  seq(max(lythrum_summary$mean_gsl), max(lythrum_summary$mean_gsl) + 20, length.out = 100)
hypothetical_lythrum_insitu_est <- lythrum_insitu.f(hypothetical_lythrum_gsl_expand)

hypothetical_lythrum.p <- ggplot() +
  geom_line(aes(x = lythrum_gsl_range, y = lythrum_insitu_est, 
                group = combined_aesthetic("flwr", "insitu"), 
                col = combined_aesthetic("flwr", "insitu"), 
                linetype = combined_aesthetic("flwr", "insitu")),
            linewidth = 5) + # in situ flowering estimate curve
  geom_line(aes(x = hypothetical_lythrum_gsl_expand, y = hypothetical_lythrum_insitu_est,
                group = "forecast",
                col = "forecast",
                linetype = "forecast"),
            linewidth = 5,
            alpha = 0.4) +
  geom_line(aes(x = lythrum_gsl_range, y = lythrum_flowering.f(lythrum_gsl_range), 
                group = combined_aesthetic("flwr", "cg"), 
                col = combined_aesthetic("flwr", "cg"), 
                linetype = combined_aesthetic("flwr", "cg")),
            linewidth = 2.5) + # CG flowering time (Iwasa-Cohen fitted curve)
  annotate(geom = "text", label = round(max(hypothetical_lythrum_insitu_est)),
           x = min(hypothetical_lythrum_gsl_expand),
           y = max(hypothetical_lythrum_insitu_est) + 6,
           size = 6) +
  annotate(geom = "text", label = round(min(hypothetical_lythrum_insitu_est)),
           x = max(hypothetical_lythrum_gsl_expand),
           y = min(hypothetical_lythrum_insitu_est) - 4,
           size = 6)+
  
  scale_color_manual(values = c(
    "flwr_cg" = lythrum_colour,
    "flwr_insitu" = lythrum_colour,
    "forecast" = lythrum_colour
  ), labels = custom_labels_flowering) +
  scale_linetype_manual(values = c(
    "flwr_cg" = 11,
    "flwr_insitu" = "solid",
    "forecast" = "solid"
  ), labels = custom_labels_flowering) +
  guides(col = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) + 
  
  ylim(160, max(lythrum_flowering.f(lythrum_gsl_range))) +
  labs(x = expression(paste("Growing Season Length (", italic("T"), ") [days]")),
       y = "Optimal flowering day of year",
       tag = "a)")+
  theme_bw(base_size = 12) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, "cm"),
        panel.grid = element_blank(),
        plot.margin= unit(c(0, 0, 0, 0), "pt"))

full_hypothetical_lythrum <- ggdraw(hypothetical_lythrum.p) +
  draw_plot(lythrum_insitucompar.p, x = 0.145, y = 0.12, width = 0.35, height = 0.45)

hypothetical_solidago_gsl_expand <-  seq(max(solidago_summary$mean_gsl), max(solidago_summary$mean_gsl) + 20, length.out = 100)
hypothetical_solidago_insitu_est <- solidago_insitu.f(hypothetical_solidago_gsl_expand)

hypothetical_solidago.p <- ggplot() +
  geom_line(aes(x = solidago_gsl_range, y = solidago_insitu_est, 
                group = combined_aesthetic("flwr", "insitu"), 
                col = combined_aesthetic("flwr", "insitu"), 
                linetype = combined_aesthetic("flwr", "insitu")),
            linewidth = 5) + # in situ flowering estimate curve
  geom_line(aes(x = hypothetical_solidago_gsl_expand, y = hypothetical_solidago_insitu_est,
                group = "forecast",
                col = "forecast",
                linetype = "forecast"),
            linewidth = 5,
            alpha = 0.4) +
  geom_line(aes(x = solidago_gsl_range, y = solidago_flowering.f(solidago_gsl_range), 
                group = combined_aesthetic("flwr", "cg"), 
                col = combined_aesthetic("flwr", "cg"), 
                linetype = combined_aesthetic("flwr", "cg")),
            linewidth = 2.5) + # CG flowering time (Iwasa-Cohen fitted curve)
  annotate(geom = "text", label = round(max(hypothetical_solidago_insitu_est)),
           x = min(hypothetical_solidago_gsl_expand),
           y = max(hypothetical_solidago_insitu_est) + 10,
           size = 6) +
  annotate(geom = "text", label = round(min(hypothetical_solidago_insitu_est)),
           x = max(hypothetical_solidago_gsl_expand),
           y = min(hypothetical_solidago_insitu_est) - 8,
           size = 6) +
  
  scale_color_manual(values = c(
    "flwr_cg" = solidago_colour,
    "flwr_insitu" = solidago_colour,
    "forecast" = solidago_colour
  ), labels = custom_labels_flowering) +
  scale_linetype_manual(values = c(
    "flwr_cg" = 11,
    "flwr_insitu" = "solid",
    "forecast" = "solid"
  ), labels = custom_labels_flowering) +
  guides(col = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) + 
  
  ylim(130, 290) +
  labs(x = expression(paste("Growing Season Length (", italic("T"), ") [days]")),
       y = "Optimal flowering day of year",
       tag = "b)") +
  theme_bw(base_size = 12) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.height = unit(1, "cm"),
        panel.grid = element_blank(),
        plot.margin= unit(c(0, 0, 0, 0), "pt"))

full_hypothetical_solidago <- ggdraw(hypothetical_solidago.p) +
  draw_plot(solidago_insitucompar.p, x = 0.145, y = 0.12, width = 0.35, height = 0.45)

ggsave(full_hypothetical_lythrum + full_hypothetical_solidago,
       filename = "output/hypothetical_expansion_combined_20241105.jpeg",
       width = 40, height = 14, units = "cm", dpi = 700)


       