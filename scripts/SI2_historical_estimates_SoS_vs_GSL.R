#####################################################
##                                                 ##
##     In situ historical flowering estimates      ##
##                                                 ##
##                                                 ##
##                                                 ##
##               JP - 13th June 2024               ##
##                                                 ##
#####################################################

rm(list = ls())
options(width = 100)

library(tidyverse)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

# Load historical temp data
lythrum_pops_temps <- read.csv("data/SMHI_pthbv_t_1961_2022_daily_4326-lythrumpops.csv")
solidago_pops_temps <- read.csv("data/SMHI_pthbv_t_1961_2022_daily_4326-Solidagopops.csv")

lythrum_pops <- c(1:15) # no. lythrum pops
solidago_pops <- c(1:24)
solidago_boreal_pops <- c(16, 24, 3, 9, 6, 19, 13, 10) # pop id's that correspond to boreal subsample, consistent with main analysis

year_range <- c(1961, 2022)

#-----------------------------------------------------------------------------#
# Delineating by consecutive thermal boundaries (>threshold_temp for >5 days) #
#-----------------------------------------------------------------------------#

threshold_temp <- 5 
cap_length <- 5 # no. days of consecutive threshold temp to surpass

season_consec_caps.f <- function(dat, pops, minyear, maxyear, threshold_temp){
  names(dat) <- NULL
  dat <- as.vector(dat[,1])
  df <- colsplit(dat, ";", c("date", c(paste0("pop", pops))))
  df <- separate(df, col = date, into = c("year", "month", "day"), sep = "\\-")
  
  ## Subsetting years in the dataset
  df <- subset(df, year >= minyear & year <= maxyear)
  
  ## EXTRACTION of growing season bounds  & GSL (growing season length)
  
  years <- unique(df$year)
  
  spring_start <- matrix(NA, ncol = 1 + length(pops), nrow = length(years)) # empty results mat for start of spring -- "1 +" because first column is years
  autumn_start <- matrix(NA, ncol = 1 + length(pops), nrow = length(years)) # empty results mat for start of spring -- "1 +" because first column is years
  gsl <- matrix(NA, ncol = 1 + length(pops), nrow = length(years)) # same as above, but for growing season length (autumn start - spring start)
  
  for (i in 1:length(years)){
    
    annum_temp <- subset(df, df$year == years[i]) # subset df into year
    
    spring_start[i,1] <- years[i]
    autumn_start[i,1] <- years[i]
    
    gsl[i,1] <- years[i]
    
    cap_length <- cap_length # number of consecutive days we want surpassing threshold temp
    
    # function to find start of spring
    find_spring_cap <- function(temp) { # first date in the spring when consecutive N (=5) days are over X (>5) degrees C
      spring_caps <- c()
      for (i in 1:(length(temp)-cap_length)){
        spring_caps[i] <- sum(temp[i:(i+cap_length-1)] > threshold_temp)
      }
      ss <- min(which(spring_caps == cap_length))
      return(ss)
    }
    find_autumn_cap <- function(temp) { # first time in the year 
      autumn_caps <- c()
      for (i in 1:(length(temp)-cap_length)){
        autumn_caps[i] <- sum(temp[i:(i+cap_length-1)] < threshold_temp)
      }
      as <- which(autumn_caps[(182+1):length(autumn_caps)] == cap_length)[1] + 182 # first date after middle of year (day 182) when consecutive N (=5) days are under X (<5) degrees C
      return(as)
    }
    
    for (p in 1:length(pops)){
      spring_start[i,1+p] <- find_spring_cap(annum_temp[,3+p])
      autumn_start[i,1+p] <- find_autumn_cap(annum_temp[,3+p])
      gsl[i,1+p] <- find_autumn_cap(annum_temp[,3+p]) - find_spring_cap(annum_temp[,3+p])
      
    }
  }
  spring_start <- as.data.frame(spring_start); colnames(spring_start) <- c("year", c(pops))
  autumn_start <- as.data.frame(autumn_start); colnames(autumn_start) <- c("year", c(pops))
  gsl <- as.data.frame(gsl); colnames(gsl) <- c("year", c(pops))
  data.frame(apply(spring_start, 2, function(x) as.numeric(as.numeric(x))))
  
  return(list(spring_start, autumn_start, gsl))
}

lythrum_seasoncaps_predict.list <- season_consec_caps.f(lythrum_pops_temps, lythrum_pops, min(year_range), max(year_range), threshold_temp)

solidago_seasoncaps_predict.list <- season_consec_caps.f(solidago_pops_temps, solidago_pops, min(year_range), max(year_range), threshold_temp)

#-----------------------------------------------------------------------------#
#         Organize yearly SoS and GSL calculations for both species           #
#-----------------------------------------------------------------------------#

lythrum_sos <- lythrum_seasoncaps_predict.list[[1]]
lythrum_gsl <- lythrum_seasoncaps_predict.list[[3]]

solidago_sos <- solidago_seasoncaps_predict.list[[1]]
solidago_sos <- solidago_sos[,c(1, solidago_boreal_pops + 1)] # +1 because we are keeping in the "year" column
solidago_gsl <- solidago_seasoncaps_predict.list[[3]]
solidago_gsl <- solidago_gsl[,c(1, solidago_boreal_pops + 1)]

# long
lythrum_sos_long <- lythrum_sos %>%
  pivot_longer(cols = -year, names_to = "population", values_to = "sos")
lythrum_gsl_long <- lythrum_gsl %>%
  pivot_longer(cols = -year, names_to = "population", values_to = "gsl")
solidago_sos_long <- solidago_sos %>%
  pivot_longer(cols = -year, names_to = "population", values_to = "sos")
solidago_gsl_long <- solidago_gsl %>%
  pivot_longer(cols = -year, names_to = "population", values_to = "gsl")

# merge

lythrum_merged <- merge(lythrum_sos_long, lythrum_gsl_long, by = c("year", "population"))
lythrum_merged$year <- as.numeric(lythrum_merged$year)
lythrum_merged$sos <- as.numeric(lythrum_merged$sos)
lythrum_merged$gsl <- as.numeric(lythrum_merged$gsl)

solidago_merged <- merge(solidago_sos_long, solidago_gsl_long, by = c("year", "population"))
solidago_merged$year <- as.numeric(solidago_merged$year)
solidago_merged$sos <- as.numeric(solidago_merged$sos)
solidago_merged$gsl <- as.numeric(solidago_merged$gsl)

#-----------------------------------------------------------------------------#
#                                 Plot                                        #
#-----------------------------------------------------------------------------#
lythrum_palette <- c("#8C4681", "#A6449F", "#172601", "#3E5902", "#567639")
solidago_palette <- c("#0e2d10", "#214d08", "#90980c", "#d2c200", "#fce708")

##--SoS vs. GSL correlation over the years--##

# calculate linear regressions to plot
lythrum_populations <- unique(lythrum_merged$population)
lythrum_slopes <- sapply(lythrum_populations, function(lythrum_populations) {
  lm(gsl ~ sos, data = lythrum_merged[lythrum_merged$population == lythrum_populations, ])$coefficients[2]
})
lythrum_slope_data <- data.frame(population = lythrum_populations, slope = lythrum_slopes)
lythrum_slope_labels <- lythrum_slope_data
lythrum_slope_labels$label <- sprintf("Slope: %.2f", lythrum_slope_labels$slope)

solidago_populations <- unique(solidago_merged$population)
solidago_slopes <- sapply(solidago_populations, function(solidago_populations) {
  lm(gsl ~ sos, data = solidago_merged[solidago_merged$population == solidago_populations, ])$coefficients[2]
})
solidago_slope_data <- data.frame(population = solidago_populations, slope = solidago_slopes)
solidago_slope_labels <- solidago_slope_data
solidago_slope_labels$label <- sprintf("Slope: %.2f", solidago_slope_labels$slope)

lythrum_sos_vs_gsl.p <- ggplot(lythrum_merged, aes(x = sos, y = gsl, color = year)) +
  geom_point()+
  scale_colour_gradient(low = lythrum_palette[3], high = lythrum_palette[2]) +
  geom_smooth(method = "lm", se = FALSE, aes(group = population), col = "red") +
  facet_wrap(~population, nrow = 3, ncol = 5) +
  labs(x = "Start of spring",
       y = "Growing season length",
       colour = "Year",
       shape = "Population",
       tag = expression(paste("a) ", italic(Lythrum)))) +
  geom_text(data = lythrum_slope_labels, aes(label = label, x = Inf, y = Inf), hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
  theme_bw() +
  theme(strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank()) 
  
solidago_sos_vs_gsl.p <- ggplot(solidago_merged, aes(x = sos, y = gsl, color = year)) +
  geom_point()+
  scale_colour_gradient(low = solidago_palette[1], high = solidago_palette[4]) +
  geom_smooth(method = "lm", se = FALSE, aes(group = population), col = "red") +
  facet_wrap(~population, nrow = 2, ncol =4) +
  labs(x = "Start of spring",
       y = "Growing season length",
       colour = "Year",
       shape = "Population",
       tag = expression(paste("b) ", italic(Solidago)))) +
  geom_text(data = solidago_slope_labels, aes(label = label, x = Inf, y = Inf), hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +  theme_bw() +
  theme_bw() +
  theme(strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank()) 

sos_vs_gsl.p <- grid.arrange(lythrum_sos_vs_gsl.p, solidago_sos_vs_gsl.p, ncol = 1)

ggsave(sos_vs_gsl.p,
       filename = "output/sos_vs_gsl_20241010.jpeg",
       width = 22, height = 17, units = "cm", dpi = 700)

##-- Fall start vs. GSL ? --##
# extract start of autumn estimates from same algorithm
lythrum_soa <- lythrum_seasoncaps_predict.list[[2]] 
solidago_soa <- solidago_seasoncaps_predict.list[[2]]
solidago_soa <- solidago_soa[,c(1, solidago_boreal_pops + 1)] 

# long & append to long-merged df
lythrum_soa_long <- lythrum_soa %>%
  pivot_longer(cols = -year, names_to = "population", values_to = "soa")
solidago_soa_long <- solidago_soa %>%
  pivot_longer(cols = -year, names_to = "population", values_to = "soa")

lythrum_merged$soa <- as.numeric(lythrum_soa_long$soa)
solidago_merged$soa <- as.numeric(solidago_soa_long$soa)

# linear regressions to plot
lythrum_autumn_slopes <- sapply(lythrum_populations, function(lythrum_populations) {
  lm(gsl ~ soa, data = lythrum_merged[lythrum_merged$population == lythrum_populations, ])$coefficients[2]
})
lythrum_autumn_slope_data <- data.frame(population = lythrum_populations, slope = lythrum_autumn_slopes)
lythrum_autumn_slope_labels <- lythrum_autumn_slope_data
lythrum_autumn_slope_labels$label <- sprintf("Slope: %.2f", lythrum_autumn_slope_labels$slope)


solidago_autumn_slopes <- sapply(solidago_populations, function(solidago_populations) {
  lm(gsl ~ soa, data = solidago_merged[solidago_merged$population == solidago_populations, ])$coefficients[2]
})
solidago_autumn_slope_data <- data.frame(population = solidago_populations, slope = solidago_autumn_slopes)
solidago_autumn_slope_labels <- solidago_autumn_slope_data
solidago_autumn_slope_labels$label <- sprintf("Slope: %.2f", solidago_autumn_slope_labels$slope)

# lythrum_soa_vs_gsl.p <- ggplot(lythrum_merged, aes(x = soa, y = gsl, color = year)) +
#   geom_point()+
#   scale_colour_gradient(low = lythrum_palette[3], high = lythrum_palette[2]) +
#   geom_smooth(method = "lm", se = FALSE, aes(group = population), col = "red") +
#   facet_wrap(~population, nrow = 3, ncol = 5) +
#   labs(x = "Start of Autumn",
#        y = "Growing Season Length",
#        colour = "Year",
#        shape = "Population",
#        tag = "a) Lythrum") +
#   geom_text(data = lythrum_autumn_slope_labels, aes(label = label, x = Inf, y = Inf), hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +  theme_bw() +
#   theme_bw() +
#   theme(strip.text = element_blank(),
#         panel.grid = element_blank(),
#         axis.ticks = element_blank()) 

solidago_soa_vs_gsl.p <- ggplot(solidago_merged, aes(x = soa, y = gsl, color = year)) +
  geom_point()+
  scale_colour_gradient(low = solidago_palette[1], high = solidago_palette[4]) +
  geom_smooth(method = "lm", se = FALSE, aes(group = population), col = "red") +
  facet_wrap(~population, nrow = 2, ncol =4) +
  labs(x = "Start of Autumn",
       y = "Growing Season Length",
       colour = "Year",
       shape = "Population",
       tag = "b) Solidago") +
  geom_text(data = solidago_autumn_slope_labels, aes(label = label, x = Inf, y = Inf), hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +  theme_bw() +
  theme(strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank()) 

soa_vs_gsl.p <- grid.arrange(lythrum_soa_vs_gsl.p, solidago_soa_vs_gsl.p, ncol = 1)

ggsave(soa_vs_gsl.p,
       filename = "output/soa_vs_gsl_20240617.jpeg",
       width = 22, height = 17, units = "cm", dpi = 700)


##-- Historical timeseries of SoS and GSL separately --##


lythrum_spring_historicaltrends.p <- ggplot(data = lythrum_merged, aes(x = year, y = sos, colour = population)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, aes(group = population), linewidth = 0.3) +
  ylim(min(lythrum_merged$sos),max(c(lythrum_merged$gsl))) +
  labs(x = "Year",
       y = "Start of spring",
       shape = "Population") +
  theme_bw() +
  theme(strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 
# lythrum_autumn_historicaltrends.p <- ggplot(data = lythrum_merged, aes(x = year, y = soa, colour = population)) +
#   geom_point(size = 0.5) +
#   geom_smooth(method = "lm", se = FALSE, aes(group = population), linewidth = 0.3) +
#   ylim(min(lythrum_merged$sos),max(c(lythrum_merged$soa))) +
#   labs(x = "Year",
#        y = "Start of Autumn",
#        shape = "Population") +
#   theme_bw() +
#   theme(strip.text = element_blank(),
#         panel.grid = element_blank(),
#         axis.ticks = element_blank(),
#         legend.position = "none") 
lythrum_gsl_historicaltrends.p <- ggplot(data = lythrum_merged, aes(x = year, y = gsl, colour = population)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, aes(group = population), linewidth = 0.3) +
  ylim(min(lythrum_merged$sos),max(c(lythrum_merged$gsl))) +
  labs(x = "Year",
       y = "Growing season length",
       shape = "Population") +
  theme_bw() +
  theme(strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 

lythrum_trends.p <- grid.arrange(lythrum_spring_historicaltrends.p, lythrum_gsl_historicaltrends.p,
                                 ncol = 2)

solidago_spring_historicaltrends.p <- ggplot(data = solidago_merged, aes(x = year, y = sos, colour = population)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, aes(group = population), linewidth = 0.3) +
  ylim(min(solidago_merged$sos),max(c(solidago_merged$gsl))) +
  labs(x = "Year",
       y = "Start of spring",
       shape = "Population") +
  theme_bw() +
  theme(strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 
# solidago_autumn_historicaltrends.p <- ggplot(data = solidago_merged, aes(x = year, y = soa, colour = population)) +
#   geom_point(size = 0.5) +
#   geom_smooth(method = "lm", se = FALSE, aes(group = population), linewidth = 0.3) +
#   ylim(min(solidago_merged$sos),max(c(solidago_merged$soa))) +
#   labs(x = "Year",
#        y = "Start of Autumn",
#        shape = "Population") +
#   theme_bw() +
#   theme(strip.text = element_blank(),
#         panel.grid = element_blank(),
#         axis.ticks = element_blank(),
#         legend.position = "none") 
solidago_gsl_historicaltrends.p <- ggplot(data = solidago_merged, aes(x = year, y = gsl, colour = population)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, aes(group = population), linewidth = 0.3) +
  ylim(min(solidago_merged$sos),max(c(solidago_merged$gsl))) +
  labs(x = "Year",
       y = "Growing season length",
       shape = "Population") +
  theme_bw() +
  theme(strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 

solidago_trends.p <- grid.arrange(solidago_spring_historicaltrends.p, solidago_gsl_historicaltrends.p,
                                 ncol = 2)

historical_trends.p <- grid.arrange(arrangeGrob(lythrum_trends.p, top = textGrob(expression(paste("a) ", italic(Lythrum))), x = 0, hjust = 0)), 
                                    arrangeGrob(solidago_trends.p, top = textGrob(expression(paste("b) ", italic(Solidago))), x= 0, hjust = 0)))

ggsave(historical_trends.p,
       filename = "output/historical_trends_20241010.jpeg",
       width = 15, height = 12, units = "cm", dpi = 700)

