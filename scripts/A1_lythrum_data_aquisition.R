#####################################################
##                                                 ##
##         Lythrum Optimal Flowering Time          ##
##                                                 ##
##               Data acquisition                  ##
##                                                 ##
##            JP & JJ - 16th July 2024             ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(reshape2)

#_______________________________________________________________________________
#### 1. Population coordinates ####

lythrum_coords <- read_csv("data/PopCoordinatesGatheredMay2023.csv") 
  
## Coordinates as decimals
lythrum_coordinates <- lythrum_coords %>% 
  separate(`WGS84 N`, paste("lat",c("d","m","s"), sep="_") ) %>%
  separate(`WGS84 E`, paste("long",c("d","m","s"), sep="_" ) ) %>%
  dplyr::select(pop = 1, 8:13) %>% 
  mutate(across(2:7, as.numeric)) %>%
  mutate(lat_dec=lat_d + lat_m/60 + lat_s/60^2,
         long_dec=long_d + long_m/60 + long_s/60^2)

pops <- lythrum_coordinates$pop
latitudes_pops <- lythrum_coordinates$lat_dec

#_______________________________________________________________________________
#### 2. Temperature data ####

pops_temps <- read.csv("data/SMHI_pthbv_t_1961_2022_daily_4326-lythrumpops.csv")
umea_cg_temps <- read.csv("data/umea_commongarden_1999_2001_daily_4326.csv")
year_range <- c(1961, 1999)
threshold_temp <- 5

#_______________________________________________________________________________
#### 3. Delineating by consecutive thermal boundaries (>threshold_temp for >5 days) ####

season_consec_caps.f <- function(dat, latitudes, pops, minyear, maxyear, threshold_temp){
  names(dat) <- NULL
  dat <- as.vector(dat[,1])
  df <- colsplit(dat, ";", c("date", c(paste0("pop", pops))))
  df <- separate(df, col = date, into = c("year", "month", "day"), sep = "\\-")
  
  ## Subsetting years in the dataset
  df <- subset(df, year >= minyear & year <= maxyear)
  
  ## EXTRACTION of growing season bounds  & GSL (growing season length)
 
  years <- unique(df$year)
  
  spring_start <- matrix(NA, ncol = 1 + length(latitudes), nrow = length(years)) # empty results mat for start of spring -- "1 +" because first column is years
  autumn_start <- matrix(NA, ncol = 1 + length(latitudes), nrow = length(years)) # empty results mat for start of spring -- "1 +" because first column is years
  gsl <- matrix(NA, ncol = 1 + length(latitudes), nrow = length(years)) # same as above, but for growing season length (autumn start - spring start)
  
  for (i in 1:length(years)){
    
    annum_temp <- subset(df, df$year == years[i]) # subset df into year

    spring_start[i,1] <- years[i]
    autumn_start[i,1] <- years[i]
    
    gsl[i,1] <- years[i]
    
    cap_length <- 5 # number of consecutive days we want surpassing threshold temp
    
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
    
    for (l in 1:length(latitudes)){
      spring_start[i,1+l] <- find_spring_cap(annum_temp[,3+l])
      autumn_start[i,1+l] <- find_autumn_cap(annum_temp[,3+l])
      gsl[i,1+l] <- find_autumn_cap(annum_temp[,3+l]) - find_spring_cap(annum_temp[,3+l])
      
    }
  }
  spring_start <- as.data.frame(spring_start); colnames(spring_start) <- c("year", c(pops))
  autumn_start <- as.data.frame(autumn_start); colnames(autumn_start) <- c("year", c(pops))
  gsl <- as.data.frame(gsl); colnames(gsl) <- c("year", c(pops))
  data.frame(apply(spring_start, 2, function(x) as.numeric(as.numeric(x))))
  
  return(list(spring_start, autumn_start, gsl))
}

#_______________________________________________________________________________
#### 4. Thermal boundaries for populations ####

lythrum_seasoncaps_predict.list <- season_consec_caps.f(pops_temps, latitudes_pops, pops, min(year_range), max(year_range), threshold_temp)
save(lythrum_seasoncaps_predict.list, file = "data/lythrum_seasoncaps_predict.list.RData") # save to access "spring start" estimates for in situ flowering estimates

umea_predict.list <- season_consec_caps.f(umea_cg_temps, 1, 1, 1999, 2001, threshold_temp)
save(umea_predict.list, file = "data/umea_predict.list.RData") # save to access common garden SOS for comparison with in situ flowering estimates

cap_gsl <- lythrum_seasoncaps_predict.list[[3]][,-1]
cap_sos <- lythrum_seasoncaps_predict.list[[1]][c(37),-1] # row 37 = year 1997, the collection year

cap_gsl_pops <- tibble(pop = as.numeric(colnames(cap_gsl)),
                   capped_season_length =
                     apply(cap_gsl, 2, function(x){mean(as.numeric(x))}),
                   capped_sos = 
                     apply(cap_sos, 2, function(x){mean(as.numeric(x))}))

#_______________________________________________________________________________
#### 5. Merging season length to the flowering data ####

lythrum99 <- read_csv("../lythrum_phenology_99.csv", na = c("", " ", ".")) %>% 
  dplyr::select(pop, plant, morph, flowering_start = flwrstar, flowering_length = flwrleng) %>% 
  na.omit() %>% 
  left_join(y = cap_gsl_pops, by = "pop")

## Save
save(lythrum99, file = "data/lythrum_99.RData")

#_______________________________________________________________________________
#### 6. Spatial data ####

lythrum_coords<- lythrum_coordinates %>%
  left_join(x = ., y = cap_gsl_pops, by = "pop")

save(lythrum_coords, file = "data/lythrum_population_coordinates.RData")

