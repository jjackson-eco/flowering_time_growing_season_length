#####################################################
##                                                 ##
##        SOLIDAGO Optimal Flowering Time          ##
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
library(lubridate)

#_______________________________________________________________________________
#### 1. Population coordinates ####

solidago_coords <- read_csv("data/Solidago_pops_GPS-coordinates.csv")

## Coordinates as decimals
solidago_coords <- solidago_coords %>% 
  slice(1:24) %>%
  select(pop = 1, 3:4) %>% 
  rename(habitat = ...3) %>%
  separate(`latitude`, paste("lat",c("d","m","s"), sep="_")) %>%
  mutate(across(starts_with("lat_"), as.numeric)) %>%  
  mutate(lat_dec = lat_d + lat_m / 60 + lat_s / 60^2)

#_______________________________________________________________________________
#### 2. Solidago flowering data ####

# Read raw data & summarize mean and sd by population
solidago_dat <- read.csv("data/Solidago_FloweringStartUppsalaGarden2005.csv")
solidago_dat <- solidago_dat[, c("pop", "ind", "FlStart")] %>%
  setNames(c("pop", "ind", "flowering_start"))
solidago_dat$flowering_start <- as.Date(paste0("2005-", solidago_dat$flowering_start), format = "%Y-%d-%b")
solidago_dat$flowering_start <- as.numeric(format(solidago_dat$flowering_start, "%j"))
solidago_dat <- solidago_dat[!is.na(solidago_dat$flowering_start),]

solidago_pops <- sort(unique(solidago_dat$pop))

#_______________________________________________________________________________
#### 3. Temperature data ####

# SMHI data for the 24 Solidago source populations used in Uppsala 2005 common garden
# 1961-2022 daily temps
solidago_pops_temps <- read.csv("data/SMHI_pthbv_t_1961_2022_daily_4326-Solidagopops.csv")

# Uppsala temps
uppsala_cg_temps <- read.csv("data/uppsala_commongarden_1961_2005_daily_4326.csv") # !! temp data should be in order of pop index (1-24) in solidago data, but check !! 
year_range <- c(1961, 2003)
threshold_temp <- 5

#_______________________________________________________________________________
#### 4. Delineating by consecutive thermal boundaries (>threshold_temp for >5 days) ####

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

#_______________________________________________________________________________
#### 5. Thermal boundaries for populations ####

solidago_seasoncaps_predict.list <- season_consec_caps.f(solidago_pops_temps, solidago_pops, min(year_range), max(year_range), threshold_temp)
save(solidago_seasoncaps_predict.list, file = "data/solidago_seasoncaps_predict.list.RData") # save to access "spring start" estimates for in situ flowering estimates

uppsala_predict.list <- season_consec_caps.f(uppsala_cg_temps, 1, 1999, 2005, threshold_temp)
save(uppsala_predict.list, file = "data/uppsala_predict.list.RData") # save to access common garden SOS for comparison with in situ flowering estimates

solidago_cap_gsl <- solidago_seasoncaps_predict.list[[3]][,-1]
solidago_cap_sos <- solidago_seasoncaps_predict.list[[1]][c(41),-1] # row 42 = year 2001, the collection year

solidago_cap_gsl_pops <- tibble(pop = as.numeric(colnames(solidago_cap_gsl)),
                                capped_season_length =
                                  apply(solidago_cap_gsl, 2, function(x){mean(as.numeric(x))}),
                                capped_sos = 
                                  apply(solidago_cap_sos, 2, function(x){mean(as.numeric(x))}))

#_______________________________________________________________________________
#### 5. Merging season length to the flowering data ####

pop_info <- read.csv("data/Solidago_pops_GPS-coordinates.csv")
names(pop_info)[names(pop_info) == "Population"] <- "pop"
pop_info <- pop_info[1:length(solidago_pops), c("pop","X.1")]

solidago05 <- solidago_dat %>%
  left_join(solidago_cap_gsl_pops, by = "pop") %>%
  left_join(pop_info, by = "pop") %>%
  rename(habitat = X.1) %>%
  filter(habitat == "boreal") %>% # subset to only data from boreal populations
  mutate(flowering_start = flowering_start - min(flowering_start, na.rm = T) + 1) # make relative to first record of flowering, to match structure of Lythrum data

## Save
save(solidago05, file = "data/solidago_05.RData")

#_______________________________________________________________________________
#### 6. Spatial data ####

solidago_coordinates <- read_csv("data/Solidago_pops_GPS-coordinates.csv") %>% 
  separate(`latitude`, paste("lat",c("d","m","s"), sep="_")) %>%
  separate(`longitude`, paste("lon",c("d","m","s"), sep="_")) %>%
  mutate(across(lat_d:lon_s, as.numeric)) %>%  
  mutate(lat_dec = lat_d + lat_m / 60 + lat_s / 60^2,
         lon_dec = lon_d + lon_m / 60 + lon_s / 60^2) %>% 
  filter(...3 == "boreal") %>% 
  dplyr::select(Population, lat_dec, lon_dec) %>% 
  left_join(x = ., y = dplyr::select(solidago_cap_gsl_pops, 1:2),
            by = c("Population" = "pop"))

save(solidago_coordinates, file = "data/solidago_population_coordinates.RData")
  

