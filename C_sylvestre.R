#################
## Analyses and plots for Chordeuma sylvestre note
## 
## author: Willson Gaul willson.gaul@ucdconnect.ie
## created: 28 March 2020
## last modified: 7 April 2020
#################
library(wgutil)
library(Hmisc)
library(tidyverse)
library(lubridate)
library(SpadeR)

setwd("~/Documents/UCD/PhD_Project/Chordeuma_sylvestre/")

n_perm <- 100 # number of permutations to use for sp accumulation curves

### load millipede data
# NBDC training and test data
load("~/Documents/Data_Analysis/UCD/NBDC/data/6_Oct_download/NBDC_training_data.Rdata")
load("~/Documents/Data_Analysis/UCD/NBDC/data/6_Oct_download/test_vault/NBDC_test_data.Rdata")

gbif <- read_delim("./GBIF/occurrence.txt", delim = "\t")
gbif <- gbif[gbif$taxonRank == "SPECIES" & 
               gbif$hasGeospatialIssues == "false" & 
               !is.na(gbif$eventDate), ]
gbif <- select(gbif, gbifID, basisOfRecord, occurrenceID, sex, 
               associatedSequences, eventDate, year, month, 
               day, locationID, countryCode, locality, decimalLatitude, 
               decimalLongitude, coordinateUncertaintyInMeters,
               identifiedBy,
               taxonID, scientificName, acceptedScientificName)

hec_names <- read_csv("~/Documents/Data_Analysis/UCD/maps_of_ignorance/millipede_maps_of_ignorance/data/Irish_land_hectads.csv")
# hopefully make hectads line up with grid cells of predictor raster brick
# hec_names$eastings = hec_names$eastings - 4900 
# hec_names$northings = hec_names$northings - 4900

# combine test and training data for all groups
test_data <- test_data$millipede
training_data <- training_data$millipede
mill <- bind_rows(test_data, training_data)
rm(test_data, training_data)

### prepare millipede data ----------------------------------------------------
## make checlist ID variable
mill$checklist_ID <- paste0(mill$RecorderHash, mill$StartDate, mill$EndDate, 
                            mill$SiteName, mill$Latitude, mill$Longitude)
## calculate list lengths for each checklist
mill$list_length <- NA
for(i in 1:length(unique(mill$checklist_ID))) {
  id <- unique(mill$checklist_ID)[i]
  ll <- length(which(mill$checklist_ID == id))
  mill$list_length[mill$checklist_ID == id] <- ll
}

## make Julian day and year variables
mill$year <- year(mill$StartDate)
mill$day_of_year <- yday(mill$StartDate)
mill$temp_resolution <- mill$EndDate - mill$StartDate

# order by date
mill <- mill[order(mill$StartDate), ]

# remove non-species taxa
mill <- mill[-which(is.na(mill$Species)), ]

m1 <- mill[mill$year < 1985, ]
m2 <- mill[mill$year > 1985 & mill$year < 2006, ]

m1_wide <- select(m1, checklist_ID, Genus_species) %>%
  mutate(present = 1) %>%
  unique() %>%
  spread(key = Genus_species, value = present, fill = 0) %>%
  select(-checklist_ID)

m2_wide <- select(m2, checklist_ID, Genus_species) %>%
  mutate(present = 1) %>%
  unique() %>%
  spread(key = Genus_species, value = present, fill = 0) %>%
  select(-checklist_ID)

table(mill$Genus_species)[order(as.numeric(table(mill$Genus_species)))]


# species richness estimates for the 2 time periods using records --------------
m1_sp <- ChaoSpecies(t(m1_wide), datatype = "incidence_raw", k = 10)
m2_sp <- ChaoSpecies(t(m2_wide), datatype = "incidence_raw", k = 10)

## species richness estimates with increasing sample size
n_records <- seq(from = 200, to = 1000, by = 200)
m1_rare <- lapply(n_records, FUN = function(nr, dat) {
  # select subset of records
  dat <- dat[sample(1:nrow(dat), size = nr, replace = FALSE), ] 
  # spread
  dat <- select(dat, checklist_ID, Genus_species) %>%
    mutate(present = 1) %>%
    unique() %>%
    spread(key = Genus_species, value = present, fill = 0) %>%
    select(-checklist_ID)
  # estimate species richness
  dat_sp <- ChaoSpecies(t(dat), datatype = "incidence_raw", k = 10)
  dat_sp <- data.frame(dat_sp$Species_table)
  dat_sp$model <- rownames(dat_sp)
  dat_sp$nrecords <- nr
  dat_sp
}, dat = m1)
m1_rare <- bind_rows(m1_rare)
m1_rare_wide <- pivot_longer(m1_rare, cols = 1:4, names_to = "metric", 
                             values_to = "value") %>%
  filter(metric != "X...s.e." & !grepl(".*Homogeneous Model.*|.*iChao2.*", 
                                       model))
# for 2nd time period
n_records <- seq(from = 200, to = 3200, by = 200)
m2_rare <- lapply(n_records, FUN = function(nr, dat) {
  # select subset of records
  dat <- dat[sample(1:nrow(dat), size = nr, replace = FALSE), ] 
  # spread
  dat <- select(dat, checklist_ID, Genus_species) %>%
    mutate(present = 1) %>%
    unique() %>%
    spread(key = Genus_species, value = present, fill = 0) %>%
    select(-checklist_ID)
  # estimate species richness
  dat_sp <- ChaoSpecies(t(dat), datatype = "incidence_raw", k = 10)
  dat_sp <- data.frame(dat_sp$Species_table)
  dat_sp$model <- rownames(dat_sp)
  dat_sp$nrecords <- nr
  dat_sp
}, dat = m2)
m2_rare <- bind_rows(m2_rare)
m2_rare_wide <- pivot_longer(m2_rare, cols = 1:4, names_to = "metric", 
                             values_to = "value") %>%
  filter(metric != "X...s.e." & !grepl(".*Homogeneous Model.*|.*iChao2.*", 
                                       model))

## make df of species accumulations using records
sp_acc_m1 <- data.frame(record = 1:nrow(m1), obs_n_sp = NA)
for(i in 1:nrow(sp_acc_m1)) {
  sp_acc_m1$obs_n_sp[i] = length(unique(m1$Genus_species[1:i]))
}
# add random accumulations
for(i in 1:n_perm) {
  perm_sp <- m1$Genus_species # get original species records
  perm_sp <- sample(perm_sp, replace = F) # randomly permute species records
  perm_accum <- c() # make empty vector to hold cumulative number of species
  for(j in 1:length(perm_sp)) {
    # calculate cumulative number of species
    perm_accum[j] <- length(unique(perm_sp[1:j]))
  }
  # add permutation accumulation as a column in sp_acc_m1 data frame
  sp_acc_m1[, ncol(sp_acc_m1) + 1] <- perm_accum
}
sp_acc_m1_long <- pivot_longer(sp_acc_m1, 
                               cols = 2:ncol(sp_acc_m1), 
                               names_to = "iteration", 
                               values_to = "cum_sp")

sp_acc_m2 <- data.frame(record = 1:nrow(m2), obs_n_sp = NA)
for(i in 1:nrow(sp_acc_m2)) {
  sp_acc_m2$obs_n_sp[i] = length(unique(m2$Genus_species[1:i]))
}
# add random accumulations
for(i in 1:n_perm) {
  perm_sp <- m2$Genus_species # get original species records
  perm_sp <- sample(perm_sp, replace = F) # randomly permute species records
  perm_accum <- c() # make empty vector to hold cumulative number of species
  for(j in 1:length(perm_sp)) {
    # calculate cumulative number of species
    perm_accum[j] <- length(unique(perm_sp[1:j]))
  }
  # add permutation accumulation as a column in sp_acc_m1 data frame
  sp_acc_m2[, ncol(sp_acc_m2) + 1] <- perm_accum
}
sp_acc_m2_long <- pivot_longer(sp_acc_m2, 
                               cols = 2:ncol(sp_acc_m2), 
                               names_to = "iteration", 
                               values_to = "cum_sp")

## get only records from the area within 20km of where I found C. sylvestre
nearby_m1 <- m1[m1$eastings > 318000 - 20000 &
                  m1$eastings < 318000 + 20000 & 
                  m1$northings > 230000 - 20000 & 
                  m1$northings < 230000 + 20000, ]
nearby_m2 <- m2[m2$eastings > 318000 - 20000 &
                  m2$eastings < 318000 + 20000 & 
                  m2$northings > 230000 - 20000 & 
                  m2$northings < 230000 + 20000, ]

nearby_m1_wide <- select(nearby_m1, checklist_ID, Genus_species) %>%
  mutate(present = 1) %>%
  unique() %>%
  spread(key = Genus_species, value = present, fill = 0) %>%
  select(-checklist_ID)
nearby_m2_wide <- select(nearby_m2, checklist_ID, Genus_species) %>%
  mutate(present = 1) %>%
  unique() %>%
  spread(key = Genus_species, value = present, fill = 0) %>%
  select(-checklist_ID)

# look at species richness estimates for 20km radius area
nearby_m1_sp <- ChaoSpecies(t(nearby_m1_wide), datatype = "incidence_raw", 
                            k = 10)
nearby_m2_sp <- ChaoSpecies(t(nearby_m2_wide), datatype = "incidence_raw", 
                            k = 10)
nearby_m1_sp
nearby_m2_sp

## do random sp accumulation curves for nearby records -------------------------
sp_acc_nearby_m1 <- data.frame(record = 1:nrow(nearby_m1), obs_n_sp = NA)
for(i in 1:nrow(sp_acc_nearby_m1)) {
  sp_acc_nearby_m1$obs_n_sp[i] = length(unique(nearby_m1$Genus_species[1:i]))
}
# add random accumulations
for(i in 1:n_perm) {
  perm_sp <- nearby_m1$Genus_species # get original species records
  perm_sp <- sample(perm_sp, replace = F) # randomly permute species records
  perm_accum <- c() # make empty vector to hold cumulative number of species
  for(j in 1:length(perm_sp)) {
    # calculate cumulative number of species
    perm_accum[j] <- length(unique(perm_sp[1:j]))
  }
  # add permutation accumulation as a column in sp_acc_nearby_m1 data frame
  sp_acc_nearby_m1[, ncol(sp_acc_nearby_m1) + 1] <- perm_accum
}

sp_acc_nearby_m2 <- data.frame(record = 1:nrow(nearby_m2), obs_n_sp = NA)
for(i in 1:nrow(sp_acc_nearby_m2)) {
  sp_acc_nearby_m2$obs_n_sp[i] = length(unique(nearby_m2$Genus_species[1:i]))
}
# add random accumulations
for(i in 1:n_perm) {
  perm_sp <- nearby_m2$Genus_species # get original species records
  perm_sp <- sample(perm_sp, replace = F) # randomly permute species records
  perm_accum <- c() # make empty vector to hold cumulative number of species
  for(j in 1:length(perm_sp)) {
    # calculate cumulative number of species
    perm_accum[j] <- length(unique(perm_sp[1:j]))
  }
  # add permutation accumulation as a column in sp_acc_m1 data frame
  sp_acc_nearby_m2[, ncol(sp_acc_nearby_m2) + 1] <- perm_accum
}

# pivot long
sp_acc_nearby_m1_long <- pivot_longer(sp_acc_nearby_m1, 
                                      cols = 2:ncol(sp_acc_nearby_m1), 
                                      names_to = "iteration", 
                                      values_to = "cum_sp")
sp_acc_nearby_m2_long <- pivot_longer(sp_acc_nearby_m2, 
                                      cols = 2:ncol(sp_acc_nearby_m2), 
                                      names_to = "iteration", 
                                      values_to = "cum_sp")
sp_acc_nearby_m1_long$period <- 1
sp_acc_nearby_m2_long$period <- 2
sp_acc_nearby_all <- bind_rows(sp_acc_nearby_m1_long, sp_acc_nearby_m2_long)

hist(nearby_m1$day_of_year, breaks = 20)
hist(nearby_m2$day_of_year, breaks = 20)

### species accumulation and species richness by adding grid squares ----------
m1_bySite_wide <- filter(m1, Precision == 1000) %>%
  select(checklist_ID, GridReference, Genus_species) %>%
  mutate(present = 1) %>%
  unique() %>%
  spread(key = Genus_species, value = present, fill = 0) %>%
  select(-checklist_ID) %>%
  group_by(GridReference) %>%
  summarise_at(vars(names(.)[2]:names(.)[length(names(.))]), sum)

m2_bySite_wide <- filter(m2, Precision == 1000) %>%
  select(checklist_ID, GridReference, Genus_species) %>%
  mutate(present = 1) %>%
  unique() %>%
  spread(key = Genus_species, value = present, fill = 0) %>%
  select(-checklist_ID) %>%
  group_by(GridReference) %>%
  summarise_at(vars(names(.)[2]:names(.)[length(names(.))]), sum)

# reduce counts to p/a
m1_bySite_wide[, 2:ncol(m1_bySite_wide)] <- pa(
  m1_bySite_wide[, 2:ncol(m1_bySite_wide)])
m2_bySite_wide[, 2:ncol(m2_bySite_wide)] <- pa(
  m2_bySite_wide[, 2:ncol(m2_bySite_wide)])

# estimate species richness when adding sites (1 km grid squares)
m1_bySite_sp <- ChaoSpecies(t(m1_bySite_wide[, 2:ncol(m1_bySite_wide)]), 
                            datatype = "incidence_raw", k = 10)
m2_bySite_sp <- ChaoSpecies(t(m2_bySite_wide[, 2:ncol(m2_bySite_wide)]), 
                            datatype = "incidence_raw", k = 10)

### make df of species accumulations for graphing 
m1_bySite_long <- filter(m1, Precision == 1000) %>%
  select(GridReference, Genus_species) %>% 
  unique()

# get a data frame with a row for each grid cell
sp_acc_bySite_m1 <- data.frame(
  record = 1:length(unique(m1_bySite_long$GridReference)), 
  GridReference = unique(m1_bySite_long$GridReference), 
  obs_n_sp = NA)
# count number of species as each row (grid cell) is added
for(i in 1:nrow(sp_acc_bySite_m1)) {
  sp_acc_bySite_m1$obs_n_sp[i] = length(unique(
    m1_bySite_long$Genus_species[m1_bySite_long$GridReference %in% 
                                   sp_acc_bySite_m1$GridReference[1:i]]))
}

# add random accumulations
for(i in 1:n_perm) {
  perm_site <- sp_acc_bySite_m1$GridReference # get original GridReference order
  perm_site <- sample(perm_site, replace = F) # randomly permute sites
  perm_accum <- c() # make empty vector to hold cumulative number of species
  for(j in 1:length(perm_site)) {
    # calculate cumulative number of species
    perm_accum[j] <- length(unique(
      m1_bySite_long$Genus_species[m1_bySite_long$GridReference %in% 
                                     perm_site[1:j]]))
  }
  # add permutation accumulation as a column in sp_acc_bySite_m1 data frame
  sp_acc_bySite_m1[, ncol(sp_acc_bySite_m1) + 1] <- perm_accum
}
sp_acc_bySite_m1_long <- pivot_longer(sp_acc_bySite_m1, 
                               cols = 3:ncol(sp_acc_bySite_m1), 
                               names_to = "iteration", 
                               values_to = "cum_sp")

# period 2 
m2_bySite_long <- filter(m2, Precision == 1000) %>%
  select(GridReference, Genus_species) %>% 
  unique()

# get a data frame with a row for each grid cell
sp_acc_bySite_m2 <- data.frame(
  record = 1:length(unique(m2_bySite_long$GridReference)), 
  GridReference = unique(m2_bySite_long$GridReference), 
  obs_n_sp = NA)
# count number of species as each row (grid cell) is added
for(i in 1:nrow(sp_acc_bySite_m2)) {
  sp_acc_bySite_m2$obs_n_sp[i] = length(unique(
    m2_bySite_long$Genus_species[m2_bySite_long$GridReference %in% 
                                   sp_acc_bySite_m2$GridReference[1:i]]))
}

# add random accumulations
for(i in 1:n_perm) {
  perm_site <- sp_acc_bySite_m2$GridReference # get original GridReference order
  perm_site <- sample(perm_site, replace = F) # randomly permute sites
  perm_accum <- c() # make empty vector to hold cumulative number of species
  for(j in 1:length(perm_site)) {
    # calculate cumulative number of species
    perm_accum[j] <- length(unique(
      m2_bySite_long$Genus_species[m2_bySite_long$GridReference %in% 
                                     perm_site[1:j]]))
  }
  # add permutation accumulation as a column in sp_acc_bySite_m1 data frame
  sp_acc_bySite_m2[, ncol(sp_acc_bySite_m2) + 1] <- perm_accum
}
sp_acc_bySite_m2_long <- pivot_longer(sp_acc_bySite_m2, 
                                      cols = 3:ncol(sp_acc_bySite_m2), 
                                      names_to = "iteration", 
                                      values_to = "cum_sp")

### end species accumulation by adding grid squares ---------------------------


#### plots ------------------------------------------------------------------
t_size <- 25

## histogram of record dates
year_hist <- ggplot(data = mill,
                    aes(x = year)) + 
  geom_histogram(bins = 45) + 
  xlab("Year") + ylab("Number of records") + 
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(1.8*t_size, "points"))
year_hist

## species accumulation by records in both time periods
sp_acc_m1_long$period <- 1
sp_acc_m2_long$period <- 2
sp_acc_all_long <- bind_rows(sp_acc_m1_long, sp_acc_m2_long)

ggplot(data = sp_acc_all_long, aes(x = record, y = cum_sp)) + 
  geom_line(aes(group = iteration)) + 
  geom_smooth(data = sp_acc_all_long, aes(x = record, y = cum_sp)) +
  geom_line(data = sp_acc_all_long[sp_acc_all_long$iteration == "obs_n_sp", ], 
            aes(x = record, y = cum_sp), color = "orange", size = 1.1) + 
  facet_wrap(~factor(as.character(period), levels = c("1", "2"), 
                     labels = c("1971 to 1984", "1986 to 2005"))) + 
  xlab("Number of records") + ylab("Cumulative number of species") + 
  theme_bw() + 
  theme(legend.position = "none")

# ## These show that the sp. richness estimates are not "creeping up" with sample
# # size, but seem to be settled, at least in time period 2
# ggplot(data = m1_rare_wide, aes(x = nrecords, y = value, color = metric)) + 
#   geom_point() + 
#   facet_wrap(~model) + 
#   ylim(c(0, 50))
# 
# ggplot(data = m2_rare_wide, aes(x = nrecords, y = value, color = metric)) + 
#   geom_point() + 
#   facet_wrap(~model) + 
#   ylim(c(0, 100))

# # plot location of surveys by year in period 2
# ggplot(data = m2, aes(x = eastings, y = northings, 
#                       color = as.numeric(as.character(year)))) + 
#   geom_point()

# Nearby species accumulation curves
ggplot(data = sp_acc_nearby_all, aes(x = record, y = cum_sp)) + 
  geom_line(aes(group = iteration)) + 
  geom_smooth(data = sp_acc_nearby_all, aes(x = record, y = cum_sp)) + 
  geom_line(data = sp_acc_nearby_all[sp_acc_nearby_all$iteration == 
                                       "obs_n_sp", ], 
            aes(x = record, y = cum_sp), color = "orange", size = 1.1) + 
  facet_wrap(~factor(as.character(period), levels = c("1", "2"), 
                     labels = c("1971 to 1984", "1986 to 2005"))) + 
  xlab("Number of records") + ylab("Cumulative number of species") + 
  ggtitle("Recording within 20 km of new C. sylvestre record") + 
  theme_bw() + 
  theme(legend.position = "none")


## species accumulation by adding sites (grid cells) for both time periods
sp_acc_bySite_m1_long$period <- "1"
sp_acc_bySite_m2_long$period <- "2"

sp_acc_bySite_all_long <- bind_rows(sp_acc_bySite_m1_long, sp_acc_bySite_m2_long)

ggplot(data = sp_acc_bySite_all_long[sp_acc_bySite_all_long$iteration != 
                                       "obs_n_sp", ], 
       aes(x = record, y = cum_sp)) + 
  geom_line(aes(group = iteration)) + 
  geom_smooth(aes(x = record, y = cum_sp)) +
  facet_wrap(~factor(period, levels = c("1", "2"), 
                     labels = c("1971 to 1985", "1986 to 2005"))) + 
  xlab("Number of 1 km grid cells sampled") + 
  ylab("Cumulative number of species") +
  theme(legend.position = "none")

### numbers for paper -----------------------------------------------------
nrow(mill) # number of records
length(unique(mill$hectad)) # number of hectads with at least one record
summary(as.numeric(table(mill$hectad)))
hist(as.numeric(table(mill$hectad)), breaks = 50) 

table(mill$year)

# species richness in Ireland in two time periods
m1_sp # species richness 1971 to 1984
m2_sp # species richness 1986 to 2005

## How many doubletons and singletons?  (To make sure Chao can work ok)
table(m1$Genus_species)[order(table(m1$Genus_species))] # time period 1
table(m2$Genus_species)[order(table(m2$Genus_species))] # time period 2


table(gbif$countryCode)


### print plots ---------------------------------------------------------------
ggsave("year_hist.jpg", year_hist, device = "jpeg", width = 25, height = 25, 
       units = "cm")
