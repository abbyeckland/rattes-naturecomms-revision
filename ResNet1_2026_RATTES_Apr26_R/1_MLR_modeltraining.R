
# Import packages ---------------------------------------------------------

library(leaps)
library(readr)
library(corrplot)
library(ggplot2)
library(GGally)
library(caret)
library(MASS)
library(factoextra)
library(cluster)
library(reshape2)
library(vip)
library(stringr)
library(data.table)
library(vip)
library(dplyr)

# Version control  ------------------------------------------------------

# update version of code below!
# options are rattes_v1p2_silt and rattes_v1p2_3param_silt, since we train the model using the silt assumption
version <- "rattes_v1p2_silt" # main code version / used to generate results with top 4 predictors

version_settings <- list(
  "rattes_v1p2_silt" = list(
    rattes_version = "rattes_v1p2_silt",
    pred_version   = "rattes_v1p2_silt",
    sy_input_file  = "SYfiles/MLR_Input_silt.csv",
    sa_input_file  = "SYfiles/MLR_SAatDamInitial_silt.csv", 
    predictors     = c("lSR", "yrs_tot", "AveTrap", "wSAatdam", "fdar:indicator_fdar")),
  
  "rattes_v1p2_clay" = list(
    rattes_version = "rattes_v1p2_clay",
    pred_version   = "rattes_v1p2_clay",
    sy_input_file  = "SYfiles/MLR_Input_clay.csv",
    sa_input_file  = "SYfiles/MLR_SAatDamInitial_clay.csv", 
    predictors     = c("lSR", "yrs_tot", "AveTrap", "wSAatdam", "fdar:indicator_fdar")),
  
  "rattes_v1p2_sand" = list(
    rattes_version = "rattes_v1p2_sand",
    pred_version   = "rattes_v1p2_sand",
    sy_input_file  = "SYfiles/MLR_Input_sand.csv",
    sa_input_file  = "SYfiles/MLR_SAatDamInitial_sand.csv", 
    predictors     = c("lSR", "yrs_tot", "AveTrap", "wSAatdam", "fdar:indicator_fdar")),
  
  "rattes_v1p2_3param_silt" = list(
    rattes_version = "rattes_v1p2_3param_silt",
    pred_version   = "rattes_v1p2_3param_silt",
    sy_input_file  = "SYfiles/MLR_Input_silt.csv",
    sa_input_file  = "SYfiles/MLR_SAatDamInitial_silt.csv", 
    predictors     = c("lSR", "AveTrap", "wSAatdam", "fdar:indicator_fdar")))

# Stop early if version is not recognized
if (!version %in% names(version_settings)) {
  stop(paste0(
    "Version '", version, "' not recognized. Valid options are: ",
    paste(names(version_settings), collapse = ", ")
  ))
}

# Pull settings for selected version
settings <- version_settings[[version]]

rattes_version <- settings$rattes_version
pred_version   <- settings$pred_version
sy_input_file  <- settings$sy_input_file


# Print selected settings so user knows what is running
message("Running version: ", rattes_version)
message("Using SY input file: ", sy_input_file)
message("Using predictors: ", paste(settings$predictors, collapse = ", "))


# Data import and merging ------------------------------------------------------

### IMPORT SITE DATA
# set wd to the current folder
setwd(".")

# import study site data for mlr model building (we will use this dataframe later when making predictions)
ss_data <- read.csv(sy_input_file, header=TRUE, stringsAsFactors = FALSE)  
ss_data <- ss_data %>% rename("MAQ_NHDcms" = "MAQ_NHDcfs")  # ensure label is correct for proper conversions (units are indeed cms)

####### make sure 'readin_lakecat_data.R' is updated with current model version and site data file before continuing

# then check if the lakecat CSV has already been processed and file exists in folder directory
run_readin_lakecat_data <- function() {
  file_path <- paste0("lakecat_ss_data_", version, ".csv")
  
  # check if file exists already
  if (!file.exists(file_path)) {
    # file does not exist, so run the script
    message("CSV file is not found. Running 'readin_lakecat_data.R'...")
    source("readin_lakecat_data.R")
  } else {
    message("CSV file already exists. Skipping 'readin_lakecat_data.R'.")
  }
}

# call the function
run_readin_lakecat_data()



### IMPORT LAKECAT DATA
# import lakecat_ss_data which binded lakecat data to ss_data
file_path <- paste0("lakecat_ss_data_", version, ".csv")
lakecat_ss_data <- read.csv(file_path, header=TRUE, stringsAsFactors = FALSE)

# import ss_wbcomid, which binded study site sid and wbcomid
file_path <- paste0("ss_wbcomid_data_", version, ".csv")
ss_wbcomid <- read.csv(file_path, header=TRUE, stringsAsFactors = FALSE)

lakecat_ss_data <- lakecat_ss_data %>% rename("MAQ_NHDcms" = "MAQ_NHDcfs")  # ensure label is correct for proper conversions

### IMPORT FCPG DATA
# import fcpg data for model building
fcpg_data <- read.csv("fcpg_data/FCPG_snap_output_12_5_24.csv", header=TRUE, stringsAsFactors = FALSE)
fcpg_data <- fcpg_data[fcpg_data$IsRiverMth==0,]  # remove river mouths
fcpg_cols<- c("slope", "prcp", "elev", "tmin", "tmax", "NID")  # keep only relevant data columns (ignore land cover fcpg data, too much is missing)
fcpg_subset <- fcpg_data[names(fcpg_data) %in% fcpg_cols]
print(nrow(fcpg_subset))

# merge fcpg data with study site and lakecat data
lakecat_ss_fcpg_data <- merge(lakecat_ss_data, fcpg_subset, by="NID", all=TRUE)
lakecat_ss_fcpg_data <- lakecat_ss_fcpg_data[!is.na(lakecat_ss_fcpg_data$SID),]   # filter out rows with nan sid

# filter out prediction sites and lock dams
site_data <- lakecat_ss_fcpg_data[
  lakecat_ss_fcpg_data$issite == 1 & 
    (lakecat_ss_fcpg_data$islock == 0 | (lakecat_ss_fcpg_data$islock == 1 & lakecat_ss_fcpg_data$PermStorag == 1)), 
]

# save queries to folder "info"
site_data <- site_data[!is.na(site_data$SID),]
print(nrow(site_data)) # number of sites
# create info folder if doesn't already exist
if (!dir.exists("info")) {
  dir.create("info")
}
# with locks removed, now see how many sites have lakecat data
site_data_nan_check <- site_data %>% filter(!is.na(WsAreaSqKm))  # check this with lakecat parameter 'WsAreaSqKm'
sink("info/number_of_study_sites_with_lakecat.txt", split = TRUE)
print(nrow(site_data_nan_check)) # number of sites with lakecat data
sink()


# Data pre-processing -----------------------------------------------------


### DEAL WITH MISSING OR BAD PARAMETER DATA
# find the starting index based on column name
start_index <- which(colnames(site_data) == "PctAg2006Slp20Ws") # beginning index of lakecat column data (first parameter in list should be "PctAg2006Slp20Ws")
# create list of columns from that point onward
param_list_lakecat <- colnames(site_data)[start_index:ncol(site_data)]
print(param_list_lakecat)

# add relevant additional parameter columns from 'ss_data' df for pre-processing
additional_cols <- c('damlength_m', 'SA_m2', 'damH', 'NrX_Final', 'NrY_Final', 'DA', 'MAQ_NHDcms', 'D50', 'AveTrap', 'wSAatdam', 'pathx', 'yrc')
param_list <- c(additional_cols, param_list_lakecat)

# print the number of nan values in each column (ideal if no nans in the ss_data columns)
site_data_param_list <- site_data[, param_list]
sapply(site_data_param_list, function(x) sum(is.na(x)))

# drop sites from dataframe that contain all nans, but only for columns within param_list
site_data_nonans <- site_data[!apply(site_data[param_list], 1, function(row) all(is.na(row))), ]
print(nrow(site_data_nonans))

### OUTLIER REMOVAL
# however, do not remove outliers from following columns (have alrady been qa/qc'ed)
no_outlier_removal_from <- c("AveTrap", "wSAatdam", "pathx", "yrc", "DA", "NrX_Final", "NrY_Final")
param_list_outlier_check <- param_list[!param_list %in% no_outlier_removal_from]

# loop over each variable in param_list
sink("info/parameter_outlier_removal_study_sites.txt", split = TRUE)
for (variable in param_list_outlier_check) {
  # exclude NA and zero values
  valid_data <- site_data_nonans[[variable]][!is.na(site_data_nonans[[variable]]) & site_data_nonans[[variable]] != 0]

  # standardize the variable (only valid data)
  z_scores <- scale(valid_data)

  # identify outliers based on z-scores
  outliers <- which(abs(z_scores) > 5)  # make data nan if value is > 5 z-scores from the mean

  # print the number of outliers for the current variable
  cat("Number of outliers in", variable, ":", length(outliers), "\n")

  # set the identified outlier values to NA in the original dataset
  # map the outlier indices back to the original dataset
  all_indices <- which(!is.na(site_data_nonans[[variable]]) & site_data_nonans[[variable]] != 0)
  site_data_nonans[all_indices[outliers], variable] <- NA
}
sink()

### CHECK FOR NEGATIVE DATA
# check if any columns have negative values (avoids any later transformation issues)
negative_columns <- names(site_data_nonans)[sapply(site_data_nonans, function(x) any(x < 0 & !is.na(x)))]
# print the column names with negative values
print(negative_columns)
negative_columns <- negative_columns[negative_columns %in% param_list]
print(negative_columns)     # code below will handle expected negatives in NrX_Final and temperature data (Tmin8110Ws, Tmean8110Ws, and Tmax8110Ws)

# rescale "Tmin8110Ws" column to avoid negatives and save to new column (e.g., "Tmin8110Ws_norm")
min_max_normalize <- function(x) {    # function for min-max normalization
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  # apply normalization only to non-NA values
  normalized <- ifelse(!is.na(x), (x - min_val) / (max_val - min_val), NA)
  normalized <- (x - min_val) / (max_val - min_val)
  normalized <- normalized * (10 - 1.01) + 1.01 # Scale to desired range [1.01-10]
  return(normalized)
}

# apply min-max normalization to temp columns and save as new column
site_data_nonans$Tmin8110Ws_norm <- min_max_normalize(site_data_nonans$Tmin8110Ws)  # lakecat Tmin
site_data_nonans$Tmean8110Ws_norm <- min_max_normalize(site_data_nonans$Tmean8110Ws)  # lakecat Tmean
site_data_nonans$Tmax8110Ws_norm <- min_max_normalize(site_data_nonans$Tmax8110Ws)  # lakecat Tmax
site_data_nonans$tmin_norm <- min_max_normalize(site_data_nonans$tmin)  # fcpg tmin
site_data_nonans$tmax_norm <- min_max_normalize(site_data_nonans$tmax)  # fcpg tmax

# reverse the sign of NrX_Final (longitude)
site_data_nonans$NrX_Final <- abs(site_data_nonans$NrX_Final)


### CHECK FOR MISSING DATA
# check if any columns have a placeholder for missing data. Missing data is represented by -9998 or -9999
count_missing <- colSums(site_data_nonans == -9998 | site_data_nonans == -9999, na.rm=TRUE)
cols_with_missing <- names(count_missing[count_missing > 0])
print(count_missing)
print(cols_with_missing)

# filter param_list to include only columns that are in cols_with_missing
cols_to_check <- param_list[param_list %in% cols_with_missing]  # cols_to_check may be empty, keep code here anyways in case missing data in future
# remove missing data values with NA
site_data_nonans[cols_to_check] <- lapply(site_data_nonans[cols_to_check], function(col) {
  col[col == -9998 | col == -9999] <- NA
  return(col)
})

# remove any parameters that have constant values (standard dev is 0); this is necessary prior to data transformations
# select only numeric columns
numeric_columns <- site_data_nonans[, sapply(site_data_nonans, is.numeric)]
# calculate the standard deviation for each column
std_devs <- sapply(numeric_columns, sd, na.rm = TRUE)
# identify columns with a standard deviation of zero
constant_columns <- names(std_devs[std_devs == 0])
constant_columns <- constant_columns[constant_columns %in% param_list]
print(constant_columns)
# drop the identified columns from the dataframe
site_data_nonans <- site_data_nonans[, !(names(site_data_nonans) %in% constant_columns)]
# update parameter list after removing constant column(s) parameters
param_list <- param_list[!(param_list %in% constant_columns)]
param_list_lakecat <- param_list_lakecat[!(param_list_lakecat %in% constant_columns)]


### CONFIRM NO MISSING DATA IN REQUIRED COLUMNS
# Calculate the percentage of missing data in each column
missing_percent <- colMeans(is.na(site_data_nonans)) * 100
cols_with_missing <- names(missing_percent[missing_percent > 0])
print(cols_with_missing)
# ensure no mandatory columns are missing survey data
ss_cols <- c("cap1", "cap2",  "yr1", "yr2")
cols_with_missing <- cols_with_missing[cols_with_missing %in% ss_cols]
print(cols_with_missing) # should be 0
if (length(cols_with_missing) >= 1 ) {
  stop("Warning: required columns have NAs.")
}

# remove rows with missing data only within cols_with_missing
site_data_nonans <- site_data_nonans[complete.cases(site_data_nonans[cols_with_missing]), ]
sink("info/number_of_study_sites_nolocks_permstorage.txt", split = TRUE)
print(nrow(site_data_nonans))  # final number of sites to build model (904-includes perm storage, no lock dams)
sink()



# Data transformations -----------------------------------------------

### DEVELOP RESPONSE VARIABLE (lSR) (logged sedimentation rate)
# create sedimentation rate parameter
site_data_nonans$yrs_tot <- site_data_nonans$yr2 - site_data_nonans$yr1
site_data_nonans$SV <- site_data_nonans$cap1 - site_data_nonans$cap2
site_data_nonans$SR <- site_data_nonans$SV / site_data_nonans$yrs_tot
site_data_nonans$lSR <- log(site_data_nonans$SR)  # log here is natural log
site_data_nonans$AveTrap <- site_data_nonans$AveTrap*100  # convert avetrap to percentage, to help with data transformations

### PREPARE INDEPENDENT VARIABLES
# create/combine parameters from existing parameters

# combine all forest types into one forest param
site_data_nonans$pct_frst_all_2006 <- rowSums(site_data_nonans[,c("PctDecid2006Ws","PctConif2006Ws","PctMxFst2006Ws")], na.rm=TRUE)
site_data_nonans$pct_frst_all_2006[is.nan(site_data_nonans$PctDecid2006Ws) & is.nan(site_data_nonans$PctConif2006Ws) & is.nan(site_data_nonans$PctMxFst2006Ws)] <- NaN  # if all 3 columns are nan, making resulting column nan
site_data_nonans$pct_frst_all_2011 <- rowSums(site_data_nonans[,c("PctDecid2011Ws","PctConif2011Ws","PctMxFst2011Ws")], na.rm=TRUE)
site_data_nonans$pct_frst_all_2011[is.nan(site_data_nonans$PctDecid2011Ws) & is.nan(site_data_nonans$PctConif2011Ws) & is.nan(site_data_nonans$PctMxFst2011Ws)] <- NaN

# combine all wetland types into one wetland param
site_data_nonans$pct_wetl_all_2006 <- rowSums(site_data_nonans[,c("PctWdWet2006Ws","PctHbWet2006Ws")], na.rm=TRUE)
site_data_nonans$pct_wetl_all_2006[is.nan(site_data_nonans$PctWdWet2006Ws) & is.nan(site_data_nonans$PctHbWet2006Ws)] <- NaN
site_data_nonans$pct_wetl_all_2011 <- rowSums(site_data_nonans[,c("PctWdWet2011Ws","PctHbWet2011Ws")], na.rm=TRUE)
site_data_nonans$pct_wetl_all_2011[is.nan(site_data_nonans$PctWdWet2011Ws) & is.nan(site_data_nonans$PctHbWet2011Ws)] <- NaN

# combine percent glacial till, glacial lake sediments, and percent ice into a single param
site_data_nonans$pct_icy_sed <- rowSums(site_data_nonans[,c("PctGlacTilClayWs","PctGlacTilLoamWs","PctGlacTilCrsWs","PctGlacLakeCrsWs","PctGlacLakeFineWs")], na.rm=TRUE)
site_data_nonans$pct_icy_sed[is.nan(site_data_nonans$PctGlacTilClayWs) & is.nan(site_data_nonans$PctGlacTilLoamWs) & is.nan(site_data_nonans$PctGlacTilCrsWs) & is.nan(site_data_nonans$PctGlacLakeCrsWs) & is.nan(site_data_nonans$PctGlacLakeFineWs)] <- NaN

# combine all eolian lithology types into one eolian param
site_data_nonans$pct_eolian <- rowSums(site_data_nonans[,c("PctEolCrsWs","PctEolFineWs")], na.rm=TRUE)
site_data_nonans$pct_eolian[is.nan(site_data_nonans$PctEolCrsWs) & is.nan(site_data_nonans$PctEolFineWs)] <- NaN

# combine all agriculture on slope params into one ag on slope param
site_data_nonans$pct_ag_on_slp <- rowSums(site_data_nonans[,c("PctAg2006Slp20Ws","PctAg2006Slp10Ws")], na.rm=TRUE)
site_data_nonans$pct_ag_on_slp[is.nan(site_data_nonans$PctAg2006Slp20Ws) & is.nan(site_data_nonans$PctAg2006Slp10Ws)] <- NaN

# combine annual forest fire area and forest loss into a single decadal parameter
site_data_nonans$pct_fire_0010Ws <- rowSums(site_data_nonans[,c("PctFire2000Ws","PctFire2001Ws","PctFire2002Ws","PctFire2003Ws","PctFire2004Ws","PctFire2005Ws","PctFire2006Ws","PctFire2007Ws","PctFire2008Ws","PctFire2009Ws","PctFire2010Ws")], na.rm=TRUE)
site_data_nonans$pct_fire_0010Ws[is.nan(site_data_nonans$PctFire2000Ws) & is.nan(site_data_nonans$PctFire2001Ws) & is.nan(site_data_nonans$PctFire2002Ws) & is.nan(site_data_nonans$PctFire2003Ws) & is.nan(site_data_nonans$PctFire2004Ws) & is.nan(site_data_nonans$PctFire2005Ws) & is.nan(site_data_nonans$PctFire2006Ws) & is.nan(site_data_nonans$PctFire2007Ws) & is.nan(site_data_nonans$PctFire2008Ws) & is.nan(site_data_nonans$PctFire2009Ws) & is.nan(site_data_nonans$PctFire2010Ws)] <- NaN
site_data_nonans$pct_frstloss_0113Ws <- rowSums(site_data_nonans[,c("PctFrstLoss2001Ws","PctFrstLoss2002Ws","PctFrstLoss2003Ws","PctFrstLoss2004Ws","PctFrstLoss2005Ws","PctFrstLoss2006Ws","PctFrstLoss2007Ws","PctFrstLoss2008Ws","PctFrstLoss2009Ws","PctFrstLoss2010Ws","PctFrstLoss2011Ws","PctFrstLoss2012Ws","PctFrstLoss2013Ws")], na.rm=TRUE)
site_data_nonans$pct_frstloss_0113Ws[is.nan(site_data_nonans$PctFrstLoss2001Ws) & is.nan(site_data_nonans$PctFrstLoss2002Ws) & is.nan(site_data_nonans$PctFrstLoss2003Ws) & is.nan(site_data_nonans$PctFrstLoss2004Ws) & is.nan(site_data_nonans$PctFrstLoss2005Ws) & is.nan(site_data_nonans$PctFrstLoss2006Ws) & is.nan(site_data_nonans$PctFrstLoss2007Ws) & is.nan(site_data_nonans$PctFrstLoss2008Ws) & is.nan(site_data_nonans$PctFrstLoss2009Ws) & is.nan(site_data_nonans$PctFrstLoss2010Ws) & is.nan(site_data_nonans$PctFrstLoss2011Ws)& is.nan(site_data_nonans$PctFrstLoss2012Ws) & is.nan(site_data_nonans$PctFrstLoss2013Ws)] <- NaN

# calculate temperature difference for Tmax8110Ws and Tmin8110Ws (lakecat param), only if both are not NaN and not zero
site_data_nonans$tdiff_8110Ws <- ifelse(
  !is.na(site_data_nonans$Tmax8110Ws) & !is.na(site_data_nonans$Tmin8110Ws) & 
    site_data_nonans$Tmax8110Ws != 0 & site_data_nonans$Tmin8110Ws != 0,
  site_data_nonans$Tmax8110Ws - site_data_nonans$Tmin8110Ws,
  NA
)

# calculate temp difference for tmax and tmin (fcpg param), only if both are not NaN and not zero
site_data_nonans$tdiff_fcpg <- ifelse(
  !is.na(site_data_nonans$tmax) & !is.na(site_data_nonans$tmin) & 
    site_data_nonans$tmax != 0 & site_data_nonans$tmin != 0,
  site_data_nonans$tmax - site_data_nonans$tmin,
  NA
)

# convert mean annual streamflow in m3/s to m3/yr, only if MAQ_NHDcms is not NaN and greater than zero
mean_ann_sf <- ifelse(
  !is.na(site_data_nonans$MAQ_NHDcms) & site_data_nonans$MAQ_NHDcms > 0,
  site_data_nonans$MAQ_NHDcms * 31536000,
  NA
)

# calculate average of capacity, only if cap1 and cap2 are not NaN and greater than 0
mean_cap <- ifelse(
  !is.na(site_data_nonans$cap1) & !is.na(site_data_nonans$cap2) & 
    site_data_nonans$cap1 > 0 & site_data_nonans$cap2 > 0,
  (site_data_nonans$cap1 + site_data_nonans$cap2) / 2,
  NA
)

# compute water residence time (wrt) in years, only if mean_ann_sf and mean_cap are not NaN and greater than 0
site_data_nonans$wrt <- ifelse(
  !is.na(mean_ann_sf) & mean_ann_sf > 0 & !is.na(mean_cap) & mean_cap > 0,
  mean_cap / mean_ann_sf,
  NA
)

# compute flow to drainage area ratio, or specific discharge, only if DA and MAQ_NHDcms are not NaN and greater than 0
site_data_nonans$fdar <- ifelse(
  !is.na(site_data_nonans$DA) & !is.na(site_data_nonans$MAQ_NHDcms) &
    site_data_nonans$DA > 0 & site_data_nonans$MAQ_NHDcms > 0,
  (site_data_nonans$MAQ_NHDcms/site_data_nonans$DA)*100,
  NA
)

# calculate wSAratio, only if wSAatdam and DA are not NaN and greater than 0
site_data_nonans$wSAratio <- ifelse(
  !is.na(site_data_nonans$wSAatdam) & !is.na(site_data_nonans$DA) & 
    site_data_nonans$wSAatdam > 0 & site_data_nonans$DA > 0,
  (site_data_nonans$wSAatdam / site_data_nonans$DA) * 100,
  NA
)

# calculate Relief as the difference between ElevWs (entire watershed) and ElevCat (local catchment elevation), only if both are not NaN and greater than 0
site_data_nonans$relief <- ifelse(
  !is.na(site_data_nonans$ElevWs) & !is.na(site_data_nonans$ElevCat) & 
    site_data_nonans$ElevWs > 0 & site_data_nonans$ElevCat > 0 & site_data_nonans$ElevWs > site_data_nonans$ElevCat,
  site_data_nonans$ElevWs - site_data_nonans$ElevCat,
  NA
)

# check how many nans in each param
sapply(site_data_nonans, function(x) sum(is.na(x)))

#now replace any remaining 0 values with small number (transformations don't like if value of 0)
site_data_nonans[site_data_nonans == 0] <- 1e-6


### now remove lSR outliers
# make copy of site_data_nonans df
lSR_outliers_removed <- site_data_nonans

# calculate the z-scores for the 'lSR' column
lSR_outliers_removed <- lSR_outliers_removed %>%
  mutate(lSR_z_score = (lSR - mean(lSR, na.rm = TRUE)) / sd(lSR, na.rm = TRUE))

# filter the dataframe to remove rows where the z-score is greater than absolute value of 3
filtered_site_data_nonans <- lSR_outliers_removed %>%
  filter(abs(lSR_z_score) <= 3)

# remove the temporary z-score column if desired
filtered_site_data_nonans <- filtered_site_data_nonans %>%
  dplyr::select(-lSR_z_score)

# view the filtered dataframe
print(nrow(filtered_site_data_nonans))  # if number is 904, no outliers identified


### PLOT DISTRIBUTION OF lSR VALUES TO CHECK FOR NORMALITY
# q-q plot checks normality of distribution within quantiles
if (!dir.exists("figures/model_building/lsr_normality")) {
  dir.create("figures/model_building/lsr_normality", recursive = TRUE)
}
png("figures/model_building/lsr_normality/QQ_plot_lSR.png", width = 1600, height = 1600, res = 300)
qqnorm(filtered_site_data_nonans$lSR)
qqline(filtered_site_data_nonans$lSR, col = "red")
dev.off()
# print to screen
qqnorm(filtered_site_data_nonans$lSR, main = "QQ Plot of lSR")
qqline(filtered_site_data_nonans$lSR, col = "red")


# histogram also shows normality of distribution via bins
png("figures/model_building/lsr_normality/histogram_plot_lSR.png", width = 1200, height = 1200, res = 150)
p <- ggplot(filtered_site_data_nonans, aes(x = lSR)) +
  geom_histogram(#binwidth = .25,   # adjust binwidth if needed
    fill = "lightblue",
    color = "black") +
  labs(title = "Histogram of logged sedimentation rates (lSR) for study sites",
       x = "lsr",
       y = "Frequency") +
  theme_minimal()  # Use a minimal theme for a clean look
# print to screen
print(p)
dev.off()
print(p)



### PREPARE DATAFRAME FOR TRANSFORMATIONS
# make new dataframe with only desired columns and response variable, lSR
y_cols <- c("lSR")
new_cols <- c("pct_frst_all_2006", "pct_wetl_all_2006",  "pct_icy_sed", "pct_eolian", "pct_ag_on_slp", "pct_fire_0010Ws", "pct_frstloss_0113Ws", "tmin_norm", "tmax_norm", "tdiff_fcpg", "wrt", "yrs_tot", "wSAratio", "relief", "fdar")  # all the derived columns
all_cols <- c(y_cols, param_list, new_cols)
# prepare dataframe
df_ready <- filtered_site_data_nonans[, all_cols]
# define some columns to drop for data we don't need
cols_to_drop_major <- c("COMID_y", "COMID", # don't need these now
                  "DA",  # model is better with wSAatdam, and this introduces colinearities, so removing
                  "Tmin8110Ws", "Tmean8110Ws", "Tmax8110Ws",# using "Tmin8110Ws_norm" et al instead since no negative numbers
                  "tmin", "tmax", # using normalized params instead
                  "PctAg2006Slp20Ws", "PctAg2006Slp10Ws",  # combined these into one column
                  "PctOw2011Ws", "PctIce2011Ws", "PctUrbOp2011Ws", "PctUrbLo2011Ws", "PctUrbMd2011Ws", "PctUrbHi2011Ws", "PctBl2011Ws", "PctDecid2011Ws", "PctConif2011Ws", "PctMxFst2011Ws", "PctShrb2011Ws", "PctGrs2011Ws", "PctHay2011Ws", "PctCrop2011Ws", "PctWdWet2011Ws", "PctHbWet2011Ws", "PctImp2011Ws", # already have these params for 2011. 2006 is better since more likely survey period will cover it
                  "PctDecid2006Ws","PctConif2006Ws","PctMxFst2006Ws", # combined all forest into a single parameter above called "pct_frst_all_2006"
                  "PctDecid2011Ws","PctConif2011Ws","PctMxFst2011Ws", # combined all forest into a single parameter above called "pct_frst_all_2011"
                  "PctWdWet2006Ws","PctHbWet2006Ws", # combined all wetland into a single parameter above called "pct_wetl_all_2006"
                  "PctWdWet2011Ws","PctHbWet2011Ws", # combined all wetland into a single parameter above called "pct_wetl_all_2011"
                  "PctGlacTilClayWs","PctGlacTilLoamWs","PctGlacTilCrsWs","PctGlacLakeCrsWs","PctGlacLakeFineWs", # combined into one icy_sed param
                  "PctEolCrsWs","PctEolFineWs", # combined into one eolian param
                  "PctFire2000Ws", "PctFire2001Ws", "PctFire2002Ws", "PctFire2003Ws", "PctFire2004Ws", "PctFire2005Ws", "PctFire2006Ws", "PctFire2007Ws", "PctFire2008Ws", "PctFire2009Ws", "PctFire2010Ws", # combined into a single decadal param
                  "PctFrstLoss2001Ws", "PctFrstLoss2002Ws", "PctFrstLoss2003Ws", "PctFrstLoss2004Ws", "PctFrstLoss2005Ws", "PctFrstLoss2006Ws", "PctFrstLoss2007Ws", "PctFrstLoss2008Ws", "PctFrstLoss2009Ws", "PctFrstLoss2010Ws", "PctFrstLoss2011Ws", "PctFrstLoss2012Ws", "PctFrstLoss2013Ws",  # combined into a single decadal param
                  "NABD_DensWs", "NABD_NIDStorWs", "NABD_NrmStorWs", # this is basically a repeat of NID that has been slightly filtered
                  "DamDensWs", "DamNIDStorWs", "DamNrmStorWs", # don't want anything else about upstream damming since we are including this separately and accounting for changes with time
                  "CoalMineDensWs", # there is another param called "MineDensWs" that has a better corr with lSR
                  "PctOw2006Ws", "PctOw2011Ws", # there is another param called "PctWaterWs" that has a better corr with lSR
                  "REGION", "LRR_SYMBOL",  # from earlier versions
                  "Precip8110Ws", "ElevWs", "ElevCat") # already have complete datasets of elevation and precip from fcpg (these are the lakecat ones, with some missing data)
df_ready <- df_ready[, !names(df_ready) %in% cols_to_drop_major]

# (optional) save df_ready prior to transformations to file
file_path <- paste0("site_data_pre-transformations_", version, ".csv")
write.csv(df_ready, file=file_path, row.names=FALSE)  # this dataframe contains lSR, and all the IVs prior to data transformations, but does not include any other site data


### REMOVE BAD / MISSING DATA
# check again for any bad data that may have been introduced in derived columns above
# remove any parameters that have constant values (standard dev. is 0)
# select only numeric columns
numeric_columns <- df_ready[, sapply(df_ready, is.numeric)]
# calculate the standard deviation for each column
std_devs <- sapply(numeric_columns, sd, na.rm = TRUE)
# identify columns with a standard deviation of zero
constant_columns <- names(std_devs[std_devs == 0])
constant_columns <- constant_columns[constant_columns %in% param_list]
print(constant_columns) # aim for empty list

# identify columns where all data is NAN 
nan_columns <- names(which(sapply(numeric_columns, function(x) all(is.na(x)))))
print(nan_columns) # aim for empty list

# identify columns with very few unique values (e.g., less than 30 unique values excluding NaN)
few_unique_columns <- names(which(sapply(numeric_columns, function(x) {
  unique_values <- unique(na.omit(x))  # Exclude NaN values
  length(unique_values) <= 30  # Used 30 so when split data into training and test sets, there is enough data within the test set to still perform well
})))
print(few_unique_columns)

# combine all problematic columns
problematic_columns <- unique(c(constant_columns, nan_columns, few_unique_columns))
print(problematic_columns)
# drop the identified columns from the dataframe
df_ready <- df_ready[, !(names(df_ready) %in% problematic_columns)]
# update parameter lists after removing problematic column(s) parameters
param_list <- param_list[!(param_list %in% problematic_columns)]
param_list_lakecat <- param_list_lakecat[!(param_list_lakecat %in% problematic_columns)]


### WRITE FUNCTIONS FOR DATA TRANSFORMATIONS
# generate candidate transformations for lambda values (LOG TRANSFORMATION CONDUCTED, BUT ARCHITECTURE RETAINED IN CASE WANT TO DO BOX-COX IN FUTURE)
lambda1_values <- seq(0, 0, by = 0)      # lambda values between 0 and 0 (forcing a log transformation instead of box-cox for this iteration)
lambda2_values <- seq(0, 0, by = 0)      # potential for creating 'shifting' parameters to handle 0 data, not including for this iteration

# create a grid of all combinations of lambda_1 and lambda_2
transforms_grid <- expand.grid(lambda1 = lambda1_values, lambda2 = lambda2_values)

# generate transformation functions for each combination of lambda_1 and lambda_2
transforms <- lapply(1:nrow(transforms_grid), function(i) {
  lambda1 <- transforms_grid$lambda1[i]
  lambda2 <- transforms_grid$lambda2[i]

  transform <- function(x) {
  if (lambda1 == 0) {
    return(log(x + lambda2))
  } else if (lambda1 == 1) {
    return(x)
  } else {
    return(((x + lambda2)^lambda1 - 1) / lambda1)
    }
  }

  return(list(name = if (lambda1 == 1) "Linear" else "Box-Cox",
              lambda1 = lambda1,
              lambda2 = lambda2,
              transform = transform))
})

# function to compute AIC for a model
compute_AIC <- function(model) {
  if (!is.null(model)) {
    return(AIC(model))
  } else {
    return(NA)
  }
}

# create list of lock SIDs for later that don't contain permanent storage
locks <- ss_data[ss_data$islock==1 & ss_data$PermStorag == 0, ]
lock_sids <- unique(locks$SID)


### IDENTIFY BEST DATA TRANSFORMATION (for this iteration, forcing to be log)
# initialize dataframe to store best identified transformation
output_df_lsr <- data.frame(variable = character(), best_transformation = character(),
                            lambda1 = numeric(), lambda2 = numeric(), AIC = numeric())

# iterate over each predictor variable
for (col_name in names(df_ready)[-c(1)]) {
  cat("Variable:", col_name, "\n")

  # fit models with different transformations for the current predictor variable
  models <- lapply(transforms, function(transform) {
    transform_function <- transform$transform
    # filter out nan values
    non_na_data <- na.omit(df_ready[[col_name]])
    # filter out zero values                                              
    non_na_data <- non_na_data[non_na_data > 1e-6]
    # apply transformations to non-na data
    transformed_data <- transform_function(non_na_data)
    # create a subset of lSR corresponding to non-na data       
    valid_lSR <- df_ready$lSR[!is.na(df_ready[[col_name]]) & (df_ready[[col_name]] > 1e-6)]
    # fit the linear regression model
    lm(valid_lSR ~ transformed_data)
  })

  # compute AIC for each model
  aics <- sapply(models, compute_AIC)
  # find the index of the transformation with the minimum AIC
  best_transform_index <- which.min(aics)
  # extract the corresponding transformation name, lambda value, and AIC
  best_transform <- transforms[[best_transform_index]]$name
  best_lambda1 <- transforms[[best_transform_index]]$lambda1
  best_lambda2 <- transforms[[best_transform_index]]$lambda2
  best_AIC <- aics[best_transform_index]

  # store the results in the output dataframe
  output_df_lsr <- bind_rows(output_df_lsr, data.frame(variable = col_name,
                                                       best_transformation = best_transform,
                                                       lambda1 = best_lambda1,
                                                       lambda2 = best_lambda2,
                                                       AIC = best_AIC))
}

# view the output dataframe with identified best transformations
print(output_df_lsr)

# save a copy of this for later use and reference
output_df_lsr_usa <- output_df_lsr
# save transformation list to file
if (!dir.exists("transformations")) {
  dir.create("transformations")
}
file_path <- paste0("transformations/site_data_transformation_type_", version, ".csv")
write.csv(output_df_lsr_usa, file=file_path, row.names=FALSE)


### APPLY TRANSFORMATIONS TO DATA
# initialize dataframe to store transformed data
transformed_df_lsr <- df_ready
# iterate over each row in output_df_lsr
for (i in 1:nrow(output_df_lsr)) {
  # get the variable name and its corresponding lambda value
  variable_name <- output_df_lsr$variable[i]
  lambda1_value <- output_df_lsr$lambda1[i]
  lambda2_value <- output_df_lsr$lambda2[i]
  # filter out nan values for the current variable
  non_na_values <- na.omit(df_ready[[variable_name]])
  # filter out zero values -                                                
  non_na_values <- non_na_values[non_na_values > 1e-6]

  # apply the corresponding transformation using the lambda value to non-nan values
  if (lambda1_value == 0) {
    transformed_values <- log(non_na_values + lambda2_value)
  } else if (lambda1_value == 1) {
    transformed_values <- non_na_values  # linear transformation
  } else {
    transformed_values <- ((non_na_values + lambda2_value)^lambda1_value - 1)/lambda1_value
  }

  # create a vector with the same length as the original variable, filled with NA
  transformed_vector <- rep(NA, length(df_ready[[variable_name]]))
  # replace the NA positions with the transformed values                   
  valid_indices <- which(!is.na(df_ready[[variable_name]]) & df_ready[[variable_name]] > 1e-6)
  transformed_vector[valid_indices] <- transformed_values

  # update the dataframe with the transformed values
  transformed_df_lsr[[variable_name]] <- transformed_vector
}

# view the transformed dataframe
print(transformed_df_lsr)

# save copy of transformed dataframe
transformed_df_lsr_ <- transformed_df_lsr
# save transformed data to folder
file_path <- paste0("transformations/site_data_transformed_", version, ".csv")
write.csv(transformed_df_lsr_, file=file_path, row.names=FALSE)




# MLR model building, lSR  -------------------------------------------------------

### CREATE INTERACTION TERMS FOR PARAMETERS WITH MISSING DATA AND/OR ONLY APPLY TEMPORALLY
# first, create interaction terms for temporally-applicable columns
# add yr2 column back to dataframe
transformed_df_lsr_$yr2 <- filtered_site_data_nonans$yr2
# rename column with partial year values to full year values (ignore temp and precip params, assuming similar trends over time)
columns_to_rename <- c("Pestic97Ws" = "Pestic1997Ws", "pct_fire_0010Ws" = "pct_fire_decade_2010Ws", "pct_frstloss_0113Ws" = "pct_frstloss_decade_2013Ws")
# rename the columns based on the mapping
names(transformed_df_lsr_)[names(transformed_df_lsr_) %in% names(columns_to_rename)] <- columns_to_rename[names(transformed_df_lsr_)[names(transformed_df_lsr_) %in% names(columns_to_rename)]]

# get the column names containing years using a general regex for 4-digit years
column_names_with_years <- grep("[0-9]{4}", colnames(transformed_df_lsr_), value = TRUE)
# drop column_names_with_years if the year does not begin with a '1' or '2'
column_names_with_years <- column_names_with_years[grepl("_1[0-9]{3}|_2[0-9]{3}|1[0-9]{3}|2[0-9]{3}", column_names_with_years)]
print(column_names_with_years) # parameters with temporal indicator term

# create indicator functions based on the year in the column name and the 'yr2' column
for (col_header in column_names_with_years) {
  # extract the year from the column name
  year_in_col <- as.numeric(sub(".*([0-9]{4}).*", "\\1", col_header))

  # create the indicator variable based on the condition
  indicator_col_name <- paste0("indicator_", col_header)
  transformed_df_lsr_[[indicator_col_name]] <- ifelse(year_in_col <= transformed_df_lsr_$yr2, 1, 0)
}

# drop yr2 from the transformed_df_lsr_ df (no longer need and has not been transformed)
transformed_df_lsr_ <- subset(transformed_df_lsr_, select = -yr2)


# second, add indicator/update indicator for params with missing data
# initialize a vector to store column names with NaN values
columns_with_issues <- c()

# create list of columns to check for nans and nas (ignore lSR and indicator columns)
df_ready_cols <- colnames(transformed_df_lsr_)
cols_to_drop <- c("lSR")
indicator_cols <- grep("indicator", df_ready_cols, value = TRUE)
cols_to_drop <- unique(c(cols_to_drop, indicator_cols))
df_ready_cols <- df_ready_cols[!df_ready_cols %in% cols_to_drop]

# iterate over each column specified in the list
for (col in df_ready_cols) {
  # check if the column contains NaN, Na
  if (any(is.nan(transformed_df_lsr_[[col]])) ||
      any(is.na(transformed_df_lsr_[[col]]))) {
    # if NaN values are found, store the column name
    columns_with_issues <- c(columns_with_issues, col)
  }
}

# print column names with NaN values
print(columns_with_issues)

# create interaction terms for columns with missing data
for (col in columns_with_issues) {
  # create indicator variable for current lc column with missing data
  transformed_df_lsr_[[paste0("indicator2_", col)]] <- as.integer(!is.na(transformed_df_lsr_[[col]]))
}


# third, if multiple indicator functions for a single parameter, combine into a single indicator function:
# identify columns with "indicator_" and "indicator2_"
indicator_cols <- grep("indicator_", colnames(transformed_df_lsr_), value = TRUE)
indicator2_cols <- grep("indicator2_", colnames(transformed_df_lsr_), value = TRUE)
all_indicator_cols <- c(indicator_cols, indicator2_cols)

# print the string after 'indicator_' or 'indicator2'
params <- sub("^.*?_", "", all_indicator_cols)  # Extract substring after first underscore

# if parameter has 2 indicator columns, remove them and create single new indicator columns
unique_params <- unique(params)

# initialize a vector to collect columns to remove from df
# THIS WILL NOT WORK WELL IF PARAMS HAVE SIMILAR WORD IN THEM (I.E., SLOPE)
columns_to_remove <- c()
for (param in unique_params) {
  indicator_cols <- grep(paste0("_", param, "$"), all_indicator_cols, value = TRUE)
  if (length(indicator_cols) > 1) {
    # create new indicator column
    new_indicator_col <- paste0("indicatorc_", param)
    transformed_df_lsr_[[new_indicator_col]] <- ifelse(rowSums(transformed_df_lsr_[indicator_cols]) == 2, 1, 0)
    # collect columns to remove later
    columns_to_remove <- c(columns_to_remove, indicator_cols)
  }
}

# Remove all collected columns at once after the loop
transformed_df_lsr_[columns_to_remove] <- NULL


# now rename all indicator column headers with 'indicator' instead of 'indicator2' and 'indicatorc'
names(transformed_df_lsr_) <- gsub("indicator2|indicatorc", "indicator", names(transformed_df_lsr_))

# list of interaction terms
indicator_cols <- grep("indicator", names(transformed_df_lsr_), value = TRUE)
indicator_cols_trimmed <- sub("^.*?_", "", indicator_cols)  # extract substring after first underscore

# now obtain IVs that don't have an interaction term
IVs <- names(transformed_df_lsr_)[!grepl("indicator", names(transformed_df_lsr_))]
IVs <- setdiff(IVs, indicator_cols_trimmed)
# drop 'lSR' from list since its the dependent variable
IVs <- IVs[-c(1)]
print(IVs) # check to ensure no lakecat/fcpg data in this list, likely these have some missing data


### CREATE MLR MODEL USING ALL POSSIBLE PARAMETERS AND INTERACTION TERMS
valid_indicators <- paste0(indicator_cols_trimmed, ':', indicator_cols)

# make copy of transformed_df_lsr_ prior to changing nans to 0
transformed_df_lsr_nanskept <- transformed_df_lsr_
# convert nans in transformed_df_lsr_ to a number that can be multiplied by the indicator function
transformed_df_lsr_[is.na(transformed_df_lsr_)] <- 0

# construct MLR formula for all data
formula <- as.formula(paste("lSR ~", paste(c(IVs, valid_indicators), collapse = " + ")))
# create MLR model and save model summary to folder
if (!dir.exists("model_evaluation/model_summaries")) {
  dir.create("model_evaluation/model_summaries", recursive = TRUE)
}
filename <- paste0("model_evaluation/model_summaries/lsr_summary_alldata_nopca_", version, ".txt")
sink(filename)
full_model_transform_lsr <- lm(formula, data = transformed_df_lsr_)
print(summary(full_model_transform_lsr))
sink()
print(summary(full_model_transform_lsr))


### PLOT X-Y RELATIONSHIPS
# specify the y parameter
y_parameter_lsr <- "lSR"
# create a list to store plots
plots_lsr <- list()
# create directory to save plot figures in
plot_folder_lsr <- "figures/model_building/X-Y_plots"
dir.create(plot_folder_lsr, showWarnings=FALSE)

# iterate over each x variable
for (col_name in names(transformed_df_lsr_)[-c(1)]) {
  # check if the column name contains 'indicator'
  if (!grepl("indicator", col_name)) {
    # only plot actual data (not 0 or NAN)
    filtered_data <- transformed_df_lsr_[!is.nan(transformed_df_lsr_[[col_name]]) & transformed_df_lsr_[[col_name]] != 0, ]
    cat("Number of rows in", col_name, ":", nrow(filtered_data), "\n")
    # create a scatterplot for the current x variable against the y parameter
    p <- ggplot(filtered_data, aes_string(x = col_name, y = y_parameter_lsr)) +
      geom_point() +
      labs(title = paste("Scatterplot of", col_name, "vs", y_parameter_lsr, "(transformed)")) +
      annotate("text", x = Inf, y = -Inf, label = paste("n =", nrow(filtered_data)),
               hjust = 1, vjust = -1)
    # construct filename based on x and y axes
    filename <- paste(plot_folder_lsr, "/", paste("plot_", paste(col_name, "vs", y_parameter_lsr, "transformed"), ".png", sep = ""), sep = "")
    # save the plot
    ggsave(filename, plot = p, width = 8, height = 6)
    # store the plot in the list
    plots_lsr[[col_name]] <- p
  }
}

# print all plots to screen
print(plots_lsr)


### COMPUTE RELATIVE R2 CONTRIBUTIONS FOR EACH PARAMETER
# fit the initial model with only the response variable
initial_model <- lm(lSR ~ 1, data = transformed_df_lsr_)
# get coefficients from the model
coefficients <- coef(full_model_transform_lsr)
# remove the intercept coefficient
coefficients_no_intercept <- coefficients[-1]  # Remove the first element (intercept)
# initialize a vector to store the relative R^2 contributions
r2_contributions <- numeric(length = length(coefficients_no_intercept))
# iterate through each parameter to be added one by one
for (i in seq_along(coefficients_no_intercept)) {
  # use factors by ensuring they are correctly specified in the formula
  current_vars <- c(IVs, valid_indicators)[1:i]
  # construct the formula with factors properly represented
  formula_r2 <- as.formula(paste("lSR ~", paste(current_vars, collapse = " + ")))
  # fit a new model with the current set of parameters plus the parameter being added
  new_model <- lm(formula_r2, data = transformed_df_lsr_)
  # calculate the change in R^2 between the new model and the initial model
  delta_r2 <- summary(new_model)$r.squared - summary(initial_model)$r.squared
  # store the change in R^2 for the current parameter
  r2_contributions[i] <- delta_r2
}

# print the relative R^2 contributions
print(r2_contributions)
# save contributions to folder
coeff_names <- names(coefficients_no_intercept)
r2_contr <- data.frame(coeff_names, r2_contributions)
r2_contr <- t(r2_contr)
if (!dir.exists("model_evaluation/r2_contributions")) {
  dir.create("model_evaluation/r2_contributions", recursive = TRUE)
}
write.csv(r2_contr, paste("model_evaluation/r2_contributions/coefficients_r2_contributions_alldata_nopca_",version, ".csv"), row.names=FALSE)

# access intercept and coefficients
intercept_usa_nopca <- coef(full_model_transform_lsr)[1]
coefficients_usa_nopca <- coef(full_model_transform_lsr)[-1]
print(intercept_usa_nopca)
print(coefficients_usa_nopca)
# save int and coefs to file
if (!dir.exists("model_evaluation/coefficients_intercepts")) {
  dir.create("model_evaluation/coefficients_intercepts", recursive = TRUE)
}
write.csv(coefficients_usa_nopca, paste("model_evaluation/coefficients_intercepts/coefficients_alldata_nopca_", version, ".csv"), row.names=TRUE)
write.csv(intercept_usa_nopca, paste("model_evaluation/coefficients_intercepts/intercept_alldata_nopca_", version, ".csv"), row.names=TRUE)


### MAKE MODEL EVALUATION PLOTS
# standardized residuals plot
residuals <-rstandard(full_model_transform_lsr)
final_data <- cbind(transformed_df_lsr_, residuals)
final_data[order(-residuals),]
# save to file
file_path <- "figures/model_building/residuals/alldata_nopca"
dir.create(file_path, recursive=TRUE, showWarnings=FALSE)
plot_filename <- paste(file_path, paste("Standardized_Residuals_alldata_nopca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1600, height=1400, res=250, units="px")
plot(final_data$lSR, residuals, ylab='Standardized Residuals', xlab='x', main = paste("Standardized Residuals vs lSR (transformed, all params)"))
abline(0,0)
abline(h=3, col="red", lty=2)
abline(h=-3, col="red", lty=2)
dev.off()

# print to screen
plot(final_data$lSR, residuals, ylab='Standardized Residuals', xlab='x', main = paste("Standardized Residuals vs lSR (transformed, all params)"))
abline(0,0)
abline(h=3, col="red", lty=2)
abline(h=-3, col="red", lty=2)

# number of parameters included in lm model (minus the intercept)
n_params <- length(coef(full_model_transform_lsr))-1

# Q-Q plots (shows normality of residuals)
plot_filename <- paste(file_path, paste("q-q_alldata_nopca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1200, height=1000, res=200)
qqnorm(residuals, main = paste('Q-Q Plot, lSR (transformed, all params),', n_params, "Params"), cex=0.5)
qqline(residuals)
dev.off()

# print to screen
qqnorm(residuals, main = paste('Q-Q Plot, lSR (transformed, all params),', n_params, "Params"), cex=0.5)
qqline(residuals)


# histogram of residuals
plot_filename <- paste(file_path, paste("histogram_alldata_nopca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1400, height=1200, res=200)
hist(residuals, breaks=30, main=paste("Histogram of Residuals, lSR (transformed, all params),", n_params, "Params"), xlab="Residuals", ylab="Frequency")
dev.off()

#print to screen
hist(residuals, breaks=30, main=paste("Histogram of Residuals, lSR (transformed, all params),", n_params, "Params"), xlab="Residuals", ylab="Frequency")




# Principal component analysis (PCA) -------------------------------------------------

### WRITE PCA FUNCTION
# function to create principal components and save variance explained
perform_pca <- function(data, pca_cols, pc_name) {
  # start with new dataframe containing only specified columns for PCA
  pca_df <- subset(data, select = pca_cols)
  # identify rows with complete cases
  complete_cases <- complete.cases(pca_df)
  # scale the data for complete cases
  pca_scaled <- scale(pca_df[complete_cases, ])
  # perform PCA on the scaled complete cases data
  pca_output <- prcomp(pca_scaled, center = TRUE, scale. = TRUE)
  # extract the loadings (coefficients) of the first principal component
  loadings <- pca_output$rotation[, 1]
  # create the principal component scores for complete cases
  pc_scores <- pca_scaled %*% loadings
  # extract the variance explained by the first principal component
  variance_explained <- pca_output$sdev[1]^2 / sum(pca_output$sdev^2)
  # initialize the result vector with NA
  result_vector <- rep(NA, nrow(data))
  # assign the computed PCA scores to the complete cases
  result_vector[complete_cases] <- pc_scores
  
  # convert result to dataframe
  result_df <- data.frame(result_vector)
  colnames(result_df) <- pc_name

  # add the variance explained as an attribute to the dataframe
  attr(result_df, "variance_explained") <- variance_explained

  return(result_df)
}


# make copy of transformed_df_lsr_
transformed_df_lsr_fill <- transformed_df_lsr_nanskept

### PERFORM PCA FOR SIMILAR/COLINEAR PARAMETERS
# Clay, Sand, and OM all around 540 data points
# soil characteristics (pc_soil)
pca_cols_lsr_soil <- c("ClayWs", "SandWs", "PermWs", "OmWs")
data_with_pca_lsr_soil <- perform_pca(transformed_df_lsr_fill, pca_cols_lsr_soil, "pc_soil")
variance_pc_soil <- attr(data_with_pca_lsr_soil, "variance_explained")
print(variance_pc_soil)

# Dam length has 831, dam height has 850
# dam dimensions (pc_damdim)
pca_cols_lsr_damdim <- c("damlength_m", "damH")
data_with_pca_lsr_damdim <- perform_pca(transformed_df_lsr_fill, pca_cols_lsr_damdim, "pc_damdim")
variance_pc_damdim <- attr(data_with_pca_lsr_damdim, "variance_explained")
print(variance_pc_damdim)

# All temp params around 877
# temperature (pc_temp_fcpg)
pca_cols_lsr_temp_fcpg <- c("tmin_norm", "tmax_norm", "NrY_Final")
data_with_pca_lsr_temp_fcpg <- perform_pca(transformed_df_lsr_fill, pca_cols_lsr_temp_fcpg, "pc_temp_fcpg")
variance_pc_temp_fcpg <- attr(data_with_pca_lsr_temp_fcpg, "variance_explained")
print(variance_pc_temp_fcpg)

# CBNFWs has 504, AgKffactWs has 504, Pestic has 527, Fert has 504
# general agriculture (pc_ag)
pca_cols_lsr_ag <- c("CBNFWs", "AgKffactWs", "Pestic1997Ws", "FertWs", "PctCrop2006Ws")
data_with_pca_lsr_ag <- perform_pca(transformed_df_lsr_fill, pca_cols_lsr_ag, "pc_ag")
variance_pc_ag <- attr(data_with_pca_lsr_ag, "variance_explained")
print(variance_pc_ag)

# Superfund has 158, TRIDens has 213, MineDens has 153
# pollution (pc_poll)
pca_cols_lsr_poll <- c("SuperfundDensWs", "TRIDensWs", "MineDensWs", "NPDESDensWs")
data_with_pca_lsr_poll <- perform_pca(transformed_df_lsr_fill, pca_cols_lsr_poll, "pc_poll")
variance_pc_poll <- attr(data_with_pca_lsr_poll, "variance_explained")
print(variance_pc_poll)

# PctImp has 535, PctUrbOp has 538, PctUrbLo has 532, PctUrbMd has 518, PctUrbHi has 454, RdDens has 541, RdCrs has 505, HuDen has 542, PopDen has 541
# removing PctUrbHi since this reduces data points by 50+
# urbanization (pc_urban_2006)
pca_cols_lsr_urban_2006 <- c("PctImp2006Ws", "PctUrbOp2006Ws", "PctUrbLo2006Ws", "PctUrbMd2006Ws", "RdDensWs", "RdCrsWs", "HUDen2010Ws", "PopDen2010Ws")
data_with_pca_lsr_urban_2006 <- perform_pca(transformed_df_lsr_fill, pca_cols_lsr_urban_2006, "pc_urban_2006")
variance_pc_urban_2006 <- attr(data_with_pca_lsr_urban_2006, "variance_explained")
print(variance_pc_urban_2006)

# carbonate residuals / non carbonate residuals (pc_carb)
pca_cols_lsr_carb <- c("PctCarbResidWs", "PctNonCarbResidWs")
data_with_pca_lsr_carb <- perform_pca(transformed_df_lsr_fill, pca_cols_lsr_carb, "pc_carb")
variance_pc_carb <- attr(data_with_pca_lsr_carb, "variance_explained")
print(variance_pc_carb)

# water (pc_water)
pca_cols_lsr_water <- c("SA_m2", "MAQ_NHDcms")
data_with_pca_lsr_water <- perform_pca(transformed_df_lsr_fill, pca_cols_lsr_water, "pc_water")
variance_pc_water <- attr(data_with_pca_lsr_water, "variance_explained")
print(variance_pc_water)




# Correlation matrix with PCs  ----------------------------------------------------

# create custom color palette
my_palette <- colorRampPalette(c("mediumpurple3", "white", "darkorange3"))(100)

# merge pc data into transformed_df_lsr_fill df
transformed_df_lsr_pca <- cbind(transformed_df_lsr_nanskept, 
                                      data_with_pca_lsr_soil, 
                                      data_with_pca_lsr_damdim, 
                                      data_with_pca_lsr_temp_fcpg, 
                                      data_with_pca_lsr_ag, 
                                      data_with_pca_lsr_poll, 
                                      data_with_pca_lsr_urban_2006,
                                      data_with_pca_lsr_water,
                                      data_with_pca_lsr_carb)

# drop indicator columns since will recreate these
cols_to_drop <- names(transformed_df_lsr_pca)[grepl("indicator", names(transformed_df_lsr_pca))]
cols_to_drop <-c(cols_to_drop, "LRR_SYMBOL")
transformed_df_lsr_pca <- transformed_df_lsr_pca[, -which(names(transformed_df_lsr_pca) %in% cols_to_drop)]

### MAKE CORRELATION MATRIX
# specify file path to save fig
file_path <- "figures/model_building/correlation_matrices"
dir.create(file_path, showWarnings=FALSE)

cor_matrix_USA_lsr_pca <- cor(transformed_df_lsr_pca, use = "pairwise.complete.obs")
# save just the first row of the corr matrix (with lSR)
USA_lsr_cors_pca <- cor_matrix_USA_lsr_pca[1, ]
round_cor_matrix_USA_lsr_pca <-round(cor_matrix_USA_lsr_pca, 1)

plot_filename <- paste(file_path, paste("correlation_matrix_alldata_pca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1600, height=1400, res=300, units="px")

plot <- corrplot(round_cor_matrix_USA_lsr_pca,
                 mar=c(0,0,3,0),
                 method="color",
                 type="lower",
                 tl.col="black",
                 tl.srt = 45,
                 addCoef.col = "black",
                 number.cex = 0.15,
                 col = my_palette,
                 tl.cex = 0.29,  # Adjust parameter for smaller font size
                 cl.cex = 0.5)  # Adjust the colorbar font size
mtext(side = 2, text = "Parameters", line = 1.5, cex = 0.65)
mtext(side = 1, text = "r", line = 4, cex = 0.65)
mtext("Correlation Matrix, USA, lSR (transformed, with PCs)", side = 3, line = 0.5, cex = 0.8)
dev.off()


plot <- corrplot(round_cor_matrix_USA_lsr_pca,
                 title = "Correlation Matrix, USA, lSR (transformed, with PCs)",
                 mar=c(0,0,1,0),
                 method="color",
                 type="lower",
                 tl.col="black",
                 tl.srt = 45,
                 addCoef.col = "black",
                 number.cex = 0.4,
                 col = my_palette,
                 tl.cex = 0.6,  # Adjust parameter for smaller font size
                 width = 16,  # Increase plot width
                 height = 16)  # Increase plot height
mtext(side = 2, text = "Parameters", line = -4)

mtext(side = 1, text = "r", line = 3.5)
print(plot)




# MLR model building, lSR (with PCA params) -----------------------------------------------

# remove variables that we made PCs out of
pc_cols_to_drop <- unique(c(pca_cols_lsr_soil, 
                            pca_cols_lsr_damdim, 
                            pca_cols_lsr_temp_fcpg, 
                            pca_cols_lsr_ag, 
                            pca_cols_lsr_poll, 
                            pca_cols_lsr_urban_2006,
                            pca_cols_lsr_water,
                            pca_cols_lsr_carb))
cols_to_keep <- setdiff(names(transformed_df_lsr_pca), pc_cols_to_drop)
# remove the columns using negative indexing
transformed_df_lsr_pca <- transformed_df_lsr_pca[, colnames(transformed_df_lsr_pca) %in% cols_to_keep]

# create copy of pca and transformed dataframe
mlr_usa_lsr_df <- transformed_df_lsr_pca


### CREATE INTERACTION TERMS FOR PARAMETERS WITH MISSING DATA AND/OR ONLY APPLY TEMPORALLY
# first, create interaction terms for temporally-applicable columns (need to add yr2 to dataframe)
# add yr2 column back to dataframe
mlr_usa_lsr_df$yr2 <- filtered_site_data_nonans$yr2
# rename column with only partial years to full years (ignore renaming temp and precip params - assuming similar magnitude over time)
columns_to_rename <- c("Pestic97Ws" = "Pestic1997Ws", "pct_fire_0010Ws" = "pct_fire_decade_2010Ws", "pct_frstloss_0113Ws" = "pct_frstloss_decade_2013Ws")
# rename the columns based on the mapping
names(mlr_usa_lsr_df)[names(mlr_usa_lsr_df) %in% names(columns_to_rename)] <- columns_to_rename[names(mlr_usa_lsr_df)[names(mlr_usa_lsr_df) %in% names(columns_to_rename)]]

# get the column names containing years using a general regex for 4-digit years
column_names_with_years <- grep("[0-9]{4}", colnames(mlr_usa_lsr_df), value = TRUE)
# drop column_names_with_years if the year does not begin with a '1' or '2'
column_names_with_years <- column_names_with_years[grepl("_1[0-9]{3}|_2[0-9]{3}|1[0-9]{3}|2[0-9]{3}", column_names_with_years)]

# create indicator functions based on the year in the column name and the 'yr2' column
for (col_header in column_names_with_years) {
  # extract the year from the column name
  year_in_col <- as.numeric(sub(".*([0-9]{4}).*", "\\1", col_header))
  # create the indicator variable based on the condition
  indicator_col_name <- paste0("indicator_", col_header)
  mlr_usa_lsr_df[[indicator_col_name]] <- ifelse(year_in_col <= mlr_usa_lsr_df$yr2, 1, 0)
}

# print the transformed data frame with indicators
print(mlr_usa_lsr_df)

# drop yr2 from the mlr_usa_lsr_df df (no longer need and has not been transformed)
mlr_usa_lsr_df <- subset(mlr_usa_lsr_df, select = -yr2)


# second, add indicator/update indicator params with missing data
# initialize a vector to store column names with NaN values
columns_with_issues <- c()

# create list of columns to check for nans, na or values >= 0 and <= 1e-6
df_ready_cols <- colnames(mlr_usa_lsr_df)
cols_to_drop <- c("lSR")
indicator_cols <- grep("indicator", df_ready_cols, value = TRUE)
cols_to_drop <- unique(c(cols_to_drop, indicator_cols))
df_ready_cols <- df_ready_cols[!df_ready_cols %in% cols_to_drop]

# iterate over each column specified in the list
for (col in df_ready_cols) {
  # check if the column contains NaN, Na, or values >= 0 and < 1e-6
  if (any(is.nan(mlr_usa_lsr_df[[col]])) ||
      any(is.na(mlr_usa_lsr_df[[col]]))) {
    # if NaN values are found, store the column name
    columns_with_issues <- c(columns_with_issues, col)
  }
}

# print column names with NaN values
print(columns_with_issues)

# create interaction terms for columns with zero data
for (col in columns_with_issues) {
  # create indicator variable for current lc column with missing data
  mlr_usa_lsr_df[[paste0("indicator2_", col)]] <- as.integer(!is.na(mlr_usa_lsr_df[[col]]))
}
print(colnames(mlr_usa_lsr_df))


# third, if multiple indicator functions for a single parameter, combine into a single indicator function:
# identify columns with "indicator_" and "indicator2_"
indicator_cols <- grep("indicator_", colnames(mlr_usa_lsr_df), value = TRUE)
indicator2_cols <- grep("indicator2_", colnames(mlr_usa_lsr_df), value = TRUE)
all_indicator_cols <- c(indicator_cols, indicator2_cols)

# print the string after 'indicator_' or 'indicator2'
params <- sub("^.*?_", "", all_indicator_cols)  # Extract substring after first underscore

# remove duplicates and create new indicator columns
unique_params <- unique(params)
for (param in unique_params) {
  indicator_cols <- grep(paste0("_", param, "$"), all_indicator_cols, value = TRUE)
  if (length(indicator_cols) > 1) {
    # Create new indicator column
    new_indicator_col <- paste0("indicatorc_", param)
    mlr_usa_lsr_df[[new_indicator_col]] <- ifelse(rowSums(mlr_usa_lsr_df[indicator_cols]) == 2, 1, 0)
    # Remove original columns
    mlr_usa_lsr_df[indicator_cols] <- NULL
  }
}
colnames(mlr_usa_lsr_df)

# now rename all indicator column headers with 'indicator' instead of 'indicator2' and 'indicatorc'
names(mlr_usa_lsr_df) <- gsub("indicator2|indicatorc", "indicator", names(mlr_usa_lsr_df))

# list of interaction terms
indicator_cols <- grep("indicator", names(mlr_usa_lsr_df), value = TRUE)
indicator_cols_trimmed <- sub("^.*?_", "", indicator_cols)  # Extract substring after first underscore

# now obtain IVs that don't have an interaction term
IVs <- names(mlr_usa_lsr_df)[!grepl("indicator", names(mlr_usa_lsr_df))]
IVs <- setdiff(IVs, indicator_cols_trimmed)
# drop 'lsr' from list since the dependent variable
IVs <- IVs[-c(1)] # check for accuracy


### CREATE MLR MODEL USING ALL POSSIBLE PARAMETERS (INCLUDING PCS), INTERACTION TERMS
valid_indicators <- paste0(indicator_cols_trimmed, ':', indicator_cols)

# make copy of mlr_usa_lsr_df prior to changing nans to 0
mlr_usa_lsr_df_nanskept <- mlr_usa_lsr_df
# convert nans in mlr_usa_lsr_df to a number that can be multiplied by the indicator function
mlr_usa_lsr_df[is.na(mlr_usa_lsr_df)] <- 0

# construct MLR formula with all data (+ pcs)
formula <- as.formula(paste("lSR ~", paste(c(IVs, valid_indicators), collapse = " + ")))
# create MLR model and save model summary to folder
filename <- paste0("model_evaluation/model_summaries/lsr_summary_alldata_pca_", version, ".txt")
sink(filename)
full_model_transform_pca_lsr <- lm(formula, data = mlr_usa_lsr_df)
print(summary(full_model_transform_pca_lsr))
sink()
print(summary(full_model_transform_pca_lsr))


### COMPUTE RELATIVE R2 CONTRIBUTIONS FOR EACH PARAMETER
# fit the initial model with only the response variable
initial_model <- lm(lSR ~ 1, data = mlr_usa_lsr_df)
# get coefficients from the model
coefficients <- coef(full_model_transform_pca_lsr)
# remove the intercept coefficient
coefficients_no_intercept <- coefficients[-1]  # Remove the first element (intercept)
# initialize a vector to store the relative R^2 contributions
r2_contributions <- numeric(length = length(coefficients_no_intercept))
# iterate through each parameter to be added one by one
for (i in seq_along(coefficients_no_intercept)) {
  # fit a new model with the current set of parameters plus the parameter being added
  formula_r2 <- as.formula(paste("lSR ~", paste(c(IVs,valid_indicators)[1:i],
                                                collapse = " + ")))
  new_model <- lm(formula_r2, data = mlr_usa_lsr_df)
  # calculate the change in R^2 between the new model and the initial model
  delta_r2 <- summary(new_model)$r.squared - summary(initial_model)$r.squared
  # store the change in R^2 for the current parameter
  r2_contributions[i] <- delta_r2
}

# print the relative R^2 contributions
print(r2_contributions)
# save contributions to folder
coeff_names <- names(coefficients_no_intercept)
r2_contr <- data.frame(coeff_names, r2_contributions)
r2_contr <- t(r2_contr)
write.csv(r2_contr, paste("model_evaluation/r2_contributions/coefficients_r2_contributions_alldata_pca_",version, ".csv"), row.names=FALSE)

# access intercept and coefficients
intercept_usa_pca <- coef(full_model_transform_pca_lsr)[1]
coefficients_usa_pca <- coef(full_model_transform_pca_lsr)[-1]
print(intercept_usa_pca)
print(coefficients_usa_pca)
# save int and coefs to file
write.csv(coefficients_usa_pca, paste("model_evaluation/coefficients_intercepts/coefficients_alldata_pca_", version, ".csv"), row.names=TRUE)
write.csv(intercept_usa_pca, paste("model_evaluation/coefficients_intercepts/intercept_alldata_pca_", version, ".csv"), row.names=TRUE)


### MAKE MODEL EVALUATION PLOTS
# standardized residuals plot
residuals <-rstandard(full_model_transform_pca_lsr)
final_data <- cbind(mlr_usa_lsr_df, residuals)
final_data[order(-residuals),]

# save to file
file_path <- "figures/model_building/residuals/alldata_pca"
dir.create(file_path, showWarnings=FALSE)

plot_filename <- paste(file_path, paste("Standardized_Residuals_alldata_pca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1600, height=1400, res=250, units="px")
plot(final_data$lSR, residuals, ylab='Standardized Residuals', xlab='x', main = paste("Standardized Residuals vs lSR (transformed, all params)"))
abline(0,0)
abline(h=3, col="red", lty=2)
abline(h=-3, col="red", lty=2)
dev.off()

# print to screen
plot(final_data$lSR, residuals, ylab='Standardized Residuals', xlab='x', main = paste("Standardized Residuals vs lSR (transformed, all params)"))
abline(0,0)
abline(h=3, col="red", lty=2)
abline(h=-3, col="red", lty=2)

# number of parameters included in lm model (minus the intercept)
n_params <- length(coef(full_model_transform_pca_lsr))-1


# Q-Q plots (shows normality of residuals)
# Calculate residuals and their variances
residuals <- resid(full_model_transform_pca_lsr)

plot_filename <- paste(file_path, paste("q-q_allparams_pca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1200, height=1000, res=200)
qqnorm(residuals, main = paste('Q-Q Plot, lSR (transformed, all params),', n_params, "Params"), cex=0.5)
qqline(residuals)
dev.off()

# print to screen
qqnorm(residuals, main = paste('Q-Q Plot, lSR (transformed, all params),', n_params, "Params"), cex=0.5)
qqline(residuals)


# histogram of residuals
plot_filename <- paste(file_path, paste("histogram_allparams_pca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1400, height=1200, res=200)
hist(residuals, breaks=30, main=paste("Histogram of Residuals, lSR (transformed, all params),", n_params, "Params"), xlab="Residuals", ylab="Frequency")
dev.off()

#print to pane
hist(residuals, breaks=30, main=paste("Histogram of Residuals, lSR (transformed, all params),", n_params, "Params"), xlab="Residuals", ylab="Frequency")



### CREATE MLR MODEL USING BEST, SELECTED PARAMETERS (NO LIMIT ON NUMBER)
no_colinearities <- NULL
force_in <- c("yrs_tot")
colinearities <- c("pc_water:indicator_pc_water") # update with any identified colinear variables
# perform subset selection
subset_USA_lsr <-
  regsubsets(formula,
             data = mlr_usa_lsr_df,
             nvmax = NULL,    # NULL for no limit on number of variables
             force.in = force_in, force.out = colinearities,
             method = "forward")


# summary of subset selection
subset_summary <- summary(subset_USA_lsr)

# choose the best subset based on different criteria
best_adjr2 <- which.max(subset_summary$adjr2)
best_aic <- which.min(subset_summary$aic)
best_bic <- which.min(subset_summary$bic)
best_cp <- which.min(subset_summary$cp)

# select predictors for the best models based on each criterion
selected_adjr2 <- names(coef(subset_USA_lsr, id = best_adjr2))[-1]
#selected_aic <- names(coef(subset_USA_lsr, id = best_aic))[-1]
selected_bic <- names(coef(subset_USA_lsr, id = best_bic))[-1]
selected_cp <- names(coef(subset_USA_lsr, id = best_cp))[-1]

# choose best subset based on adjusted R-squared
best_subset_USA_lsr <- best_bic
selected_predictors_USA_lsr_null_pca <- names(coef(subset_USA_lsr, id=best_subset_USA_lsr))[-1]    # print names of best variables
selected_predictors_USA_lsr_null_pca <- c("lSR", selected_predictors_USA_lsr_null_pca)     # adds lsr term
selected_predictors_USA_lsr_null_pca <- gsub("\"", "", selected_predictors_USA_lsr_null_pca)    # removes quotation marks which impede parsing
print(selected_predictors_USA_lsr_null_pca)


# if an interaction term is included in the list of selected predictors, then i need to split the terms
updated_selected_predictors_USA <- list()

# Iterate through each value in selected_predictors_USA_lsr_null_pca
for (value in selected_predictors_USA_lsr_null_pca) {
  # Check if the value contains a colon
  if (grepl(":", value)) {
    # If a colon is found, split the value into two values
    split_values <- unlist(strsplit(value, ":", fixed = TRUE))
    # Add the split values to the updated list
    updated_selected_predictors_USA <- c(updated_selected_predictors_USA, split_values)
  } else {
    # If no colon is found, add the value to the updated list without modification
    updated_selected_predictors_USA <- c(updated_selected_predictors_USA, value)
  }
}

updated_selected_predictors_USA <- gsub("\"", "", updated_selected_predictors_USA)    # removes quotation marks which impede parsing
print(updated_selected_predictors_USA)

# create new MLR model with only selected predictors and interaction terms
mlr_usa_lsr_df_best_null_pca <- subset(mlr_usa_lsr_df, select=updated_selected_predictors_USA)

# check for any colinearity, if so then need to remove colinear parameters before building model
cols_to_include_lsr <- names(mlr_usa_lsr_df_best_null_pca)[!grepl("indicator", names(mlr_usa_lsr_df_best_null_pca))]

# make zeros nan so will be ignored by corr matrix
mlr_usa_lsr_df_best_null_pca[mlr_usa_lsr_df_best_null_pca == 0] <- NA

# generate corr matrix
cor_matrix_USA_lsr <- cor(mlr_usa_lsr_df_best_null_pca[cols_to_include_lsr], use = "pairwise.complete.obs")
round_cor_matrix_USA_lsr <-round(cor_matrix_USA_lsr, 2)

file_path <- "figures/model_building/correlation_matrices"
dir.create(file_path, showWarnings=FALSE)

# number of parameters included in lm model (minus the response term)
n_params <- updated_selected_predictors_USA[!grepl("indicator", updated_selected_predictors_USA)]
n_params <- length(n_params)-1

plot_filename <- paste(file_path, paste("correlation_matrix_subset_selected_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1600, height=1400, res=300, units="px")
plot <- corrplot(round_cor_matrix_USA_lsr,
                 mar=c(1,0,2,0),
                 method="color",
                 type="lower",
                 tl.col="black",
                 tl.srt = 45,
                 tl.offset = 0.75, 
                 addCoef.col = "black",
                 number.cex = .5,
                 col = my_palette,
                 tl.cex = 0.6,  # Adjust parameter for smaller font size
                 cl.cex = 0.6)  # Adjust the colorbar font size
mtext(side = 2, text = "Parameters", line = 2.5, cex = 0.9)
mtext(side = 1, text = "r", line = 3.4, cex = 0.9)
mtext(paste("Correlation Matrix, USA, lSR (transformed, with PCs), Best", n_params, "Params"), side = 3, line = 0.8, cex = 1)
dev.off()


plot <- corrplot(round_cor_matrix_USA_lsr,
                 mar=c(0,0,2,0),
                 method="color",
                 type="lower",
                 tl.col="black",
                 tl.srt = 45,
                 tl.offset = 0.7, 
                 addCoef.col = "black",
                 number.cex = 0.7,
                 col = my_palette,
                 tl.cex = 1,  # Adjust parameter for smaller font size
                 width = 16,  # Increase plot width
                 height = 16)  # Increase plot height
mtext(side = 2, text = "Parameters", line = 2.5)
mtext(side = 1, text = "r", line = 3.5)
mtext(paste("Correlation Matrix, USA, lSR (transformed, with PCs), Best", n_params, "Params"), side = 3, line = -.8, cex = .8)
print(plot)


# make list of IVs to include in model
IVs <- updated_selected_predictors_USA
IVs <- IVs[!grepl("indicator", IVs)]
IVs <- IVs[!IVs %in% params]
IVs <- IVs[-1]
interaction_terms <- c()

# iterate through each value in selected_predictors_USA_lsr_null_pca
for (value in selected_predictors_USA_lsr_null_pca) {
  # check if the value contains a colon
  if (str_detect(value, ":")) {
    # if a colon is found, save the value as an interaction term
    interaction_terms <- c(interaction_terms, value)
  }
}
print(interaction_terms)

# replace all nan values with 0, so it can be multiplied by indicator params
mlr_usa_lsr_df_best_null_pca[is.na(mlr_usa_lsr_df_best_null_pca)] <- 0

formula_null <- as.formula(paste("lSR ~", paste(c(IVs, interaction_terms), collapse = " + ")))
# create MLR model and save model summary to folder
filename <- paste0("model_evaluation/model_summaries/lsr_summary_subsetdata_pca_", version, ".txt")
sink(filename)
mlr_usa_lsr_df_best_model_null_pca <- lm(formula_null, data = mlr_usa_lsr_df_best_null_pca)
print(summary(mlr_usa_lsr_df_best_model_null_pca))
sink()
print(summary(mlr_usa_lsr_df_best_model_null_pca))


### COMPUTE RELATIVE R2 CONTRIBUTIONS FOR EACH PARAMETER
# fit the initial model with only the response variable
initial_model <- lm(lSR ~ 1, data = mlr_usa_lsr_df_best_null_pca)
# get coefficients from the model
coefficients <- coef(mlr_usa_lsr_df_best_model_null_pca)
# remove the intercept coefficient
coefficients_no_intercept <- coefficients[-1]  # Remove the first element (intercept)
# initialize a vector to store the relative R^2 contributions
r2_contributions <- numeric(length = length(coefficients_no_intercept))
# iterate through each parameter to be added one by one
for (i in seq_along(coefficients_no_intercept)) {
  # fit a new model with the current set of parameters plus the parameter being added
  formula_r2 <- as.formula(paste("lSR ~", paste(c(IVs,interaction_terms)[1:i],collapse = " + ")))
  new_model <- lm(formula_r2, data = mlr_usa_lsr_df_best_null_pca)
  # Calculate the change in R^2 between the new model and the initial model
  delta_r2 <- summary(new_model)$r.squared - summary(initial_model)$r.squared
  # Store the change in R^2 for the current parameter
  r2_contributions[i] <- delta_r2
}

# print the relative R^2 contributions
print(r2_contributions)
# save contributions to folder
coeff_names <- names(coefficients_no_intercept)
r2_contr <- data.frame(coeff_names, r2_contributions)
r2_contr <- t(r2_contr)
write.csv(r2_contr, paste("model_evaluation/r2_contributions/coefficients__r2_contributions_subsetdata_pca_",version, ".csv"), row.names=FALSE)

# access intercept and coefficients
intercept_usa_null_pca <- coef(mlr_usa_lsr_df_best_model_null_pca)[1]
coefficients_usa_null_pca <- coef(mlr_usa_lsr_df_best_model_null_pca)[-1]
print(intercept_usa_null_pca)
print(coefficients_usa_null_pca)
# save int and coefs to file
write.csv(coefficients_usa_null_pca, paste("model_evaluation/coefficients_intercepts/coefficients_subsetdata_pca_", version, ".csv"), row.names=TRUE)
write.csv(intercept_usa_null_pca, paste("model_evaluation/coefficients_intercepts/intercept_subsetdata_pca_", version, ".csv"), row.names=TRUE)


### MAKE MODEL EVALUATION PLOTS
# standardized residuals plot
residuals <-rstandard(mlr_usa_lsr_df_best_model_null_pca)
final_data <- cbind(mlr_usa_lsr_df_best_null_pca, residuals)
final_data[order(-residuals),]

# save to file
file_path <- "figures/model_building/residuals/subsetdata_pca"
dir.create(file_path, showWarnings=FALSE)

plot_filename <- paste(file_path, paste("Standardized_Residuals_subsetdata_pca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1600, height=1400, res=250, units="px")
plot(final_data$lSR, residuals, ylab='Standardized Residuals', xlab='x', main = paste("Standardized Residuals vs lSR (transformed, select params)"))
abline(0,0)
abline(h=3, col="red", lty=2)
abline(h=-3, col="red", lty=2)
dev.off()

# print to screen
plot(final_data$lSR, residuals, ylab='Standardized Residuals', xlab='x', main = paste("Standardized Residuals vs lSR (transformed, select params)"))
abline(0,0)
abline(h=3, col="red", lty=2)
abline(h=-3, col="red", lty=2)

# number of parameters included in lm model (minus the intercept)
n_params <- length(coef(mlr_usa_lsr_df_best_model_null_pca))-1


# Q-Q plots (shows normality of residuals)
# Calculate residuals and their variances
residuals <- resid(mlr_usa_lsr_df_best_model_null_pca)

plot_filename <- paste(file_path, paste("q-q_subsetparams_pca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1200, height=1000, res=200)
qqnorm(residuals, main = paste('Q-Q Plot, lSR (transformed, select params),', n_params, "Params"), cex=0.5)
qqline(residuals)
dev.off()

# print to screen
qqnorm(residuals, main = paste('Q-Q Plot, lSR (transformed, select params),', n_params, "Params"), cex=0.5)
qqline(residuals)


# histogram of residuals
plot_filename <- paste(file_path, paste("histogram_subsetparams_pca_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1400, height=1200, res=200)
hist(residuals, breaks=30, main=paste("Histogram of Residuals, lSR (transformed, select params),", n_params, "Params"), xlab="Residuals", ylab="Frequency")
dev.off()

#print to pane
hist(residuals, breaks=30, main=paste("Histogram of Residuals, lSR (transformed, select params),", n_params, "Params"), xlab="Residuals", ylab="Frequency")


### FINAL MLR MODEL WITH JUST MAIN 4 PARAMETERS THAT TOGETHER EXPLAIN OVER 80% OF THE VARIANCE IN LSR
if (version == "rattes_v1p2_3param_silt") {
  selected_predictors_USA_lsr_4 <- selected_predictors_USA_lsr_null_pca[
    selected_predictors_USA_lsr_null_pca %in% c("lSR", "AveTrap", "wSAatdam", "fdar:indicator_fdar")
  ]
} else {
  selected_predictors_USA_lsr_4 <- selected_predictors_USA_lsr_null_pca[
    selected_predictors_USA_lsr_null_pca %in% c("lSR", "yrs_tot", "AveTrap", "wSAatdam", "fdar:indicator_fdar")
  ]
}
print(selected_predictors_USA_lsr_4)

# if an interaction term is included in the list of selected predictors, then i need to split the terms
updated_selected_predictors_USA <- list()

# Iterate through each value in selected_predictors_USA_lsr_4
for (value in selected_predictors_USA_lsr_4) {
  # Check if the value contains a colon
  if (grepl(":", value)) {
    # If a colon is found, split the value into two values
    split_values <- unlist(strsplit(value, ":", fixed = TRUE))
    # Add the split values to the updated list
    updated_selected_predictors_USA <- c(updated_selected_predictors_USA, split_values)
  } else {
    # If no colon is found, add the value to the updated list without modification
    updated_selected_predictors_USA <- c(updated_selected_predictors_USA, value)
  }
}

updated_selected_predictors_USA <- gsub("\"", "", updated_selected_predictors_USA)    # removes quotation marks which impede parsing
print(updated_selected_predictors_USA)

# create new MLR model with only selected predictors and interaction terms
mlr_usa_lsr_df_best_4 <- subset(mlr_usa_lsr_df, select=updated_selected_predictors_USA)

# check for any colinearity, if so then need to remove colinear parameters before building model
cols_to_include_lsr <- names(mlr_usa_lsr_df_best_4)[!grepl("indicator", names(mlr_usa_lsr_df_best_4))]

# make zeros nan so will be ignored by corr matrix
mlr_usa_lsr_df_best_4[mlr_usa_lsr_df_best_4 == 0] <- NA

# generate corr matrix
cor_matrix_USA_lsr <- cor(mlr_usa_lsr_df_best_4[cols_to_include_lsr], use = "pairwise.complete.obs")
round_cor_matrix_USA_lsr <-round(cor_matrix_USA_lsr, 2)

file_path <- "figures/model_building/correlation_matrices"
dir.create(file_path, showWarnings=FALSE)

# number of parameters included in lm model (minus the response term)
n_params <- updated_selected_predictors_USA[!grepl("indicator", updated_selected_predictors_USA)]
n_params <- length(n_params)-1

plot_filename <- paste(file_path, paste("correlation_matrix_subset_best4_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1600, height=1400, res=300, units="px")
plot <- corrplot(round_cor_matrix_USA_lsr,
                 mar=c(1,0,2,0),
                 method="color",
                 type="lower",
                 tl.col="black",
                 tl.srt = 45,
                 tl.offset = 0.75, 
                 addCoef.col = "black",
                 number.cex = .5,
                 col = my_palette,
                 tl.cex = 0.6,  # Adjust parameter for smaller font size
                 cl.cex = 0.6)  # Adjust the colorbar font size
mtext(side = 2, text = "Parameters", line = 1.2, cex = 0.9)
mtext(side = 1, text = "r", line = 3.4, cex = 0.9)
mtext(paste("Correlation Matrix, USA, lSR (transformed, with PCs), Best", n_params, "Params"), side = 3, line = 0.8, cex = 1)
dev.off()

plot <- corrplot(round_cor_matrix_USA_lsr,
                 mar=c(0,0,2,0),
                 method="color",
                 type="lower",
                 tl.col="black",
                 tl.srt = 45,
                 tl.offset = 0.7, 
                 addCoef.col = "black",
                 number.cex = 0.7,
                 col = my_palette,
                 tl.cex = 1,  # Adjust parameter for smaller font size
                 width = 16,  # Increase plot width
                 height = 16)  # Increase plot height
mtext(side = 2, text = "Parameters", line = 2)
mtext(side = 1, text = "r", line = 3.5)
mtext(paste("Correlation Matrix, USA, lSR (transformed, with PCs), Best", n_params, "Params"), side = 3, line = -.8, cex = .8)
print(plot)


# make list of IVs to include in model
IVs <- updated_selected_predictors_USA
IVs <- IVs[!grepl("indicator", IVs)]
IVs <- IVs[!IVs %in% params]
IVs <- IVs[-1]
interaction_terms <- c()

# iterate through each value in selected_predictors_USA_lsr_null_pca
for (value in selected_predictors_USA_lsr_4) {
  # check if the value contains a colon
  if (str_detect(value, ":")) {
    # if a colon is found, save the value as an interaction term
    interaction_terms <- c(interaction_terms, value)
  }
}
print(interaction_terms)

# replace all nan values with 0, so it can be multiplied by indicator params
mlr_usa_lsr_df_best_4[is.na(mlr_usa_lsr_df_best_4)] <- 0

formula_4 <- as.formula(paste("lSR ~", paste(c(IVs, interaction_terms), collapse = " + ")))
# create MLR model and save model summary to folder
filename <- paste0("model_evaluation/model_summaries/lsr_summary_subset_best4_", version, ".txt")
sink(filename)
mlr_usa_lsr_df_best_4_model <- lm(formula_4, data = mlr_usa_lsr_df_best_4)
print(summary(mlr_usa_lsr_df_best_4_model))
sink()
print(summary(mlr_usa_lsr_df_best_4_model))


### COMPUTE RELATIVE R2 CONTRIBUTIONS FOR EACH PARAMETER
# fit the initial model with only the response variable
initial_model <- lm(lSR ~ 1, data = mlr_usa_lsr_df_best_4)
# get coefficients from the model
coefficients <- coef(mlr_usa_lsr_df_best_4_model)
# remove the intercept coefficient
coefficients_no_intercept <- coefficients[-1]  # Remove the first element (intercept)
# initialize a vector to store the relative R^2 contributions
r2_contributions <- numeric(length = length(coefficients_no_intercept))
# iterate through each parameter to be added one by one
for (i in seq_along(coefficients_no_intercept)) {
  # fit a new model with the current set of parameters plus the parameter being added
  formula_r2 <- as.formula(paste("lSR ~", paste(c(IVs,interaction_terms)[1:i],collapse = " + ")))
  new_model <- lm(formula_r2, data = mlr_usa_lsr_df_best_4)
  # Calculate the change in R^2 between the new model and the initial model
  delta_r2 <- summary(new_model)$r.squared - summary(initial_model)$r.squared
  # Store the change in R^2 for the current parameter
  r2_contributions[i] <- delta_r2
}

# print the relative R^2 contributions
print(r2_contributions)
# save contributions to folder
coeff_names <- names(coefficients_no_intercept)
r2_contr <- data.frame(coeff_names, r2_contributions)
r2_contr <- t(r2_contr)
write.csv(r2_contr, paste("model_evaluation/r2_contributions/coefficients__r2_contributions_subset_best4_",version, ".csv"), row.names=FALSE)

# access intercept and coefficients
intercept_usa_best4 <- coef(mlr_usa_lsr_df_best_4_model)[1]
coefficients_usa_best4 <- coef(mlr_usa_lsr_df_best_4_model)[-1]
print(intercept_usa_best4)
print(coefficients_usa_best4)
# save int and coefs to file
write.csv(coefficients_usa_null_pca, paste("model_evaluation/coefficients_intercepts/coefficients_subset_best4_", version, ".csv"), row.names=TRUE)
write.csv(intercept_usa_null_pca, paste("model_evaluation/coefficients_intercepts/intercept_subset_best4_", version, ".csv"), row.names=TRUE)


# NEW CONFIDENCE INTERVAL CODE BELOW
coef_beta <- coef(mlr_usa_lsr_df_best_4_model)
vcov_beta <- vcov(mlr_usa_lsr_df_best_4_model)
rse <- summary(mlr_usa_lsr_df_best_4_model)$sigma
dof <- summary(mlr_usa_lsr_df_best_4_model)$df[2]
t_critical_95 <- qt(0.975, df = dof)
# End new code


### MAKE MODEL EVALUATION PLOTS
# number of parameters included in lm model (minus the intercept)
n_params <- length(coef(mlr_usa_lsr_df_best_4_model))-1

# standardized residuals
residuals <-rstandard(mlr_usa_lsr_df_best_4_model)
final_data <- cbind(mlr_usa_lsr_df_best_4, residuals)
final_data[order(-residuals),]

# save to file
file_path <- "figures/model_building/residuals/subset_best4"
dir.create(file_path, showWarnings=FALSE)

plot_filename <- paste(file_path, paste("Standardized_Residuals_subset_best4_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1400, height=1200, res=250, units="px")
plot(final_data$lSR, residuals, ylab='Standardized Residuals', xlab='lSR', main = paste("Standardized Residuals,", n_params, "param(s) vs lSR (transformed)"), cex=0.75, cex.main=0.8, cex.axis=0.85)
abline(0,0)
abline(h=3, col="red", lty=2)
abline(h=-3, col="red", lty=2)
dev.off()

# print to screen
plot(final_data$lSR, residuals, ylab='Standardized Residuals', xlab='lSR', main = paste("Standardized Residuals,", n_params, "param(s) vs lSR (transformed)"))
abline(0,0)
abline(h=3, col="red", lty=2)
abline(h=-3, col="red", lty=2)


# Q-Q plots (shows normality of residuals)
# Calculate residuals and their variances
residuals <- resid(mlr_usa_lsr_df_best_4_model)

plot_filename <- paste(file_path, paste("q-q_subset_best4_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1200, height=1000, res=250)
qqnorm(residuals, main = paste('Q-Q Plot, lSR (transformed),', n_params, 'Param(s)'), cex=0.5, cex.main=0.8, cex.axis=0.7, cex.lab=0.8)
qqline(residuals)
dev.off()

# print to screen
qqnorm(residuals, main = paste('Q-Q Plot, lSR (transformed),', n_params, 'Param(s)'), cex=0.5)
qqline(residuals)

# normal residuals plots
# create directory to save plot figures in
plot_folder_lsr <- "figures/model_building/residuals/subset_best4/residuals"
dir.create(plot_folder_lsr, showWarnings=FALSE)

# create a list to store plots
plots <- list()

#  fit linear regression models and create residual plots for each x-y relationship
for (col_name in names(mlr_usa_lsr_df_best_4)[-which(names(mlr_usa_lsr_df_best_4) == "lSR")]) {
  # check if 'indicator' is not in col_name
  if (!grepl("indicator", col_name)) {
    # remove NAs from the current column
    mlr_usa_lsr_df_best_4 <- mlr_usa_lsr_df_best_4[!is.na(mlr_usa_lsr_df_best_4[[col_name]]) & mlr_usa_lsr_df_best_4[[col_name]] != 0, ]
    # fit linear regression model
    model <- lm(lSR ~ ., data = mlr_usa_lsr_df_best_4[, c("lSR", col_name), drop = FALSE])
    # calculate residuals
    residuals <- resid(model)
    # create residual plot
    plot_filename <- paste(plot_folder_lsr, paste("Residual_Plot_", col_name, "_lsr_usa_best4_", version, ".png", sep = ""), sep = "/")
    png(plot_filename, width = 1200, height = 1000, res = 200)
    plot(mlr_usa_lsr_df_best_4[[col_name]], residuals,
         main = paste("Residual Plot for", col_name, "vs lSR (transformed)"),
         xlab = col_name,
         ylab = "Residuals")
    abline(h = 0, col = "red", lty = 2)  # Add a horizontal line at y = 0
    dev.off()

    # print the plot to plot pane
    print(plot(mlr_usa_lsr_df_best_4[[col_name]], residuals,
               main = paste("Residual Plot for", col_name, "vs lSR (transformed)"),
               xlab = col_name,
               ylab = "Residuals"))
    abline(h = 0, col = "red", lty = 2)  # Add a horizontal line at y = 0
  }
}

# histogram of residuals
file_path <- "figures/model_building/residuals/subset_best4"
dir.create(file_path, showWarnings=FALSE)

residuals <- resid(mlr_usa_lsr_df_best_4_model)

plot_filename <- paste(file_path, paste("histogram_subset_best4_", version, ".png", sep=""), sep="/")
png(plot_filename, width=1500, height=1400, res=250)
hist(residuals, breaks=30, main=paste("Histogram of Residuals, lSR (transformed),", n_params, "Param(s)"), xlab="Residuals", ylab="Frequency" )
dev.off()

#print to pane
hist(residuals, breaks=30, main=paste("Histogram of Residuals, lSR (transformed),", n_params, "Param(s)"), xlab="Residuals", ylab="Frequency")




# MLR model evaluation, lSR ----------------------------------------------------------


### EVALUATE MODEL WITH ALL DATA AND PCA
# create testing and training sets to evaluate model performance for USA model with all params
file_path <- "figures/model_building/model_evaluation"
dir.create(file_path, showWarnings=FALSE)

# split data into 70% training and 30% test
# make sure no nans in df
mlr_usa_lsr_df[is.na(mlr_usa_lsr_df)] <- 0
set.seed(121)  # for reproducibility
train_indices <- createDataPartition(mlr_usa_lsr_df$lSR, p = 0.7, list = FALSE)
train_data <- mlr_usa_lsr_df[train_indices, ]
test_data <- mlr_usa_lsr_df[-train_indices, ]

# train linear regression model
lm_model <- lm(formula, data = train_data)
print(summary(lm_model))

# predict on test data
predictions <- predict(lm_model, newdata = test_data)

# number of parameters included in lm model (minus the intercept)
n_params <- length(coef(lm_model))-1

# plot actual versus predicted values
# save to file
plot_filename <- paste(file_path, paste("actual_vs_predicted_1-1_alldata_pca.png", sep=""), sep="/")
png(plot_filename,width=1400, height=1400, res=250)
# Remove the default grey background
par(bg = "white")
plot(test_data$lSR, predictions,
     main = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
     xlab = "Actual Values", ylab = "Predicted Values",
     cex.axis=0.95,
     cex.lab=1.1,
     col.lab='black',
     xlim = c(4, 18),   # Set x-axis limits
     ylim = c(4, 18))   # Set y-axis limits to be the same as x-axis
grid(nx = NULL, ny = NULL, col = "grey80", lty = "dashed")
abline(0, 1, col = "red")  # Add a 1:1 line
dev.off()
# print to screen
plot(test_data$lSR, predictions,
     main = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
     xlab = "Actual Values", ylab = "Predicted Values")
grid(nx = NULL, ny = NULL, col = "grey80", lty = "dashed")
abline(0, 1, col = "red")  # Add a 1:1 line

# Combine actual and predicted values into a data frame
results <- data.frame(Actual = test_data$lSR, Predicted = predictions)

# create histograms with density distributions
# save to file
plot_filename <- paste(file_path, paste("actual_vs_predicted_distplot_alldata_pca.png", sep=""), sep="/")
png(plot_filename, width=1500, height=1400, res=250)
par(bg = "white")
ggplot(results, aes(x = Actual, fill = "Actual")) +
  geom_density(alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  labs(title = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
       x = "lSR", y = "Density") +
  scale_fill_manual(values = c("Actual" = "azure4", "Predicted" = "cadetblue1")) +
  guides(fill = guide_legend(title = ""))+
  theme(
    panel.background = element_rect(fill = "white", color = "black"), 
    plot.background = element_rect(fill="white", color=NA),
    panel.grid.major = element_line(color = "grey80", linetype="dashed"),  # Darker major grid lines
    #  panel.grid.minor = element_blank,  # Darker minor grid lines
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text = element_text(color = "black", size=11),  # Ensure axis text is black for visibility
    axis.title = element_text(color = "black", size=13),  # Ensure axis title is black for visibility
    legend.position = "bottom"
  )
dev.off()
# print plot so screen
ggplot(results, aes(x = Actual, fill = "Actual")) +
  geom_density(alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  labs(title = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
       x = "lSR", y = "Density") +
  scale_fill_manual(values = c("Actual" = "darkmagenta", "Predicted" = "darkorange1")) +
  guides(fill = guide_legend(title = ""))
theme(
  panel.grid.major = element_line(color = "grey40"),  # Darker major grid lines
  panel.grid.minor = element_line(color = "grey60"),  # Darker minor grid lines
  axis.line = element_line(color = "black"),  # Add axis lines
  axis.text = element_text(color = "black"),  # Ensure axis text is black for visibility
  axis.title = element_text(color = "black")  # Ensure axis title is black for visibility
)


# compute model evaluation statistics
mse <- mean((predictions - test_data$lSR)^2)
rmse <- sqrt(mse)
mae <- mean(abs(predictions - test_data$lSR))
rsquared <- summary(lm_model)$r.squared
adjusted_r2 <- summary(lm_model)$adj.r.squared
model_summ <- summary(lm_model)
rse <- model_summ$sigma

# print model evaluation statistics
cat("Mean Squared Error (MSE):", mse, "\n")
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")
cat("Residual Standard Error (RSE):", rse, "\n")
cat("R-squared:", rsquared, "\n")
cat("Adjusted R-squared:", adjusted_r2, "\n")

# create dataframe with model eval results
model_eval_mlr_usa_lsr_df_all <- data.frame(
  which_model = 'USA-All',
  MSE = mse,
  RMSE = rmse,
  MAE = mae,
  RSE = rse,
  R_squared = rsquared,
  Adj_R_squared = adjusted_r2
)


### EVALUATE MODEL WITH SUBSET PARAMETER DATA AND PCA
# split data into 70% training and 30% test
# make sure no nans in df
mlr_usa_lsr_df_best_null_pca[is.na(mlr_usa_lsr_df_best_null_pca)] <- 0
set.seed(121)  # for reproducibility
train_indices <- createDataPartition(mlr_usa_lsr_df_best_null_pca$lSR, p = 0.7, list = FALSE)
train_data <- mlr_usa_lsr_df_best_null_pca[train_indices, ]
test_data <- mlr_usa_lsr_df_best_null_pca[-train_indices, ]

# train linear regression model
lm_model <- lm(formula_null, data = train_data)
print(summary(lm_model))

# predict on test data
predictions <- predict(lm_model, newdata = test_data)

# number of parameters included in lm model (minus the intercept)
n_params <- length(coef(lm_model))-1

# plot actual versus predicted values
# save to file
plot_filename <- paste(file_path, paste("actual_vs_predicted_1-1_subsetdata_pca.png", sep=""), sep="/")
png(plot_filename,width=1400, height=1400, res=250)
# Remove the default grey background
par(bg = "white")
plot(test_data$lSR, predictions,
     main = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
     xlab = "Actual Values", ylab = "Predicted Values",
     cex.axis=0.95,
     cex.lab=1.1,
     col.lab='black',
     xlim = c(4, 18),   # Set x-axis limits
     ylim = c(4, 18))   # Set y-axis limits to be the same as x-axis
grid(nx = NULL, ny = NULL, col = "grey80", lty = "dashed")
abline(0, 1, col = "red")  # Add a 1:1 line
dev.off()
# print to screen
plot(test_data$lSR, predictions,
     main = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
     xlab = "Actual Values", ylab = "Predicted Values")
grid(nx = NULL, ny = NULL, col = "grey80", lty = "dashed")
abline(0, 1, col = "red")  # Add a 1:1 line

# Combine actual and predicted values into a data frame
results <- data.frame(Actual = test_data$lSR, Predicted = predictions)

# create histograms with density distributions
# save to file
plot_filename <- paste(file_path, paste("actual_vs_predicted_distplot_subsetdata_pca.png", sep=""), sep="/")
png(plot_filename, width=1500, height=1400, res=250)
par(bg = "white")
ggplot(results, aes(x = Actual, fill = "Actual")) +
  geom_density(alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  labs(title = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
       x = "lSR", y = "Density") +
  scale_fill_manual(values = c("Actual" = "azure4", "Predicted" = "cadetblue1")) +
  guides(fill = guide_legend(title = ""))+
theme(
  panel.background = element_rect(fill = "white", color = "black"), 
  plot.background = element_rect(fill="white", color=NA),
  panel.grid.major = element_line(color = "grey80", linetype="dashed"),  # Darker major grid lines
#  panel.grid.minor = element_blank,  # Darker minor grid lines
  axis.line = element_line(color = "black"),  # Add axis lines
  axis.text = element_text(color = "black", size=11),  # Ensure axis text is black for visibility
  axis.title = element_text(color = "black", size=13),  # Ensure axis title is black for visibility
  legend.position = "bottom"
)
dev.off()
# print plot so screen
ggplot(results, aes(x = Actual, fill = "Actual")) +
  geom_density(alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  labs(title = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
       x = "lSR", y = "Density") +
  scale_fill_manual(values = c("Actual" = "darkmagenta", "Predicted" = "darkorange1")) +
  guides(fill = guide_legend(title = ""))
theme(
  panel.grid.major = element_line(color = "grey40"),  # Darker major grid lines
  panel.grid.minor = element_line(color = "grey60"),  # Darker minor grid lines
  axis.line = element_line(color = "black"),  # Add axis lines
  axis.text = element_text(color = "black"),  # Ensure axis text is black for visibility
  axis.title = element_text(color = "black")  # Ensure axis title is black for visibility
)


# compute model evaluation statistics
mse <- mean((predictions - test_data$lSR)^2)
rmse <- sqrt(mse)
mae <- mean(abs(predictions - test_data$lSR))
rsquared <- summary(lm_model)$r.squared
adjusted_r2 <- summary(lm_model)$adj.r.squared
model_summ <- summary(lm_model)
rse <- model_summ$sigma

# print model evaluation statistics
cat("Mean Squared Error (MSE):", mse, "\n")
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")
cat("Residual Standard Error (RSE):", rse, "\n")
cat("R-squared:", rsquared, "\n")
cat("Adjusted R-squared:", adjusted_r2, "\n")

# create dataframe with model eval results
model_eval_mlr_usa_lsr_df_best_null_pca <- data.frame(
  which_model = 'USA-Null',
  MSE = mse,
  RMSE = rmse,
  MAE = mae,
  RSE = rse,
  R_squared = rsquared,
  Adj_R_squared = adjusted_r2
)


### EVALUATE MODEL WITH BEST 4 PARAMETERS SUBSET DATA AND PCA
# split data into 70% training and 30% test
# make sure no nans in df
mlr_usa_lsr_df_best_4[is.na(mlr_usa_lsr_df_best_4)] <- 0
set.seed(121)  # for reproducibility
train_indices <- createDataPartition(mlr_usa_lsr_df_best_4$lSR, p = 0.7, list = FALSE)
train_data <- mlr_usa_lsr_df_best_4[train_indices, ]
test_data <- mlr_usa_lsr_df_best_4[-train_indices, ]

# train linear regression model
lm_model <- lm(formula_4, data = train_data)
print(summary(lm_model))

# predict on test data
predictions <- predict(lm_model, newdata = test_data)

# number of parameters included in lm model (minus the intercept)
n_params <- length(coef(lm_model))-1

# plot actual versus predicted values
# save to file
plot_filename <- paste(file_path, paste("actual_vs_predicted_1-1_subset_best4.png", sep=""), sep="/")
png(plot_filename,width=1400, height=1400, res=250)
# Remove the default grey background
par(bg = "white")
plot(test_data$lSR, predictions,
     main = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
     xlab = "Actual Values", ylab = "Predicted Values",
     cex.axis=0.95,
     cex.lab=1.1,
     col.lab='black',
     xlim = c(4, 18),   # Set x-axis limits
     ylim = c(4, 18))   # Set y-axis limits to be the same as x-axis
grid(nx = NULL, ny = NULL, col = "grey80", lty = "dashed")
abline(0, 1, col = "red")  # Add a 1:1 line
dev.off()
# print to screen
plot(test_data$lSR, predictions,
     main = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
     xlab = "Actual Values", ylab = "Predicted Values")
grid(nx = NULL, ny = NULL, col = "grey80", lty = "dashed")
abline(0, 1, col = "red")  # Add a 1:1 line

# Combine actual and predicted values into a data frame
results <- data.frame(Actual = test_data$lSR, Predicted = predictions)

# create histograms with density distributions
# save to file
plot_filename <- paste(file_path, paste("actual_vs_predicted_distplot_subset_best4.png", sep=""), sep="/")
png(plot_filename, width=1500, height=1400, res=250)
par(bg = "white")
ggplot(results, aes(x = Actual, fill = "Actual")) +
  geom_density(alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  labs(title = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
       x = "lSR", y = "Density") +
  scale_fill_manual(values = c("Actual" = "azure4", "Predicted" = "cadetblue1")) +
  guides(fill = guide_legend(title = ""))+
  theme(
    panel.background = element_rect(fill = "white", color = "black"), 
    plot.background = element_rect(fill="white", color=NA),
    panel.grid.major = element_line(color = "grey80", linetype="dashed"),  # Darker major grid lines
    #  panel.grid.minor = element_blank,  # Darker minor grid lines
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text = element_text(color = "black", size=11),  # Ensure axis text is black for visibility
    axis.title = element_text(color = "black", size=13),  # Ensure axis title is black for visibility
    legend.position = "bottom"
  )
dev.off()
# print plot so screen
ggplot(results, aes(x = Actual, fill = "Actual")) +
  geom_density(alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5, position = "identity", color = "black") + # Add color aesthetic for bin outlines
  labs(title = paste("Actual vs. Predicted Values (USA), lSR,", n_params, "Param(s)"),
       x = "lSR", y = "Density") +
  scale_fill_manual(values = c("Actual" = "darkmagenta", "Predicted" = "darkorange1")) +
  guides(fill = guide_legend(title = ""))
theme(
  panel.grid.major = element_line(color = "grey40"),  # Darker major grid lines
  panel.grid.minor = element_line(color = "grey60"),  # Darker minor grid lines
  axis.line = element_line(color = "black"),  # Add axis lines
  axis.text = element_text(color = "black"),  # Ensure axis text is black for visibility
  axis.title = element_text(color = "black")  # Ensure axis title is black for visibility
)


# compute model evaluation statistics
mse <- mean((predictions - test_data$lSR)^2)
rmse <- sqrt(mse)
mae <- mean(abs(predictions - test_data$lSR))
rsquared <- summary(lm_model)$r.squared
adjusted_r2 <- summary(lm_model)$adj.r.squared
model_summ <- summary(lm_model)
rse <- model_summ$sigma

# print model evaluation statistics
cat("Mean Squared Error (MSE):", mse, "\n")
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")
cat("Residual Standard Error (RSE):", rse, "\n")
cat("R-squared:", rsquared, "\n")
cat("Adjusted R-squared:", adjusted_r2, "\n")

# create dataframe with model eval results
model_eval_mlr_usa_lsr_df_best_4_pca <- data.frame(
  which_model = 'USA-Best4',
  MSE = mse,
  RMSE = rmse,
  MAE = mae,
  RSE = rse,
  R_squared = rsquared,
  Adj_R_squared = adjusted_r2
)


### EXPORT MODEL EVAL STATISTICS
model_eval_usa_merged <- rbind(model_eval_mlr_usa_lsr_df_all, model_eval_mlr_usa_lsr_df_best_null_pca, model_eval_mlr_usa_lsr_df_best_4_pca)
print(model_eval_usa_merged)
if (!dir.exists("model_evaluation/model_statistics")) {
  dir.create("model_evaluation/model_statistics", recursive = TRUE)
}
file_path <- paste0("model_evaluation/model_statistics/model_eval_stats_lsr_", version, ".csv")
write.csv(model_eval_usa_merged, file=file_path, row.names=FALSE)

# END MODEL TRAINING