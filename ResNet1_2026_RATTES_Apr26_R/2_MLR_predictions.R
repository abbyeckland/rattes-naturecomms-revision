
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

# Predictions - Version control -----------------------------------------

# update version of code below! 
# options are rattes_v1p2_silt, rattes_v1p2_sand, rattes_v1p2_clay, and rattes_v1p2_3param_silt
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
sa_input_file  <- settings$sa_input_file

# Print selected settings so user knows what is running
message("Running version: ", rattes_version)
message("Using SY input file: ", sy_input_file)
message("Using SA input file: ", sa_input_file)
message("Using predictors: ", paste(settings$predictors, collapse = ", "))

# Predictions - Data import and merging-----------------------------------------

### IMPORT SITE DATA
# set wd to the current folder
setwd(".")

# import study site data for mlr model building (we will use this dataframe later when making predictions)
ss_data <- read.csv(sy_input_file, header=TRUE, stringsAsFactors = FALSE)  # UPDATE WITH CURRENT VERSION
ss_data <- ss_data %>% rename("MAQ_NHDcms" = "MAQ_NHDcfs") # ensure label is correct for proper conversions (units are indeed cms)
####### make sure 'readin_lakecat_data.R' is updated with current model version and site data file before continuing

# then check if the lakecat CSV has already been processed and file exists in folder directory
run_readin_lakecat_data <- function() {
  file_path <- paste0("lakecat_ss_data_", pred_version, ".csv")
  
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
file_path <- paste0("lakecat_ss_data_", pred_version, ".csv")
lakecat_ss_data <- read.csv(file_path, header=TRUE, stringsAsFactors = FALSE)

# import ss_wbcomid, which binded study site sid and wbcomid
file_path <- paste0("ss_wbcomid_data_", pred_version, ".csv")
ss_wbcomid <- read.csv(file_path, header=TRUE, stringsAsFactors = FALSE)

lakecat_ss_data <- lakecat_ss_data %>% rename("MAQ_NHDcms" = "MAQ_NHDcfs") # ensure label is correct for proper conversions

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


### CLEAR WORKING MEMORY
# specify columns to keep
to_keep <- c("version", "rattes_version", "pred_version", "sy_input_file", "sa_input_file", "version_settings", "predictors", # version objects
             "ss_data", "ss_wbcomid", "lakecat_ss_fcpg_data", "df_ready", "constant_columns", "output_df_lsr",  # dataframes and their appropriate transformations
             "selected_predictors_USA_lsr_4", "selected_predictors_USA_lsr_null_pca",                                                     # predictor variables
             "intercept_usa_best4", "intercept_usa_null_pca",                                                                             # model intercepts
             "coefficients_usa_best4", "coefficients_usa_null_pca",                                                                       # model coefficients                                                                                      # bias correction term of residual variance
             "min_max_normalize", "perform_pca",
             "lambda_opt",
             "fcpg_subset",
             "filtered_site_data_nonans", "site_data_nonans",
             "site_data",
             "lower_bound_95", "upper_bound_95", "t_critical_95",
             "coef_beta", "vcov_beta","rse","dof","t_critical_95",  # NEW CONFIDENCE INTERVAL TESTING CODE
             "lock_sids")
# all objects in environment
all_objects <- ls()
# specify objects to remove
objects_to_remove <- setdiff(all_objects, to_keep)
# remove objects
rm(list= objects_to_remove)
# check remaining objects
ls()


### PREPARE DATAFRAMES
# subset study pred data for model building
pred_data <- lakecat_ss_fcpg_data[lakecat_ss_fcpg_data$issite==0 & 
                                  (lakecat_ss_fcpg_data$islock==0 | (lakecat_ss_fcpg_data$islock==1 & lakecat_ss_fcpg_data$PermStorag == 1)),
                                  ]

# with locks removed, now see how many sites have lakecat data
pred_data_nan_check <- pred_data %>% filter(!is.na(WsAreaSqKm))
print(nrow(pred_data_nan_check))
print(nrow(pred_data))



# Predictions - Data pre-processing ---------------------------------------

### DEAL WITH MISSING OR BAD DATA
# find the starting index based on column name
start_index <- which(colnames(pred_data) == "PctAg2006Slp20Ws") # beginning index of lakecat column data (first parameter in list should be "PctAg2006Slp20Ws")
# create list of columns from that point onward
param_list_lakecat <- colnames(pred_data)[start_index:ncol(pred_data)]
print(param_list_lakecat)

# not adding elev_m since all nan values. check in future input files though if this is fixed
additional_cols <- c('damlength_m', 'SA_m2', 'damH', 'NrX_Final', 'NrY_Final', 'DA', 'MAQ_NHDcms', 'D50', 'trapp', 'pathx', 'yrc')
param_list <- c(additional_cols, param_list_lakecat)

# print the number of nan values in each column (if no nans in the ss_data columns, this is great)
pred_data_param_list <- pred_data[, param_list]
sapply(pred_data_param_list, function(x) sum(is.na(x)))

# drop rows from dataframe that contain all nans but only for columns within param_list_lakecat - SKIP FOR NOW, CAN STILL MAKE PREDICTIONS WITHOUT THIS DATA
#pred_data_nonans <- pred_data[!apply(pred_data[param_list_lakecat], 1, function(row) all(is.na(row))), ]
pred_data_nonans <- pred_data
nrow(pred_data_nonans)

### OUTLIER REMOVAL
# identify outliers prior to transformations
no_outlier_removal_from <- c("trapp", "pathx", "yrc", "yrp", "DA", "NrX_Final", "NrY_Final")
param_list_outlier_check <- param_list[!param_list %in% no_outlier_removal_from]

# loop over each variable in param_list
sink("info/parameter_outlier_removal_prediction_sites.txt", split = TRUE)
for (variable in param_list_outlier_check) {
  # exclude NA and zero values
  valid_data <- pred_data_nonans[[variable]][!is.na(pred_data_nonans[[variable]]) & pred_data_nonans[[variable]] != 0]
  # standardize the variable (only valid data)
  z_scores <- scale(valid_data)
  # identify outliers based on z-scores
  outliers <- which(abs(z_scores) > 5)
  # print the number of outliers for the current variable
  cat("Number of outliers in", variable, ":", length(outliers), "\n")
  # map the outlier indices back to the original dataset
  all_indices <- which(!is.na(pred_data_nonans[[variable]]) & pred_data_nonans[[variable]] != 0)
  pred_data_nonans[all_indices[outliers], variable] <- NA
}
sink()

# check if any columns have negative values (avoids any later transformation issues)
negative_columns <- names(pred_data_nonans)[sapply(pred_data_nonans, function(x) any(x < 0 & !is.na(x)))]
# print the column names with negative values
print(negative_columns)
# there are some other parameter columns with negative values (like slope) -- will deal with this automatically next
negative_columns <- negative_columns[negative_columns %in% param_list]
print(negative_columns)

# apply min-max normalization to column tmin and save as new column
pred_data_nonans$Tmin8110Ws_norm <- min_max_normalize(pred_data_nonans$Tmin8110Ws)
pred_data_nonans$Tmean8110Ws_norm <- min_max_normalize(pred_data_nonans$Tmean8110Ws)
pred_data_nonans$Tmax8110Ws_norm <- min_max_normalize(pred_data_nonans$Tmax8110Ws)
pred_data_nonans$tmin_norm <- min_max_normalize(pred_data_nonans$tmin)
pred_data_nonans$tmax_norm <- min_max_normalize(pred_data_nonans$tmax)

# reverse the sign of NrX_Final
pred_data_nonans$NrX_Final <- abs(pred_data_nonans$NrX_Final)

# remove negative elev values
pred_data_nonans$elev[pred_data_nonans$elev < 0] <- NA

# check if any columns have a placeholder for missing data. Missing data is represented by -9998
count_missing <- colSums(pred_data_nonans == -9998 | pred_data_nonans == -9999, na.rm=TRUE)
cols_with_missing <- names(count_missing[count_missing > 0])
print(count_missing)
print(cols_with_missing)

# filter param_list to include only columns that are in cols_with_missing
cols_to_check <- param_list[param_list %in% cols_with_missing]
# remove rows with missing data
pred_data_nonans[cols_to_check] <- lapply(pred_data_nonans[cols_to_check], function(col) {
  col[col == -9998 | col == -9998] <- NA
  return(col)
})

# remove any parameters that have constant values (std is 0) - get from study site data pre-processing
# Drop the identified columns from the dataframe
pred_data_nonans <- pred_data_nonans[, !(names(pred_data_nonans) %in% constant_columns)]
# update parameter list after removing constant columns
param_list <- param_list[!(param_list %in% constant_columns)]
param_list_lakecat <- param_list_lakecat[!(param_list_lakecat %in% constant_columns)]


# Calculate the percentage of missing data needed to run the model in each column
missing_percent <- colMeans(is.na(pred_data_nonans)) * 100
cols_with_missing <- names(missing_percent[missing_percent > 0])
ss_cols <- c("trapp", "DA", "yrp", "capp", "kappa")
cols_with_missing <- cols_with_missing[cols_with_missing %in% ss_cols]

# Check if there are any required columns with missing data
if (length(cols_with_missing) >= 1) {
  # Print rows with missing data only in the 'ss_cols' columns
  missing_rows <- pred_data_nonans[apply(pred_data_nonans[ss_cols], 1, function(row) any(is.na(row))), ]
  print(missing_rows)
  stop("Warning: required columns have NAs.")
}

# remove rows with missing data only within cols_with_missing
pred_data_nonans <- pred_data_nonans[complete.cases(pred_data_nonans[cols_with_missing]), ]
sink("info/number_of_prediction_sites.txt", split = TRUE)
print(nrow(pred_data_nonans))
sink()



# Predictions - Data transformations -------------------------------------------

### DATA TRANSFORMATIONS
# create/combine parameters from existing parameters
# combine all forest types into one forest param
pred_data_nonans$pct_frst_all_2006 <- rowSums(pred_data_nonans[,c("PctDecid2006Ws","PctConif2006Ws","PctMxFst2006Ws")], na.rm=TRUE)
pred_data_nonans$pct_frst_all_2006[is.nan(pred_data_nonans$PctDecid2006Ws) & is.nan(pred_data_nonans$PctConif2006Ws) & is.nan(pred_data_nonans$PctMxFst2006Ws)] <- NaN
pred_data_nonans$pct_frst_all_2011 <- rowSums(pred_data_nonans[,c("PctDecid2011Ws","PctConif2011Ws","PctMxFst2011Ws")], na.rm=TRUE)
pred_data_nonans$pct_frst_all_2011[is.nan(pred_data_nonans$PctDecid2011Ws) & is.nan(pred_data_nonans$PctConif2011Ws) & is.nan(pred_data_nonans$PctMxFst2011Ws)] <- NaN

# combine all wetland types into one wetland param
pred_data_nonans$pct_wetl_all_2006 <- rowSums(pred_data_nonans[,c("PctWdWet2006Ws","PctHbWet2006Ws")], na.rm=TRUE)
pred_data_nonans$pct_wetl_all_2006[is.nan(pred_data_nonans$PctWdWet2006Ws) & is.nan(pred_data_nonans$PctHbWet2006Ws)] <- NaN
pred_data_nonans$pct_wetl_all_2011 <- rowSums(pred_data_nonans[,c("PctWdWet2011Ws","PctHbWet2011Ws")], na.rm=TRUE)
pred_data_nonans$pct_wetl_all_2011[is.nan(pred_data_nonans$PctWdWet2011Ws) & is.nan(pred_data_nonans$PctHbWet2011Ws)] <- NaN

# combine percent glacial till, glacial lake sediments, and percent ice into a single param
pred_data_nonans$pct_icy_sed <- rowSums(pred_data_nonans[,c("PctGlacTilClayWs","PctGlacTilLoamWs","PctGlacTilCrsWs","PctGlacLakeCrsWs","PctGlacLakeFineWs")], na.rm=TRUE)
pred_data_nonans$pct_icy_sed[is.nan(pred_data_nonans$PctGlacTilClayWs) & is.nan(pred_data_nonans$PctGlacTilLoamWs) & is.nan(pred_data_nonans$PctGlacTilCrsWs) & is.nan(pred_data_nonans$PctGlacLakeCrsWs) & is.nan(pred_data_nonans$PctGlacLakeFineWs)] <- NaN

# combine all eolian lithology types into one eolian param
pred_data_nonans$pct_eolian <- rowSums(pred_data_nonans[,c("PctEolCrsWs","PctEolFineWs")], na.rm=TRUE)
pred_data_nonans$pct_eolian[is.nan(pred_data_nonans$PctEolCrsWs) & is.nan(pred_data_nonans$PctEolFineWs)] <- NaN

# combine all agslope lithology types into one param
pred_data_nonans$pct_ag_on_slp <- rowSums(pred_data_nonans[,c("PctAg2006Slp20Ws","PctAg2006Slp10Ws")], na.rm=TRUE)
pred_data_nonans$pct_ag_on_slp[is.nan(pred_data_nonans$PctAg2006Slp20Ws) & is.nan(pred_data_nonans$PctAg2006Slp10Ws)] <- NaN

# combine annual forest fire area and forest loss into a single decadal parameter
pred_data_nonans$pct_fire_0010Ws <- rowSums(pred_data_nonans[,c("PctFire2000Ws","PctFire2001Ws","PctFire2002Ws","PctFire2003Ws","PctFire2004Ws","PctFire2005Ws","PctFire2006Ws","PctFire2007Ws","PctFire2008Ws","PctFire2009Ws","PctFire2010Ws")], na.rm=TRUE)
pred_data_nonans$pct_fire_0010Ws[is.nan(pred_data_nonans$PctFire2000Ws) & is.nan(pred_data_nonans$PctFire2001Ws) & is.nan(pred_data_nonans$PctFire2002Ws) & is.nan(pred_data_nonans$PctFire2003Ws) & is.nan(pred_data_nonans$PctFire2004Ws) & is.nan(pred_data_nonans$PctFire2005Ws) & is.nan(pred_data_nonans$PctFire2006Ws) & is.nan(pred_data_nonans$PctFire2007Ws) & is.nan(pred_data_nonans$PctFire2008Ws) & is.nan(pred_data_nonans$PctFire2009Ws) & is.nan(pred_data_nonans$PctFire2010Ws)] <- NaN
pred_data_nonans$pct_frstloss_0113Ws <- rowSums(pred_data_nonans[,c("PctFrstLoss2001Ws","PctFrstLoss2002Ws","PctFrstLoss2003Ws","PctFrstLoss2004Ws","PctFrstLoss2005Ws","PctFrstLoss2006Ws","PctFrstLoss2007Ws","PctFrstLoss2008Ws","PctFrstLoss2009Ws","PctFrstLoss2010Ws","PctFrstLoss2011Ws","PctFrstLoss2012Ws","PctFrstLoss2013Ws")], na.rm=TRUE)
pred_data_nonans$pct_frstloss_0113Ws[is.nan(pred_data_nonans$PctFrstLoss2001Ws) & is.nan(pred_data_nonans$PctFrstLoss2002Ws) & is.nan(pred_data_nonans$PctFrstLoss2003Ws) & is.nan(pred_data_nonans$PctFrstLoss2004Ws) & is.nan(pred_data_nonans$PctFrstLoss2005Ws) & is.nan(pred_data_nonans$PctFrstLoss2006Ws) & is.nan(pred_data_nonans$PctFrstLoss2007Ws) & is.nan(pred_data_nonans$PctFrstLoss2008Ws) & is.nan(pred_data_nonans$PctFrstLoss2009Ws) & is.nan(pred_data_nonans$PctFrstLoss2010Ws) & is.nan(pred_data_nonans$PctFrstLoss2011Ws)& is.nan(pred_data_nonans$PctFrstLoss2012Ws) & is.nan(pred_data_nonans$PctFrstLoss2013Ws)] <- NaN

# Calculate temp difference for Tmax8110Ws and Tmin8110Ws, only if both are not NaN and not zero
pred_data_nonans$tdiff_8110Ws <- ifelse(
  !is.na(pred_data_nonans$Tmax8110Ws) & !is.na(pred_data_nonans$Tmin8110Ws) & 
    pred_data_nonans$Tmax8110Ws != 0 & pred_data_nonans$Tmin8110Ws != 0,
  pred_data_nonans$Tmax8110Ws - pred_data_nonans$Tmin8110Ws,
  NA
)

# Calculate temp difference for tmax and tmin, only if both are not NaN and not zero
pred_data_nonans$tdiff_fcpg <- ifelse(
  !is.na(pred_data_nonans$tmax) & !is.na(pred_data_nonans$tmin) & 
    pred_data_nonans$tmax != 0 & pred_data_nonans$tmin != 0,
  pred_data_nonans$tmax - pred_data_nonans$tmin,
  NA
)

# Calculate specific discharge, flow to drainage area ratio, only if both are not NaN and not zero
pred_data_nonans$fdar <- ifelse(
  !is.na(pred_data_nonans$DA) & !is.na(pred_data_nonans$MAQ_NHDcms) & 
    pred_data_nonans$DA > 0 & pred_data_nonans$MAQ_NHDcms > 0,
  (pred_data_nonans$MAQ_NHDcms / pred_data_nonans$DA)*100,
  NA
)

# Calculate Relief as the difference between ElevWs and ElevCat, only if both are not NaN and not zero
pred_data_nonans$relief <- ifelse(
  !is.na(pred_data_nonans$ElevWs) & !is.na(pred_data_nonans$ElevCat) & 
    pred_data_nonans$ElevWs > 0 & pred_data_nonans$ElevCat > 0 & pred_data_nonans$ElevWs > pred_data_nonans$ElevCat,
  pred_data_nonans$ElevWs - pred_data_nonans$ElevCat,
  NA
)

# check how many nans in each param
sapply(pred_data_nonans, function(x) sum(is.na(x)))


# now replace any other 0 value with small number (transformations don't like if value of 0)
pred_data_nonans[pred_data_nonans == 0] <- 1e-6  # grab any needed extra columns from this df

# make new dataframe with only desired columns
new_cols <- c("pct_frst_all_2006", "pct_frst_all_2011", "pct_wetl_all_2006", "pct_wetl_all_2011", "pct_icy_sed", "pct_eolian", "pct_ag_on_slp", "pct_fire_0010Ws", "pct_frstloss_0113Ws", "tdiff_8110Ws", "tmin_norm", "tmax_norm", "tdiff_fcpg", "relief", "fdar")
all_cols <- c(param_list, new_cols)
# ready datafarme
pred_df_ready <- pred_data_nonans[, all_cols]
# define some columns to drop for data we don't need
cols_to_drop_major <- c("COMID_y", "COMID", # don't need these now
                        "DA",  
                        "Tmin8110Ws", "Tmean8110Ws", "Tmax8110Ws",# using "Tmin8110Ws_norm" et al instead since no negative numbers
                        "tmin", "tmax", #using norms
                        "PctAg2006Slp20Ws", "PctAg2006Slp10Ws",
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
                        "DamDensWs", "DamNIDStorWs", "DamNrmStorWs", # actually, don't want anyting else about upstream damming since we are including this separately
                        "CoalMineDensWs", # there is another param called "MineDensWs" that has a better corr with lSDR
                        "PctOw2006Ws", "PctOw2011Ws", # there is another param called "PctWaterWs" that has a better corr with lSDR
                        "REGION", "LRR_SYMBOL", # from region join. will add this later
                        "Precip8110Ws", "ElevWs", "ElevCat") # already have complete datasets from fcpg
pred_df_ready <- pred_df_ready[, !names(pred_df_ready) %in% cols_to_drop_major]

# (optional) save pred_df_ready prior to transformations to file
file_path <- paste0("prediction_data_pre-transforms_", pred_version, ".csv")
write.csv(pred_df_ready, file=file_path, row.names=FALSE)

# remove columns of box-cox transformations that we can't apply right now
# don't transform these columns:
dont_transform <- c("wSAatdam", "AveTrap", "yrs_tot", "wrt", "wSAratio") # these require updating each year in the timeloop, and will be transformed within the loop
pred_output_df_lsr <- output_df_lsr[!output_df_lsr$variable %in% dont_transform, ]


### APPLY TRANSFORMATIONS TO DATA
# Initialize dataframe to store transformed data
pred_transformed_df_lsr <- pred_df_ready
# Iterate over each row in output_df
for (i in 1:nrow(pred_output_df_lsr)) {
  # Get the variable name and its corresponding lambda value
  variable_name <- pred_output_df_lsr$variable[i]
  lambda1_value <- pred_output_df_lsr$lambda1[i]
  lambda2_value <- pred_output_df_lsr$lambda2[i]
  # Filter out nan values for the current variable
  non_na_values <- na.omit(pred_df_ready[[variable_name]])
  # Filter out zero values -                                                    COMMENT OUT TO KEEP ZEROS
  non_na_values <- non_na_values[non_na_values > 1e-6]

  # Apply the corresponding transformation using the lambda value to non-nan values
  if (lambda1_value == 0) {
    transformed_values <- log(non_na_values + lambda2_value)
  } else if (lambda1_value == 1) {
    transformed_values <- non_na_values  # linear transformation
  } else {
    transformed_values <- ((non_na_values + lambda2_value)^lambda1_value - 1)/lambda1_value
  }

  # Create a vector with the same length as the original variable, filled with NA
  transformed_vector <- rep(NA, length(pred_df_ready[[variable_name]]))
  # Replace the NA positions with the transformed values                        REMOVE '& (pred_df_ready[[col_name]] > 1e-6)' TO KEEP ZEROS
  valid_indices <- which(!is.na(pred_df_ready[[variable_name]]) & pred_df_ready[[variable_name]] > 1e-6)
  transformed_vector[valid_indices] <- transformed_values

  # Update the dataframe with the transformed values
  pred_transformed_df_lsr[[variable_name]] <- transformed_vector
}

# View the transformed dataframe
print(pred_transformed_df_lsr)



# drop LSR column
pred_transformed_df_lsr_ <- pred_transformed_df_lsr
# save transformed data to folder
file_path <- paste0("transformations/prediction_data_transformed_", pred_version, ".csv")
write.csv(pred_transformed_df_lsr_, file=file_path, row.names=FALSE)




# Predictions - PCA ------------------------------------------------------------

# Create PCs for LSR
# first rename some columns
columns_to_rename <- c("Pestic97Ws" = "Pestic1997Ws", "pct_fire_0010Ws" = "pct_fire_decade_2010Ws", "pct_frstloss_0113Ws" = "pct_frstloss_decade_2013Ws")

# Rename the columns based on the mapping
names(pred_transformed_df_lsr_)[names(pred_transformed_df_lsr_) %in% names(columns_to_rename)] <- columns_to_rename[names(pred_transformed_df_lsr_)[names(pred_transformed_df_lsr_) %in% names(columns_to_rename)]]

# Apply the function to the dataframe
pred_transformed_df_lsr_fill <- pred_transformed_df_lsr_

# soil characteristics (pc_soil)
pca_cols_lsr_soil <- c("ClayWs", "SandWs", "PermWs", "OmWs")
data_with_pca_lsr_soil <- perform_pca(pred_transformed_df_lsr_fill, pca_cols_lsr_soil, "pc_soil")
variance_pc_soil <- attr(data_with_pca_lsr_soil, "variance_explained")
print(variance_pc_soil)

# dam dimensions (pc_damdim)
pca_cols_lsr_damdim <- c("damlength_m", "damH")
data_with_pca_lsr_damdim <- perform_pca(pred_transformed_df_lsr_fill, pca_cols_lsr_damdim, "pc_damdim")
variance_pc_damdim <- attr(data_with_pca_lsr_damdim, "variance_explained")
print(variance_pc_damdim)

# temperature (pc_temp_fcpg)
pca_cols_lsr_temp_fcpg <- c("tmin_norm", "tmax_norm", "NrY_Final")
data_with_pca_lsr_temp_fcpg <- perform_pca(pred_transformed_df_lsr_fill, pca_cols_lsr_temp_fcpg, "pc_temp_fcpg")
variance_pc_temp_fcpg <- attr(data_with_pca_lsr_temp_fcpg, "variance_explained")
print(variance_pc_temp_fcpg)

# general agriculture (pc_ag)
pca_cols_lsr_ag <- c("CBNFWs", "AgKffactWs", "Pestic1997Ws", "FertWs", "PctCrop2006Ws")
data_with_pca_lsr_ag <- perform_pca(pred_transformed_df_lsr_fill, pca_cols_lsr_ag, "pc_ag")
variance_pc_ag <- attr(data_with_pca_lsr_ag, "variance_explained")
print(variance_pc_ag)

# pollution (pc_poll)
pca_cols_lsr_poll <- c("SuperfundDensWs", "TRIDensWs", "MineDensWs", "NPDESDensWs")
data_with_pca_lsr_poll <- perform_pca(pred_transformed_df_lsr_fill, pca_cols_lsr_poll, "pc_poll")
variance_pc_poll <- attr(data_with_pca_lsr_poll, "variance_explained")
print(variance_pc_poll)

# urbanization (pc_urban_2006)
pca_cols_lsr_urban_2006 <- c("PctImp2006Ws", "PctUrbOp2006Ws", "PctUrbLo2006Ws", "PctUrbMd2006Ws", "RdDensWs", "RdCrsWs", "HUDen2010Ws", "PopDen2010Ws")
data_with_pca_lsr_urban_2006 <- perform_pca(pred_transformed_df_lsr_fill, pca_cols_lsr_urban_2006, "pc_urban_2006")
variance_pc_urban_2006 <- attr(data_with_pca_lsr_urban_2006, "variance_explained")
print(variance_pc_urban_2006)

# carbonate residuals / non carbonate residuals (pc_carb)
pca_cols_lsr_carb <- c("PctCarbResidWs", "PctNonCarbResidWs")
data_with_pca_lsr_carb <- perform_pca(pred_transformed_df_lsr_fill, pca_cols_lsr_carb, "pc_carb")
variance_pc_carb <- attr(data_with_pca_lsr_carb, "variance_explained")
print(variance_pc_carb)

# water (pc_water)
pca_cols_lsr_water <- c("SA_m2", "MAQ_NHDcms")
data_with_pca_lsr_water <- perform_pca(pred_transformed_df_lsr_fill, pca_cols_lsr_water, "pc_water")
variance_pc_water <- attr(data_with_pca_lsr_water, "variance_explained")
print(variance_pc_water)

# merge pc data into pred_transformed_df_lsr_fill df
pred_transformed_df_lsr_pca <- cbind(pred_transformed_df_lsr_, 
                                     data_with_pca_lsr_soil, 
                                     data_with_pca_lsr_damdim, 
                                     data_with_pca_lsr_temp_fcpg, 
                                     data_with_pca_lsr_ag, 
                                     data_with_pca_lsr_poll, 
                                     data_with_pca_lsr_urban_2006,
                                     data_with_pca_lsr_carb,
                                     data_with_pca_lsr_water)



# Predictions - Interaction terms ----------------------------------------------


# create copy of pca and transformed dataframe
pred_mlr_usa_lsr_df <- pred_transformed_df_lsr_pca  # make copy so don't have to rerun above


# Only create interaction term for params with missing data. Temporal params will be considered in the MLR annual prediction loop

# initialize a vector to store column names with NaN values
columns_with_missing <- c()

# create list of columns to check for nans, na or values >= 0 and <= 1e-6
pred_df_ready_cols <- colnames(pred_mlr_usa_lsr_df)


# iterate over each column specified in the list
for (col in pred_df_ready_cols) {
  # check if the column contains NaN, Na, or values >= 0 and < 1e-6
  if (any(is.nan(pred_mlr_usa_lsr_df[[col]])) ||
      any(is.na(pred_mlr_usa_lsr_df[[col]]))) {
    # if NaN values are found, store the column name
    columns_with_missing <- c(columns_with_missing, col)
  }
}

# print column names with NaN values
print(columns_with_missing)

# create interaction terms for columns with missing and zero-inflated data
for (col in columns_with_missing) {
  # create indicator variable for current lc column with missing data
  pred_mlr_usa_lsr_df[[paste0("indicator_", col)]] <- as.integer(!is.na(pred_mlr_usa_lsr_df[[col]]))
}
print(colnames(pred_mlr_usa_lsr_df))

# convert nans in pred_mlr_usa_lsr_df to a number that can be multiplied by the indicator function
pred_mlr_usa_lsr_df[is.na(pred_mlr_usa_lsr_df)] <- 0




# Predictions - Model setup ----------------------------------------------------

### PREPARE DF FOR PREDICTION LOOP

# make a copy so can rerun from here
pred_mlr_usa_lsr_df_working <- pred_mlr_usa_lsr_df

# add non-transformed columns back to df
pred_mlr_usa_lsr_df_working$SID <- pred_data_nonans$SID
pred_mlr_usa_lsr_df_working$DA_mi2 <- pred_data_nonans$DA /2.59                 # convert km2 to mi2
pred_mlr_usa_lsr_df_working$DA_km2 <- pred_data_nonans$DA                       # save DA in km2
pred_mlr_usa_lsr_df_working$k <- pred_data_nonans$kappa                         # get kappa term for Brown equation
pred_mlr_usa_lsr_df_working$yrp <- pred_data_nonans$yrp                         # get year to begin predictions
pred_mlr_usa_lsr_df_working$capp <- pred_data_nonans$capp                       # get capacity at yrp
pred_mlr_usa_lsr_df_working$trapp <- pred_data_nonans$trapp                     # get trap efficiency at yrp
pred_mlr_usa_lsr_df_working$yrr <- pred_data_nonans$yrr                         # get year dam removed
pred_mlr_usa_lsr_df_working$sitetags_1 <- pred_data_nonans$sitetags_1           # get site tag, indicating if dam is inside or outside study basins
pred_mlr_usa_lsr_df_working$MAQ_NHDcms <- pred_data_nonans$MAQ_NHDcms           # get MAQ to create wrt term
pred_mlr_usa_lsr_df_working$yrc <- pred_data_nonans$yrc                         # get yrc in case differs from yrp
# Replace any value of 1e-6 with NA in the MAQ_NHDcms column
pred_mlr_usa_lsr_df_working$MAQ_NHDcms[pred_mlr_usa_lsr_df_working$MAQ_NHDcms == 1e-6] <- NA

# Sort dataframe by SID
pred_mlr_usa_lsr_df_working <- pred_mlr_usa_lsr_df_working[order(pred_mlr_usa_lsr_df_working$SID), ]


# Save intercept term from selected model
intercept <- as.numeric(intercept_usa_best4)   # intercept term in model in list
print(intercept)

# Save coefficients from selected model
coefficients <- coefficients_usa_best4
print(coefficients)

# Update coefficient names
names(coefficients) <- sub("wrt:indicator_wrt", "wrt", names(coefficients))
print(names(coefficients))


# Prep input matrices for prediction loop 

### STOP
# make sure 'compute_wSAatdam.R' is updated with current model version, site data file, and wSAatdam file
# then check if the wSAatdam CSV has already been processed and exists in folder directory
run_compute_wSAatdam <- function() {
  file_path <- paste0("wSAatdam_predSites_notTransformed_", pred_version, ".csv")
  
  # check if file exists already
  if (!file.exists(file_path)) {
    # file does not exist, so run the script
    message("CSV file is not found. Running 'compute_wSAatdam.R'...")
    source("compute_wSAatdam.R")
  } else {
    message("CSV file already exists. Skipping 'compute_wSAatdam.R'.")
  }
}

# call the function
run_compute_wSAatdam()

### PREP WSAATDAM
# Import wSAatdam csv
wSAatdam_df <- read.csv(paste0("wSAatdam_predSites_notTransformed_", pred_version, ".csv"), header=TRUE, stringsAsFactors = FALSE)
colnames(wSAatdam_df) <- sub("^X", "", colnames(wSAatdam_df))
print(colnames(wSAatdam_df))
nrow(wSAatdam_df)

# Match sids from wSAatdam df to current df
wSAatdam_df <- wSAatdam_df[wSAatdam_df$SID %in% pred_mlr_usa_lsr_df_working$SID, ]
pred_mlr_usa_lsr_df_working <- pred_mlr_usa_lsr_df_working[pred_mlr_usa_lsr_df_working$SID %in% wSAatdam_df$SID, ]

# make sure wSAatdam_df SIDs ordered
wSAatdam_df <- wSAatdam_df[order(wSAatdam_df$SID), ]

# check to ensure SID lists and order are the same
differences <- pred_mlr_usa_lsr_df_working$SID != wSAatdam_df$SID
if (any(differences)) {
  print("Differences found:")
  print(data.frame(Index = which(differences), pred_mlr_usa_lsr_df_working = pred_mlr_usa_lsr_df_working$SID[differences], wSAatdam_df = wSAatdam_df$SID[differences]))
} else {
  print("No differences found.")
}

# SUBSET FOR TESTING (ONLY RUNS 1000 PREDICTION SITES NOT ALL. COMMENT OUT NEXT TWO LINES TO RUN FULL CODE)

# flags <- c(285187, 294252, 294642, 294700, 296998, 297420, 312694, 332735, 332772, 332845, 332957, 333054, 333164, 334019, 343932, 351306, 351697, 352214, 352693, 352697, 352708, 352712, 352789, 356839, 359236)
# wSAatdam_df <- wSAatdam_df[wSAatdam_df$SID %in% flags, ]
# pred_mlr_usa_lsr_df_working <- pred_mlr_usa_lsr_df_working[pred_mlr_usa_lsr_df_working$SID %in% flags, ]

# # NEW FOR CONFIDENCE INTERVAL TESTING
# wSAatdam_df <- wSAatdam_df[1:1000,]
# pred_mlr_usa_lsr_df_working <- pred_mlr_usa_lsr_df_working[1:1000,]
# # END CONFIDENCE INTERVAL TESTING


# create variables to set up matrices
SID <- wSAatdam_df$SID
yrp <- wSAatdam_df$yrp
trapp <- (pred_mlr_usa_lsr_df_working$trapp*100)
capp <- pred_mlr_usa_lsr_df_working$capp
years <- 1699:2050
col_names <- as.character(years)

# convert wSAatdam_df into a matrix
wSAatdam_m <- as.matrix(wSAatdam_df[, -c(1,2)])
colnames(wSAatdam_m) <- col_names
rownames(wSAatdam_m) <- SID

# transform wSAatdam data based on identified best transformation
# Step 1: Extract lambda value for 'wSAatdam' from output_df_lsr
lambda1_wSAatdam <- output_df_lsr[output_df_lsr$variable == 'wSAatdam', 'lambda1']
lambda2_wSAatdam <- output_df_lsr[output_df_lsr$variable == 'wSAatdam', 'lambda2']

transformation_wSAatdam <- if (lambda1_wSAatdam == 0) {
  function(x) log(x + lambda2_wSAatdam)
} else if (lambda1_wSAatdam == 1) {
  function(x) x
} else {
  function(x) ((x + lambda2_wSAatdam)^lambda1_wSAatdam - 1) / lambda1_wSAatdam
}

# Step 3: Apply transformation to wSAatdam_df
wSAatdam_m_transformed <- wSAatdam_m

# transform all values in matrix
wSAatdam_m_transformed <- sapply(wSAatdam_m, transformation_wSAatdam)

# Convert the result back to a matrix with the original dimensions
wSAatdam_m_transformed <- matrix(wSAatdam_m_transformed, nrow = nrow(wSAatdam_m), ncol = ncol(wSAatdam_m))

# Reassign the column and row names to the transformed matrix
colnames(wSAatdam_m_transformed) <- colnames(wSAatdam_m)
rownames(wSAatdam_m_transformed) <- rownames(wSAatdam_m)

# Print the transformed matrix to verify
print(wSAatdam_m_transformed)



### PREP TE MATRIX
trap_m <- matrix(NA, nrow = length(SID), ncol = length(col_names))
colnames(trap_m) <- col_names
rownames(trap_m) <- SID

for (i in seq_along(SID)) {
  yrp_value <- yrp[i]

  # replace NA with 0 for years < yrp_value
  trap_m[i, col_names[years < yrp_value]] <- 0

  # assign trapp value to the column corresponding to the yrp_value
  trap_m[i, as.character(yrp_value)] <- trapp[i]
}


# Make a copy to calculate average transformed trap
avetrap_m <- trap_m


# transform trap_df data based on identified best transformation
# Step 1: Extract lambda value for 'AveTrap' from output_df_lsr
lambda1_avetrap <- output_df_lsr[output_df_lsr$variable == 'AveTrap', 'lambda1']           #'AveTrap' since this is the variable model built on and saved in output_df_lsr
lambda2_avetrap <- output_df_lsr[output_df_lsr$variable == 'AveTrap', 'lambda2']           #'AveTrap' since this is the variable model built on and saved in output_df_lsr

# set up trap efficiency matrix
trap_m_transformed <- matrix(NA, nrow = length(SID), ncol = length(col_names))
colnames(trap_m_transformed) <- col_names
rownames(trap_m_transformed) <- SID

transformation_avetrap <- if (lambda1_avetrap == 0) {
  function(x) log(x + lambda2_avetrap)
} else if (lambda1_avetrap == 1) {
  function(x) x
} else {
  function(x) ((x + lambda2_avetrap)^lambda1_avetrap - 1) / lambda1_avetrap
}

# Step 3: Apply transformation to trap_m
for (i in seq_along(SID)) {
  yrp_value <- yrp[i]
  # replace 0 cap with capp for yrp value
  trap_m_transformed[i, as.character(yrp_value)] <- transformation_avetrap(trapp[i])
}

# make a copy so we can update all dfs and make sure working properly
avetrap_m_transformed <- trap_m_transformed


### PREP CAP AND REMAINING MATRICES
cap_m <- matrix(0, nrow = length(SID), ncol = length(col_names))
colnames(cap_m) <- col_names
rownames(cap_m) <- SID

for (i in seq_along(SID)) {
  yrp_value <- yrp[i]
  # replace 0 cap with capp for yrp value
  cap_m[i, as.character(yrp_value)] <- capp[i]
}

# initialize matrix to keep track of cap by CI limits
cap_up_95_m <- cap_m
cap_low_95_m <- cap_m

# set up lsr matrix
lsr_m <- matrix(NA, nrow = length(SID), ncol = length(col_names))
colnames(lsr_m) <- col_names
rownames(lsr_m) <- SID

# copy lsr m to inital sr, avesr
sr_m <- lsr_m
avesr_m <- lsr_m

sr_up_95_m <- lsr_m
sr_low_95_m <- lsr_m

# set up sed vol matrix
sv_m <- matrix(0, nrow = length(SID), ncol = length(col_names))
colnames(sv_m) <- col_names
rownames(sv_m) <- SID

# initalize matrix to keep track of cap by CI limits
sv_up_95_m <- sv_m
sv_low_95_m <- sv_m


# prep transformations for yrs_tot
# Step 1: Extract lambda value for 'yrs_tot' from output_df_lsr
lambda1_yrstot <- output_df_lsr[output_df_lsr$variable == 'yrs_tot', 'lambda1']
lambda2_yrstot <- output_df_lsr[output_df_lsr$variable == 'yrs_tot', 'lambda2']

transformation_yrstot <- if (lambda1_yrstot == 0) {
  function(x) log(x + lambda2_yrstot)
} else if (lambda1_yrstot == 1) {
  function(x) x
} else {
  function(x) ((x + lambda2_yrstot)^lambda1_yrstot - 1) / lambda1_yrstot
}

# prep transformations for wrt
# Step 1: Extract lambda value for 'wrt' from output_df_lsr
lambda1_wrt <- output_df_lsr[output_df_lsr$variable == 'wrt', 'lambda1']
lambda2_wrt <- output_df_lsr[output_df_lsr$variable == 'wrt', 'lambda2']

transformation_wrt <- if (lambda1_wrt == 0) {
  function(x) log(x + lambda2_wrt)
} else if (lambda1_wrt == 1) {
  function(x) x
} else {
  function(x) ((x + lambda2_wrt)^lambda1_wrt - 1) / lambda1_wrt
}



# prep transformations for wSAratio
# Step 1: Extract lambda value for 'wSAratio' from output_df_lsr
lambda1_wSAratio <- output_df_lsr[output_df_lsr$variable == 'wSAratio', 'lambda1']
lambda2_wSAratio <- output_df_lsr[output_df_lsr$variable == 'wSAratio', 'lambda2']

transformation_wSAratio <- if (lambda1_wSAratio == 0) {
  function(x) log(x + lambda2_wSAratio)
} else if (lambda1_wSAratio == 1) {
  function(x) x
} else {
  function(x) ((x + lambda2_wSAratio)^lambda1_wSAratio - 1) / lambda1_wSAratio
}

# prep transformations for yrc/yrp
# Step 1: Extract lambda value for 'yrc' from output_df_lsr
lambda1_yrc <- output_df_lsr[output_df_lsr$variable == 'yrc', 'lambda1']
lambda2_yrc <- output_df_lsr[output_df_lsr$variable == 'yrc', 'lambda2']

transformation_yrc <- if (lambda1_yrc == 0) {
  function(x) log(x + lambda2_yrc)
} else if (lambda1_yrc == 1) {
  function(x) x
} else {
  function(x) ((x + lambda2_yrc)^lambda1_yrc - 1) / lambda1_yrc
}




# Make a copy of the data for prediction sites
pred_df <- pred_mlr_usa_lsr_df_working
if (!dir.exists("predictions/data")) {
  dir.create("predictions/data", recursive = TRUE)
}
# use to predict at all sites
file_path <- paste0("predictions/data/prediction_finaldata_4loop_df_", pred_version, ".csv")
write.csv(pred_df, file=file_path, row.names=FALSE)

pred <- as.matrix(pred_mlr_usa_lsr_df_working)



# prep variables for loop
sid <- SID
yrp <- as.vector(pred[, "yrp"])
yrc <- as.vector(pred[, "yrc"])
tmin <- as.vector(yrp + 1)
tmax <- 2050
yrr <- as.vector(pred[, "yrr"] + 1)
k <- as.vector(pred[, "k"])
DA_mi2 <- as.vector(pred[, "DA_mi2"])
DA_km2 <- as.vector(pred[, "DA_km2"])
trapp <- as.vector((pred[, "trapp"])*100)
capp <- as.vector(pred[, "capp"])
maq_cmy <- as.vector((pred[ , "MAQ_NHDcms"])*31536000)
maq_cms <- as.vector(pred[ , "MAQ_NHDcms"])

model_sums_2050 <- matrix(0, nrow = length(sid), ncol = length(coefficients) + 1)




# Predictions - Model Loop  -----------------------------------------------------

# NEW CODE BEGIN, CONFIDENCE INTERVAL TESTING
### FORWARD LOOP
# In loop, average sedimentation rate calculated from yrc to yr used to update capacity at each timestep from og cap
for (i in seq_along(sid)) {
  yr_values <- tmin[i]:tmax

  if (yrc[i] < yrp[i]) {
    yr_diff <- which(colnames(cap_m)==as.character(yrc[i])):which(colnames(cap_m) == as.character(yrp[i]))
    cap_m[i, yr_diff] <- capp[i]
    cap_up_95_m[i, yr_diff] <- capp[i]
    cap_low_95_m[i, yr_diff] <- capp[i]
    sv_m[i, yr_diff] <- 0
    sv_up_95_m[i, yr_diff] <- 0
    sv_low_95_m[i, yr_diff] <- 0
  }

  # Initialize model_sums_yr for the current site
  model_sums_yr <- matrix(0, nrow = length(yr_values), ncol = length(coefficients) + 1)
  row_index <- 1

  # set up matrices to store average sed rate
  sr_list <- list()

  # Calculate and update capacity values for each year starting from tmin + 1
  for (yr in yr_values) {
    col_name <- as.character(yr)

    if (col_name %in% colnames(cap_m)) {
      prev_yr <- as.character(yr - 1)
      
      # NEW: initialize here
      x_new <- numeric(length(coefficients)+1)                   # NEW TEST
      names(x_new) <- c("(Intercept)", names(coefficients))      # NEW TEST
      x_new["(Intercept)"] <- 1                                  # NEW TEST
      
      # Initialize a matrix to store model sums for the current year
      model_sums <- matrix(0, nrow = 1, ncol = length(coefficients))

      # Iterate over each coefficient for each parameter
      for (j in 1:length(coefficients)) {
        coef_name <- names(coefficients)[j]
        coef_value <- coefficients[j]

        # For coefficient names that are interaction terms, split parameter by the colon to access both columns in "pred" df
        if (grepl(":", coef_name)) { # Check if coefficient name contains a colon
          # If it does, split the coefficient name and extract the two column names
          parts <- unlist(strsplit(coef_name, ":"))
          column1 <- parts[1]
          column2 <- parts[2]

          # Multiply values from each interaction term column by the coefficient
          if (all(c(column1, column2) %in% colnames(pred))) {
            # To apply temporal indicator, check if: column2 is a year, if it's greater than current year, and year starts with '1' or '2'
            if (grepl("^[12]\\d{3}$", column2) && as.numeric(column2) > yr) {
              interaction_product <- 0  # Set interaction product to 0
            } else {
              interaction_product <- pred[i, column1] * pred[i, column2] * coef_value
            }
            x_new[coef_name] <- pred[i, column1] * pred[i, column2]              # NEW TEST
            model_sums[1, j] <- interaction_product
          }

          # For coefficient names whose data values are changing every timestep (wSAatdam and avetrap) and are not stored in "pred" df
        } else {
          # Fetch wSAatdam value for the current year (time-weighted and transformed above) and multiply by coefficient
          if (coef_name == "wSAatdam" && col_name %in% colnames(wSAatdam_m_transformed)) {
            x_new[coef_name] <- wSAatdam_m_transformed[i, col_name]            # NEW TEST
            model_sums[1, j] <- wSAatdam_m_transformed[i, col_name] * coef_value

            # Fetch avetrap value from the previous year
          } else if (coef_name == "AveTrap") {
            if (prev_yr %in% colnames(avetrap_m_transformed)) {
              x_new[coef_name] <- avetrap_m_transformed[i, prev_yr]           # NEW TEST
              model_sums[1, j] <- avetrap_m_transformed[i, prev_yr] * coef_value
            } else {
              x_new[coef_name] <- 0                                            # NEW TEST
              model_sums[1, j] <- 0  # Set to 0 if column does not exist (which shouldn't happen)
            }

            # create yrs_tot value in timeloop and transform
          } else if (coef_name == 'yrs_tot') {
            yrs_tot_value <- yr - yrp[i]
            yrs_tot_transformed <- transformation_yrstot(yrs_tot_value)
            x_new[coef_name] <- yrs_tot_transformed                            # NEW TEST
            model_sums[1, j] <- yrs_tot_transformed * coef_value

            # create yrp value in timeloop and transform
          } else if (coef_name == 'yrp') {
            yrp_value <- yrp[i]
            yrp_transformed <- transformation_yrc(yrp_value)
            x_new[coef_name] <- yrp_transformed                               # NEW TEST
            model_sums[1, j] <- yrp_transformed * coef_value

            # create wrt value in timeloop and transform
          } else if (coef_name == 'wrt') {
            cap2 <- cap_m[i, prev_yr]
            mean_cap <- (capp[i] + cap2) / 2

            if (is.na(maq_cmy[i]) || is.na(mean_cap)) {
              # If maq_cmy is NA, set model_sums to 0
              wrt_value_transformed <- 0
            } else {
              # Calculate wrt_value and transform it if maq_cmy is not NA
              wrt_value <- mean_cap / maq_cmy[i]
              wrt_value_transformed <- transformation_wrt(wrt_value)
            }
            # Update model_sums with the transformed value or 0
            x_new[coef_name] <- wrt_value_transformed                           # NEW TEST
            model_sums[1, j] <- wrt_value_transformed * coef_value


            # create wSAratio value in timeloop and transform
          } else if (coef_name == 'wSAratio') {
            wSAatdam_value <- wSAatdam_m[i, col_name]
            wSAratio_value <- (wSAatdam_value/DA_km2[i])*100

            if (is.na(wSAratio_value)) {
              # If wSAratio_value is NA, set model_sums to 0
              wSAratio_value_transformed <- 0
            } else {
              # Transform wSAratio
              wSAratio_value_transformed <- transformation_wSAratio(wSAratio_value)
            }

            # Update model_sums with the transformed value or 0
            x_new[coef_name] <- wSAratio_value_transformed                   # NEW TEST
            model_sums[1, j] <- wSAratio_value_transformed * coef_value


            # For coefficient names whose data values are constant through time
          } else {
            # Fetch data value from 'pred' df
            if (coef_name %in% colnames(pred)) {
              x_new[coef_name] <- pred[i, coef_name]                        # NEW TEST 
              model_sums[1, j] <- pred[i, coef_name] * coef_value
            }
          }
        }
      }

      # Store the model sums and intercept for the current year
      model_sums_yr[row_index, ] <- c(model_sums, intercept)
      row_index <- row_index + 1
      
      # NEW
      fit_log <- sum(x_new * coef_beta)                              # NEW TEST
      se_mean <- sqrt(t(x_new)%*% vcov_beta %*% x_new)               # NEW TEST
      se_pred <- sqrt(se_mean^2 + rse^2)                             # NEW TEST

      if (yr <= yrr[i]) {

        # Predict lSR using annual MLR model, update current year's lsr_m
        lsr <- intercept + sum(model_sums[1,])
        # lsr_up_95 <- lsr + (t_critical_95 * upper_bound_95)
        # lsr_low_95 <- lsr - (t_critical_95 * lower_bound_95)
        lsr_up_95 <- lsr + (t_critical_95 * se_pred)                # NEW TEST
        lsr_low_95 <- lsr - (t_critical_95 * se_pred)               # NEW TEST

        sr <- exp(lsr)
        sr_up_95 <- exp(lsr_up_95)
        sr_low_95 <- exp(lsr_low_95)

        # Save sed rate values and CI to output df
        lsr_m[i, col_name] <- lsr
        sr_m[i, col_name] <- sr

        sr_up_95_m[i, col_name] <- sr_up_95
        sr_low_95_m[i, col_name] <- sr_low_95

        # append new sedimentation rate to list of annualized sed rates
        sr_list <- c(sr_list, sr)

        # Track cumulative relative uncertainty at 95% confidence level
        ru_sr_up_95 <- sr_up_95 - sr
        ru_sr_low_95 <- sr - sr_low_95

        # Update the current year's capacity
        sv_new <- sr*length(sr_list)
        cap_new <- capp[i]-sv_new

        sv_up_95_new <- sv_new + (ru_sr_up_95*length(sr_list))
        sv_low_95_new <- sv_new - (ru_sr_low_95*length(sr_list))

        cap_up_95_new <- capp[i] - sv_low_95_new  # reversed because upper bound on capacity would be the lesser SV value
        cap_low_95_new <- capp[i] - sv_up_95_new  # ditto

        # second check to limit uncertainty if capacity uncertainty goes negative
        if(cap_low_95_new < 0) {
          cap_low_95_new <- 0
        }

        if(sv_up_95_new > capp[i]) {
          sv_up_95_new <- capp[i]
        }

        sv_m[i, col_name] <- sv_new
        sv_up_95_m[i, col_name] <- sv_up_95_new
        sv_low_95_m[i, col_name] <- sv_low_95_new

        cap_m[i, col_name] <- cap_new
        cap_up_95_m[i, col_name] <- cap_up_95_new
        cap_low_95_m[i, col_name] <- cap_low_95_new

        if (cap_new > 0 && cap_new < cap_m[i, prev_yr]) {
          # Update the current year's trap efficiency based on the new capacity
          cap_acft <- cap_new/1233.482                                          # convert cap in m3 to ac-ft
          trap <- (1-1/(1+k[i]*(cap_acft/DA_mi2[i])))*100
          trap_m[i, col_name] <- trap
          trap_m_transformed[i, col_name] <- transformation_avetrap(trap)

          # Calculate simple average trap efficiency between the first timestep and current timestep
          avetrap <- (trapp[i] + trap)/2
          avetrap_m[i, col_name] <- avetrap
          avetrap_m_transformed[i, col_name] <- transformation_avetrap(avetrap)

        } else {
          # Set all subsequent years' capacity to 0 for this site since all cap lost
          col_indices <- which(colnames(cap_m) == col_name):which(colnames(cap_m) == as.character(tmax))
          cap_m[i, col_indices] <- 0
          cap_up_95_m[i, col_indices] <- 0
          cap_low_95_m[i, col_indices] <- 0
          lsr_m[i, col_indices] <- 0
          sr_m[i, col_indices] <- 0
          trap_m[i, col_indices] <- 0
          trap_m_transformed[i, col_indices] <- 0
          sv_m[i, col_indices] <- capp[i]
          sv_up_95_m[i, col_indices] <- capp[i]
          sv_low_95_m[i, col_indices] <- capp[i]

          avetrap_m[i, col_indices] <- 0
          avetrap_m_transformed[i, col_indices] <- 0

          break  # Exit the loop for this site
        }

        # if yr is not less than yrr, then dam was removed and need to release sed storage and cap
      } else {
        col_indices <- which(colnames(cap_m) == col_name):which(colnames(cap_m) == as.character(tmax))
        cap_m[i, col_indices] <- 0
        cap_up_95_m[i, col_indices] <- 0
        cap_low_95_m[i, col_indices] <- 0
        lsr_m[i, col_indices] <- 0
        sr_m[i, col_indices] <- 0
        trap_m[i, col_indices] <- 0
        trap_m_transformed[i, col_indices] <- 0
        sv_m[i, col_indices] <- 0
        sv_up_95_m[i, col_indices] <- 0
        sv_low_95_m[i, col_indices] <- 0

        avetrap_m[i, col_indices] <- 0
        avetrap_m_transformed[i, col_indices] <- 0

        break  # Exit the loop for this site
      }
    }
    # clear temporary values explicitly to ensure not being carried through
    lsr <- NULL
    lsr_up <- NULL
    lsr_low <- NULL
    sr <- NULL
    sr_up_95 <- NULL
    sr_low_95 <- NULL
    prev_yr <- NULL
    trap <- NULL
    avetrap <- NULL
    cap_new <- NULL
    cap_up_95_new <- NULL
    cap_low_95_new <- NULL
    cap_acft <- NULL
    sv_new <- NULL
    sv_up_95_new <- NULL
    sv_low_95_new <- NULL
    ru_sr_up <- NULL
    ru_sr_low <- NULL
    col_indices <- NULL
    yrs_tot_value <- NULL
    yrs_tot_transformed <- NULL
    yrp_value <- NULL
    yrp_transformed <- NULL
    cap2 <- NULL
    mean_cap <- NULL
    wrt_value <- NULL
    wrt_value_transformed <- NULL
    wSAratio_value <- NULL
    wSAratio_value_transformed <- NULL
    
    x_new <- NULL                                             # NEW TEST
    
  }
  # Append the final model_sums for the current site to model_sums_2050
  model_sums_2050[i, ] <- c(model_sums, intercept)
}
# NEW CODE END, CONFIDENCE INTERVAL TESTING


# Check how many sites have 0 remaining capacity for each year
zero_counts <- sapply(colnames(cap_m), function(col) sum(cap_m[, col] == 0, na.rm = TRUE))
print(zero_counts)



### CUMULATIVE SUMS PLOTS
## SEDIMENT VOLUME
# Sum up the sediment volumes for each year
year_sums <- (colSums(sv_m[, -(1:2)], na.rm = TRUE))/1e9
error_up_95_sums <- (colSums(sv_up_95_m[, -(1:2)], na.rm=TRUE))/1e9
error_low_95_sums <- (colSums(sv_low_95_m[, -(1:2)], na.rm=TRUE))/1e9
years <- as.numeric(names(year_sums))
year_sum_m <- data.frame(Year = years, Sum = year_sums, 
                         Error_Up_95 = error_up_95_sums, Error_Low_95 = error_low_95_sums)
# specify file path
file_path <- "figures/predictions/cumulative_sums"
dir.create(file_path, recursive = TRUE, showWarnings=FALSE)

# save plot to file
plot_filename <- paste(file_path, paste("cumsum_sed_", pred_version, ".png", sep=""), sep="/")
png(plot_filename, width=1400, height=1200, res=200)

ggplot(year_sum_m, aes(x = Year, y = Sum)) +
  geom_ribbon(aes(ymin = Error_Low_95, ymax = Error_Up_95), fill="grey70", alpha=0.3) +
  geom_line() +
  geom_point() +
  labs(
    title = paste("Predicted cumulative sedimentation (n =", nrow(sv_m), ")", sep= ""),
    x = "Year",
    y = bquote("Sediment volume (km"^3*")")
  ) +
  theme_minimal()
dev.off()

# print plot to screen
ggplot(year_sum_m, aes(x = Year, y = Sum)) +
  geom_ribbon(aes(ymin = Error_Low_95, ymax = Error_Up_95), fill="grey70", alpha=0.3) +
  geom_line() +
  geom_point() +
  labs(
    title = paste("Predicted cumulative sedimentation (n =", nrow(sv_m), ")", sep= ""),
    x = "Year",
    y = bquote("Sediment volume (km"^3*")")
  ) +
  theme_minimal()


## CAPACITY
year_sums <- (colSums(cap_m[, -(1:2)], na.rm = TRUE))/1e9
error_up_95_sums <- (colSums(cap_up_95_m[, -(1:2)], na.rm=TRUE))/1e9
error_low_95_sums <- (colSums(cap_low_95_m[, -(1:2)], na.rm=TRUE))/1e9
years <- as.numeric(names(year_sums))
year_sum_m <- data.frame(Year = years, Sum = year_sums, 
                         Error_Up_95 = error_up_95_sums, Error_Low_95 = error_low_95_sums)

# save plot to file
plot_filename <- paste(file_path, paste("cumsum_cap_", pred_version, ".png", sep=""), sep="/")
png(plot_filename, width=1400, height=1200, res=200)
ggplot(year_sum_m, aes(x = Year, y = Sum)) +
  geom_ribbon(aes(ymin = Error_Low_95, ymax = Error_Up_95), fill="grey70", alpha=0.3) +
  geom_line() +
  geom_point() +
  labs(
    title = paste("Predicted cumulative capacity (n =", nrow(cap_m), ")", sep= ""),
    x = "Year",
    y = bquote("Capacity (km"^3*")")
  ) +
  theme_minimal()
dev.off()

# print plot to screen
ggplot(year_sum_m, aes(x = Year, y = Sum)) +
  geom_ribbon(aes(ymin = Error_Low_95, ymax = Error_Up_95), fill="grey70", alpha=0.3) +
  geom_line() +
  geom_point() +
  labs(
    title = paste("Predicted cumulative capacity (n =", nrow(cap_m), ")", sep= ""),
    x = "Year",
    y = bquote("Capacity (km"^3*")")
  ) +
  theme_minimal()


### FLAG 
# Flag any reservoirs whose capacity ever increases rather than decreases through time
# Convert column names to integers for easier comparison
year_cols <- colnames(cap_m)
years <- as.integer(gsub("X", "", year_cols))

# Initialize a list to store SIDs where the values are growing over time after yrp
flagged_sids <- list()

# Iterate through each row in the matrix
for (i in seq_len(nrow(cap_m))) {
  # Get the year from 'yrp' column
  yrp <- pred[, "yrp"][i]

  # Find the columns corresponding to years greater than yrp
  year_cols_after_yrp <- year_cols[years > yrp]

  if (length(year_cols_after_yrp) > 0) {
    # Extract the values for the timeseries after yrp
    values <- as.numeric(cap_m[i, year_cols_after_yrp])

    # Check if the values are monotonically increasing
    if (length(values) > 1 && any(diff(values) > 0)) {
      # Append the SID to the flagged_sids list if the condition is met
      flagged_sids <- c(flagged_sids, rownames(cap_m)[i])
    }
  }
}

# Convert the list to a dataframe for easier handling
flagged_sids_m <- data.frame(SID = unlist(flagged_sids), stringsAsFactors = FALSE)
print(flagged_sids_m)

# Convert trap efficiency to a decimal percent
trap_decimal_m <- trap_m/100


### SAVE RESULT DATAFRAMES

# specify file path
folder <- "predictions/mlr_outputs/"
dir.create(folder, showWarnings=FALSE)

# save outputs
file_path <- paste0(folder, "cap_", pred_version, ".csv")
write.csv(cap_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "cap_up_95_", pred_version, ".csv")
write.csv(cap_up_95_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "cap_low_95_", pred_version, ".csv")
write.csv(cap_low_95_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "trap_", pred_version, ".csv")
write.csv(trap_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "trap_decimal_", pred_version, ".csv")
write.csv(trap_decimal_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "trap_transformed_", pred_version, ".csv")
write.csv(trap_m_transformed, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "avetrap_", pred_version, ".csv")
write.csv(avetrap_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "avetrap_transformed_", pred_version, ".csv")
write.csv(avetrap_m_transformed, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "lsr_", pred_version, ".csv")
write.csv(lsr_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "sr_", pred_version, ".csv")
write.csv(sr_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "sr_up_95_", pred_version, ".csv")
write.csv(sr_up_95_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "sr_low_95_", pred_version, ".csv")
write.csv(sr_low_95_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "sv_", pred_version, ".csv")
write.csv(sv_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "sv_up_95_", pred_version, ".csv")
write.csv(sv_up_95_m, file=file_path, row.names=TRUE)

file_path <- paste0(folder, "sv_low_95_", pred_version, ".csv")
write.csv(sv_low_95_m, file=file_path, row.names=TRUE)


## END!
