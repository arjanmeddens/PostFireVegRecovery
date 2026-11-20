####################################################################
# Purpose: RandomForest models with selected variables
####################################################################
# Author: Madeline Franz (and Arjan Meddens) 
# Date: September 24, 2026
####################################################################
# Input format: Request from authors
# Output format: Figures, GeoTiffs
# R version: Meddens: R version 4.5.1 (2025-06-13) -- "Great Square Root" 
#            Franz: 4.4.1
####################################################################
#### ---- 1.1 Initiate run and set working directory ---- ####
rm(list = ls()); gc()
ptm <- Sys.time()
message(paste(Sys.time()))


####################################################################
### Title CONIFER REGENERATION (%) #####
####################################################################
#fig_dir  = "/Users/arjanmeddens/OneDrive - Washington State University (email.wsu.edu)/Meddens_Franz/1_Landsat_analysis/1_Figures/"
#data_dir = "/Users/arjanmeddens/OneDrive - Washington State University (email.wsu.edu)/Meddens_Franz/1_Landsat_analysis/4_Data_out/"
fig_dir  = "/Users/arjanmeddens/OneDrive - Washington State University (email.wsu.edu)/MeddensLab/1_Franz_etal_Chapter1/workspace/figures/"
data_dir  = "/Users/arjanmeddens/OneDrive - Washington State University (email.wsu.edu)/MeddensLab/1_Franz_etal_Chapter1/Field Data/"
setwd("/Users/arjanmeddens/OneDrive - Washington State University (email.wsu.edu)/MeddensLab/1_Franz_etal_Chapter1/workspace/")

#### ---- 1.2 Set seed and directory parameters --- ####
var = 3 # 1 <> SHRUB (seed=5), 2 <> SEEDLING (seed=1), 3 <> GRASS (seed=15), 4 <> HERB, 5 <> GRASSHERB ...
train.proportion = 0.7
seed = 15; set.seed(seed)
run_id = paste("_s",seed,sep=", ")
varimp = 0
#plot_img = 1
write_out = 1
run_reduced_predictors = 0

#### ---- 1.3 Load packages --- ####
source("~/workspace/R/R_functions/model.stats.R") # see also in our shared Workspace directory  
source("~/workspace/R/R_functions/RFE_Helper_Functions_240904.R") 
package.list <- c('hrbrthemes','ggplot2','grid','here', 'terra', 'tidyverse', 'viridis','ggsci', 
                  'rcompanion','tidymodels','themis','knitr', 'beepr','skimr','iml', 'sf', 
                  'spatialsample','doParallel','pROC', 'randomForest', 'corrplot', 'factoextra', 
                  'cluster', 'DALEXtra', 'dials', 'RColorBrewer', 'tidyterra', 'rsample','vip',
                  'rsample','recipes','ggspatial','parsnip','ranger','scico') # 'GGally',


## Install them if they aren't already on the computer
#new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
## And load them to be able to use some of the functions
for(i in package.list){library(i, character.only = T)}
#- AM 260825: Needed vip from Remotes--- R-version 4.5 
#install.packages("remotes")
#remotes::install_version("vip", version = "0.4.1", repos = "https://cloud.r-project.org")
#remotes::install_github("koalaverse/vip")

######### ---------- Part 2: Initial Data Import --- ########
#### ---- 2.1 Import and check field data ---- ####

### Import field data
data.tab = read.csv(paste0(data_dir,"AllBluesUnderstoryData_Trim.csv"))
#head(data.tab)

### Turn field data into a spatial vector and reproject into "epsg:32611"
init_spat_dat = vect(data.tab, geom=c("Longitude", "Latitude"))
crs(init_spat_dat) = ("epsg:4326")
spat_dat = project(init_spat_dat, "epsg:32611")

#import the completed all-fire raster stack
all_rast = rast(paste0(data_dir,"all_rast.tif"))
all_rast1 = all_rast

######### ---------- Part 3: Raster Data Extraction and Prep --- ########
#### ---- 3.1 Extract raster data for plot locations ---- ####

# Select only the layers that will be used
#paste(names(all_rast1), collapse = ', ')

all_rast = subset(all_rast1, c(1:10, 12:14, 16:17, 34, 36:38, 40:42, 44:45, 52:57, 62, 64:71)) #Select terrain indices, Trajectory metrics from these indices: NBR, TCW, and TCB, and Phenological metrics from these indices: NDMI drydif, NDMI dryrat, NDVI dryrat, NDVI snowdif, NDVI snowrat, and TCW dryrat)

# Change VarImp abbreviations here 
if (varimp == 1) {
  all_predictors_arr = c("ELEV","SLOPE","CONSTANT","PF_NBR_GS","SEV_NBR_GS","FY_NBR_GS",
                    "ER_NBR_GS","LR_NBR_GS","RR15_NBR_GS","FY_NDMI_DD","ER_NDMI_DD",
                    "LR_NDMI_DD","FY_NDMI_DR","ER_NDMI_DR","LR_NDMI_DR","FY_NDVI_DR",
                    "ER_NDVI_DR","LR_NDVI_DR","FY_NDVI_SND","ER_NDVI_SND","LR_NDVI_SND",
                    "FY_NDVI_SNR","ER_NDVI_SNR","LR_NDVI_SNR","PF_TCB_GS","SEV_TCB_GS",
                    "FY_TCB_GS","ER_TCB_GS","LR_TCB_GS","RR15_TCB_GS","FY_TCW_DR","ER_TCW_DR",
                    "LR_TCW_DR","PF_TCW_GS","SEV_TCW_GS","FY_TCW_GS","ER_TCW_GS","LR_TCW_GS","RR_15_TCW_GS")
} else {
  all_predictors_arr = names(all_rast) 
}  
names(all_rast) = all_predictors_arr

fire_yx = terra::extract(all_rast, spat_dat, bind=TRUE)
fire_dat_na = as.data.frame(fire_yx, geom='XY') #176 obs
fire_dat_narm = na.omit(fire_dat_na) #158 obs

if (varimp == 1) {
  fire_dat = fire_dat_narm[fire_dat_narm$SEV_NBR_GS < -0.02,] ##SEVERITY Threshold!
} else {
  fire_dat = fire_dat_narm[fire_dat_narm$sev_NBR_GS < -0.02,] ##SEVERITY Threshold!
}

fire_dat = fire_dat %>% rowwise() %>% mutate(GrassHerb.pct = sum(c(grass.pct, herb.pct)))

#### ---- 3.1 Set formulas for modeling ---- ####

if (var == 1) {out_text = c("shrub") 
text.fig = "(B)"
text.y   = "Shrub cover (%)"
select_strata <- c('shrub.pct')
y_var = c("shrub.pct")
# SHRUB #
if (run_reduced_predictors == 1) {load(file=paste(fig_dir,"selected_index_",out_text,"_s",seed,".rdata", sep=''),verbose = T)
  select_formula <- as.formula(paste0("shrub.pct ~", paste0(all_predictors_arr[index_indices], collapse="+"))) 
  vir_text = "selpred"
} else {
  select_formula <- as.formula(paste0("shrub.pct ~", paste0(names(all_rast), collapse="+")))
  vir_text = "allpred"
}}



if (var == 2) {out_text = c("seedling") 
text.fig = "(A)"
text.y   = "Seedling cover (%)"
select_strata <- c('seedling.pct')
y_var = c("seedling.pct")
# SEEDLING #
if (run_reduced_predictors == 1) {load(file=paste(fig_dir,"selected_index_",out_text,"_s",seed,".rdata", sep=''),verbose = T)
  select_formula <- as.formula(paste0("seedling.pct ~", paste0(all_predictors_arr[index_indices], collapse="+"))) 
  vir_text = "selpred"
} else {
  select_formula <- as.formula(paste0("seedling.pct ~", paste0(names(all_rast), collapse="+")))
  vir_text = "allpred"
}
}
if (var == 3) {out_text = c("grass") 
text.fig = "(C)"
text.y   = "Grass cover (%)"
select_strata <- c('grass.pct')
y_var = c("grass.pct")
# GrassHerb.pct #
if (run_reduced_predictors == 1) {load(file=paste(fig_dir,"selected_index_",out_text,"_s",seed,".rdata", sep=''),verbose = T)
  select_formula <- as.formula(paste0("grass.pct ~", paste0(all_predictors_arr[index_indices], collapse="+"))) 
  vir_text = "selpred"
} else {
  select_formula <- as.formula(paste0("grass.pct ~", paste0(names(all_rast), collapse="+")))
  vir_text = "allpred"
}
}
if (var == 4) {out_text = c("herb") 
text.fig = "(D)"
text.y   = "Herb cover (%)"
select_strata <- c('herb.pct')
y_var = c("herb.pct")
if (run_reduced_predictors == 1) {load(file=paste(fig_dir,"selected_index_",out_text,"_s",seed,".rdata", sep=''),verbose = T)
  select_formula <- as.formula(paste0("herb.pct ~", paste0(all_predictors_arr[index_indices], collapse="+"))) 
  vir_text = "selpred"
} else {
  select_formula <- as.formula(paste0("herb.pct ~", paste0(names(all_rast), collapse="+")))
  vir_text = "allpred"
}
}
if (var == 5) {out_text = c("grassHerb") 
text.fig = "(D)"
text.y     = "Grass & Herb cover (%)"
select_strata <- c('GrassHerb.pct')
y_var = c("GrassHerb.pct")
if (run_reduced_predictors == 1) {load(file=paste(fig_dir,"selected_index_",out_text,"_s",seed,".rdata", sep=''),verbose = T)
  select_formula <- as.formula(paste0("GrassHerb.pct ~", paste0(all_predictors_arr[index_indices], collapse="+"))) 
  vir_text = "selpred"
} else {
  select_formula <- as.formula(paste0("GrassHerb.pct ~", paste0(names(all_rast), collapse="+")))
  vir_text = "allpred"
}
}

out_text <- paste0(out_text,"v4")

#export
#write.csv(fire_dat, here('Data_Repository/Working/Midpoints','allFire_RF_dat.csv'), row.names = F)

######### ---------- Part 4: Random Forest Modeling --- ########
# NOTE I DELETED A BUNCH OF CODE HERE

#### ---- 4.1 Determine spatial clusters for random forest ---- ####
## Optimal # of clusters is 6!

#### ---- 4.2 Partition data into train/test and create spatial cross-validation splits ---- ####
#fire_dat_import <- read.csv('Data_Repository/Working/Midpoints/allFire_RF_dat.csv')

set.seed(seed)
dat_split <- initial_split(fire_dat, prop=train.proportion, strata=all_of(select_strata)) # See PARAM SELECTION up Top
dat_train <- training(dat_split) # create object with training data
dat_test <- testing(dat_split) # create object with testing data

# Partition training data to allow for spatial cross validation (again stratified by the response class)
set.seed(seed)

spat_dat_train <- sf::st_as_sf(dat_train, coords = c("x","y"), crs = st_crs("EPSG:32611"))
train_cv <- spatial_clustering_cv(spat_dat_train, cluster_function = "hclust", v=6)#6-fold spatial cross-validation. We split into 6 clusters, per results of earlier k-means. Not repeated b/c results do not vary.
#autoplot(train_cv)

#spat_dat_train <- sf::st_as_sf(dat_train, coords = c("x","y"), crs = st_crs("EPSG:32611"))
#train_cv <- spatial_clustering_cv(spat_dat_train, cluster_function = "hclust", v=6)#6-fold spatial cross-validation. We split into 6 clusters, per results of earlier k-means. cluster number determination not repeated b/c results do not vary.
#train_cv
#autoplot(train_cv)

#### ---- 4.3 Random Forest Recursive Feature Elimination and Tuning ---- ####
model_formula <- select_formula
model_formula

# This uses the custom rangerFuncs funcitons from RFE_Helper_Functions script.
# It is currently running 500 trees. for the absolute final run, I should change the post script from _500 to _1000 or _10000 as I've coded those num.tree options into the function script.
#rangerFuncs

# SHRUB #
if (var == 1) {ctrl <- rfeControl(functions = rangerFuncs_05_500, #random forest
                                  method = "repeatedcv", #repeated cv
                                  returnResamp = "all", 
                                  repeats = 5, # number of repeats. may want to increase to 10 for final run.
                                  number = 5, # number of folds
                                  rerank = TRUE,
                                  verbose = TRUE,
                                  saveDetails = TRUE)
}

# SEEDLING #
if (var == 2) {ctrl <- rfeControl(functions = rangerFuncs_best_500,
                                  method = "repeatedcv",
                                  returnResamp = "all",
                                  repeats = 5,
                                  number = 5,
                                  rerank = TRUE,
                                  verbose = TRUE,
                                  saveDetails = TRUE)
}

# GRASS #
if (var == 3) {ctrl <- rfeControl(functions = rangerFuncs_best_500,
                                  method = "repeatedcv",
                                  returnResamp = "all",
                                  repeats = 5,
                                  number = 5,
                                  rerank = TRUE,
                                  verbose = TRUE,
                                  saveDetails = TRUE)
}

# HERB #
if (var == 4) {ctrl <- rfeControl(functions = rangerFuncs_best_500,
                                  method = "repeatedcv",
                                  returnResamp = "all",
                                  repeats = 5,
                                  number = 5,
                                  rerank = TRUE,
                                  verbose = TRUE,
                                  saveDetails = TRUE)
}

# GRASSHERB #
if (var == 5) {ctrl <- rfeControl(functions = rangerFuncs_best_500,
                                  method = "repeatedcv",
                                  returnResamp = "all",
                                  repeats = 5,
                                  number = 5,
                                  rerank = TRUE,
                                  verbose = TRUE,
                                  saveDetails = TRUE)
}

# RFE
subsets <- c(1:39) #b/c only have 39 variables.
rfe_tune_rec <- recipe(model_formula, dat_train) %>%
  step_normalize(all_numeric(), -all_outcomes()) %>% 
  step_zv(all_numeric(), -all_outcomes())
scaled_dat_train <- prep(rfe_tune_rec, training = dat_train, retain = TRUE) %>% 
  juice() #scaling the variables! need to take this into account for accuracy measures! 
# actually, run the rfe like this, then run once more with selected model.

set.seed(seed)

#This next line is the (relative) slowest--it can take a minute or two (at 500 trees). 
rfProfile <- rfe(model_formula, 
#                 data = scaled_dat_train, 
                data = dat_train,
                 metric = "RMSE",
                 maximize = FALSE,
                 sizes = subsets, 
                 rfeControl = ctrl)

#rfProfile 
plot(rfProfile, type = c("g", "o"))
#densityplot(rfProfile)
xyplot(rfProfile)
#stripplot(rfProfile)
#histogram(rfProfile)
#xyplot(rfProfile$results$RMSE ~ rfProfile$results$Variables,
#       type = c("g", "p", "l"),
#       auto.key = TRUE)
rfProfile$variables
new_x <- predictors(rfProfile)
updated_x <- paste(new_x, collapse=" + ")
updated_x

rfProfile$fit$predictions #predicted values for each plot
rfProfile$fit$mtry
rfProfile$fit$min.node.size
rfProfile$fit$variable.importance
rfProfile$fit$prediction.error
rfProfile$fit$r.squared
vip(rfProfile$fit)

# SHRUB #
if (var == 1) {
  updated_formula <- as.formula(paste0("shrub.pct ~ ", updated_x))
}

# SEEDLING #
if (var == 2) {
  updated_formula <- as.formula(paste0("seedling.pct ~ ", updated_x))
}

# GRASS #
if (var == 3) {
  updated_formula <- as.formula(paste0("grass.pct ~ ", updated_x))
}

# HERB #
if (var == 4) {
  updated_formula <- as.formula(paste0("herb.pct ~ ", updated_x))
}

# GRASSHERB #
if (var == 5) {
  updated_formula <- as.formula(paste0("GrassHerb.pct ~ ", updated_x))
}

set.seed(seed)

# Set RF formula
rfe_rec <- recipe(updated_formula, dat_train) %>%
  step_normalize(all_numeric(), -all_outcomes()) %>% 
  step_zv(all_numeric(), -all_outcomes()) 

# Set model
rfe_spec <- rand_forest() %>%
  set_args(#min_n=5, # from parsnip: 'An integer for the minimum number of data points in a node that are required for the node to be split further.' and from ranger: 'min.node.size: Minimal node size to split at. Default 1 for classification, 5 for regression, 3 for survival, and 10 for probability.'
    #           mtry=tune(), #'Number of variables to possibly split at in each node. Default is the (rounded down) square root of the number variables. Alternatively, a single argument function returning an integer, given the number of independent variables' (from ranger) or (from parsnip) 'An integer for the number of predictors that will be randomly sampled at each split when creating the tree models' default is floor(sqrt(ncol(x)))
    trees=1000) %>% # 'An integer for the number of trees contained in the ensemble.'
  set_engine("ranger", seed = seed, importance = "permutation") %>%  
  set_mode("regression")

# Set the workflow
rfe_workflow <- workflow() %>%
  add_recipe(rfe_rec) %>%
  add_model(rfe_spec)

# Tune model
set.seed(seed)

#run initial model
rfe_res <- 
rfe_workflow %>%
  tune_grid(
    resamples = train_cv,
  )

#Check best-performing models:
#show_best(tune_res, metric = "rmse")
#autoplot(tune_res)

final_rf <- rfe_workflow %>%
  finalize_workflow(select_best(rfe_res, metric = "rmse"))

final_rf

dat_fit <- last_fit(final_rf, dat_split) #final fit of best model on training set, evaluated on test sites
#dat_fit_rmse <- last_fit(final_rf, dat_split)
dat_fit <- last_fit(rfe_workflow, dat_split)



#### ---- 4.4 Evaluate and visualize model performance ---- ####

# check predictors used
rf_features <- 
  dat_fit %>% 
  extract_fit_engine() %>% 
  predictors()  #<- the caret function

rang_info <-
  dat_fit %>% 
  extract_fit_engine()
rang_info
dat_fit$.workflow

#head(rf_features)
#rf_features

# performance stats
collect_metrics(dat_fit) #RMSE and RSQ on test set
#dat_fit$.metrics #result should be the same


# xy plot of observations vs. predictions
collect_predictions(dat_fit) %>% #also on test set
  ggplot(aes(get(y_var), .pred)) +
  geom_abline(lty = 2, color = "gray50") +
  geom_point(alpha = 0.5, color = "midnightblue") +
  coord_fixed()

final_fitted <- dat_fit$.workflow[[1]]
#predict(final_fitted, dat_test) #print numerical predictions for test data

# Variable importance plot
#imp_spec <- tune_spec %>%
#  finalize_model(select_best(tune_res, metric = "rmse")) %>%
#  set_engine("ranger", importance = "permutation")
#quartz() # window()
final_fitted %>%
  extract_fit_parsnip() %>%
  vip(num_features = length(rf_features)+1, aesthetics = list(alpha = 0.8, fill = "midnightblue"))
#str(final_fitted)

vip(final_fitted,aesthetics = list(alpha = 0.8, fill = "midnightblue")) +
   ylab("Mean Importance Score") +
   xlab("Selected variables\n") + 
   theme(axis.text.x = element_text(size = 14))  + 
   theme(axis.text.y = element_text(size = 14))  +
   theme(axis.title.x = element_text(size = 16)) +
   theme(axis.title.y = element_text(size = 16)) + theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(paste(fig_dir,"VIRplot_",out_text,"_s",seed,"_",vir_text,"_241020.png",sep=''))


# Extract variable importance scores
vip_arr = dat_fit %>% 
  extract_fit_engine() %>% 
  ranger::importance() %>% 
  sort(decreasing = TRUE)
vip_arr

p_vip_arr = dat_fit %>% 
  extract_fit_engine() %>% 
  ranger::importance_pvalues(method = "altmann", formula = model_formula, data = dat_train, num.permutations = 50)
p_vip_arr

# -------------------------------------------------------------------------

names(final_fitted$pre$mold$predictors)
length(names(all_rast))
length(rf_features)
paste(rf_features, collapse = ', ')

a = extract_fit_parsnip(dat_fit)
b = extract_recipe(dat_fit)

#index_list_arr = names(vip_arr)
#predictor.include<-names(vip_arr[1:n.pred.seedlings])
#index_indices = which(all_predictors_arr %in% predictor.include)


# Save the set of predictors to run parsimonious run
#save(index_indices,vip_arr,all_predictors_arr,file=paste(fig_dir,"selected_index_",out_text,"_s",seed,".rdata", sep=''))

## Fire info
fire_abbrv = list('SF', 'BC', 'HR', 'T', 'WP', 'EA', 'EG', 'C', 'F')

SF_full = 'School'
BC_full = 'Bridge Creek'
HR_full = 'Hash Rock'
T_full = 'Tower'
WP_full = 'Wheeler Point'
EA_full = 'Easy'
EG_full = 'Egley'
C_full = 'Calamity'
F_full = 'Flagtail'

SF_year = 2005
BC_year = 2001
HR_year = 2000
T_year = 1996
WP_year = 1996
EA_year = 2002
EG_year = 2007
C_year = 2007
F_year = 2002

SF_year_samp = 2021
BC_year_samp = 2022
HR_year_samp = 2022
T_year_samp = 2022
WP_year_samp = 2023
EA_year_samp = 2023
EG_year_samp = 2023
C_year_samp = 2023
F_year_samp = 2023

#import the completed all-fire raster stack
#all_rast = rast("Data_Repository/Working/Spectral/Fire_MultIndex_RasterStacks/all_rast.tif")


#### ---- 2.3 MTBS parameters and import ---- ####

## list of the MTBS data to import
## MTBS path endings (for the possible layers: "perim", "dnbr", "rdnbr", "sev")
shp_list = list("perim")
tif_list = list("dnbr", "rdnbr", "sev")

perim_path = "_pre_post_burn_bndy"
dnbr_path  = "_pre_post_dnbr"
rdnbr_path = "_pre_post_rdnbr"
sev_path   = "_pre_post_dnbr6"


#- Load rasters

for (i in 1:length(fire_abbrv)) {
# import fire perimeters and severity rasters
#for (fire in fire_abbrv){
  fire_year = get(paste(fire_abbrv[i], "_year", sep = ""))
  full_name = get(paste(fire_abbrv[i], "_full", sep = ""))
  MTBS_path = paste0(data_dir, toString(fire_year), "/", full_name, "/", sep = "")
  if(length(tif_list) > 0)
    for (layer in tif_list){
      layer_name = paste(fire_abbrv[i], "_", layer, sep = "")
      layer_path_name = paste(full_name, get(paste(layer, "_path", sep = "")), sep = "")
      layer_path = paste(MTBS_path, layer_path_name, ".tif", sep = "")
      imported_layer = rast(layer_path)
      assign(layer_name, imported_layer, envir = .GlobalEnv)
    }
  if(length(shp_list) > 0)
    for (layer in shp_list){
      layer_name = paste(fire_abbrv[i], "_", layer, sep = "")
      layer_path_name = paste(full_name, get(paste(layer, "_path", sep = "")), sep = "")
      layer_path = paste(MTBS_path, layer_path_name, ".shp", sep = "")
      imported_layer = vect(layer_path)
      tmp.burn.bnd = assign(layer_name, imported_layer, envir = .GlobalEnv)
    }

  #############################################################
  #Select study area
  #SA = project(SF_perim, crs(all_rast)) # HERE YOU SET THE SCHOOL FIRE PERIMS
  SA = project(tmp.burn.bnd, crs(all_rast)) # HERE YOU SET THE SCHOOL FIRE PERIMS
  #SA = project(T_perim, crs(all_rast)) # HERE YOU SET THE SCHOOL FIRE PERIMS
  
  sev_rast = all_rast[[c('sev_NBR_GS')]]

  ## crop raster stack of predictors to study area
  sev_SA = mask(sev_rast, SA)
  small_SA = trim(sev_SA, value=NA)
  #plot(small_SA)
  #sev_mask = terra::ifel(small_SA <= -0.02, small_SA, NA)
  sev_mask = terra::ifel(small_SA <= -0.015, small_SA, NA)

  #plot(sev_mask)
  new_SA = ext(sev_mask)
  index_na = which(is.na(sev_mask[]) == TRUE) 

  rast_SA = crop(all_rast, new_SA)
  rast_SA = mask(rast_SA, sev_mask)

  #create data set
  pred.df <- terra::as.data.frame(c(rast_SA), cell=TRUE)
  pred.df = round(pred.df,3);head(pred.df)

  #only complete cases
  pred.df = pred.df %>% filter(complete.cases(.))
  pred.df = pred.df[,2:19]

  #plot prediction
  pred_vals <- terra::predict(rast_SA, final_fitted, na.rm = TRUE)
  pred_vals[index_na] = NA
  #plot(pred_vals)
  #plot(pred_vals$.pred)
  #lines(SA)

 
  ################################################
  #- Plot Map
  ################################################
  x.breaks = c(-117.75,-117.5)
  y.breaks = c(46.25,46.35)
  width_map = 10.2
  location.sbar = "tl"
  
  if (fire_abbrv[i] == "BR") {  width_map = 7.2 }
  if (fire_abbrv[i] == "C")  {width_map=10.2;x.breaks=c(-118.84,-118.80);y.breaks=c(44.0,44.04);location.sbar="tl"}
  if (fire_abbrv[i] == "BC") {width_map=10.2;x.breaks=c(-119.0,-118.95);y.breaks=c(45.02,45.06);location.sbar="tl"}
  if (fire_abbrv[i] == "HR") {width_map=10.2;x.breaks=c(-120.55,-120.45);y.breaks=c(44.45,44.50);location.sbar="br"}
  if (fire_abbrv[i] == "SF") {width_map=10.2;x.breaks=c(-117.75,-117.5);y.breaks=c(46.25,46.35);location.sbar="br"}
  if (fire_abbrv[i] == "F" ) {width_map=10.2;x.breaks=c(-119.20,-119.15);y.breaks=c(44.15,44.20);location.sbar="tl"}
  if (fire_abbrv[i] == "EG") {width_map=10.2;x.breaks=c(-119.6,-119.30);y.breaks=c(43.6,43.8);location.sbar="bl"}
  if (fire_abbrv[i] == "EA") {width_map=10.2;x.breaks=c(-118.52,-118.48);y.breaks=c(44.46,44.52);location.sbar="tl"}
  if (fire_abbrv[i] == "WP") {width_map=10.2;x.breaks=c(-119.9,-119.8);y.breaks=c(44.9,44.95);location.sbar="tl"}
  if (fire_abbrv[i] == "T")  {width_map=10.2;x.breaks=c(-118.80,-118.9);y.breaks=c(44.03,45.05);location.sbar="tl"}
  
  fullrange = c(0,70)
  if (var == 5) {title.leg = "Grass & Herb \n cover (%)" #text.y 
    } else { title.leg =text.y }
  p2 = ggplot() + geom_spatraster(data = pred_vals) +
    geom_sf(SA, mapping = aes(), color = "black", fill = NA,lwd=1.2) + 
    theme(legend.position="none",  axis.text = element_text(size = 30)) +
    scale_y_continuous(breaks = y.breaks) +
    scale_x_continuous(breaks = x.breaks) +
    scale_fill_scico(palette = "bamako", direction = -1, na.value = "transparent", limits = fullrange) +
    #scale_fill_gradientn(
    #  colours = colorRampPalette(c('#D55E00', '#E69F00', '#56B4E9', '#009E73'))(100),
    #  na.value = "transparent",
    #  limits = fullrange) + 
    #scale_fill_gradientn(colours = colorRampPalette(c('red', 'orange', 'lightgreen','darkgreen'))(100),na.value = "transparent",limits = fullrange) +
    #------------------------------
    #- USE THIS FOR THE LEGEND!!!!
    #theme(legend.position = c(0.875, 0.15)) + labs(fill=title.leg) +
    #------------------------------
  ggspatial::annotation_scale(width_hint = 0.2,location = location.sbar,text_cex =3,pad_y = unit(0.5, "cm"),pad_x = unit(0.5, "cm")) + 
    ggspatial::annotation_north_arrow(which_north = "true",location=c("tr"),height = unit(1.8, "cm"), width = unit(1.4, "cm"), 
                    pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"), style = north_arrow_orienteering(text_size = 15)) 
  p2
  
  #ggsave(p2,filename=paste(fig_dir,"pred_map_",fire_abbrv[i],"_",out_text,"_s",seed,".png", sep=''),width = dim(pred_vals)[1]/35, height = dim(pred_vals)[2]/35, dpi = 300)
  ggsave(p2,filename=paste(fig_dir,"pred_map_",fire_abbrv[i],"_",out_text,"_s",seed,".png", sep=''),
         width = width_map, height =   (10.2*dim(pred_vals)[1])/dim(pred_vals)[2], dpi = 300)
 
  # write out predictions
  if (write_out) {
    writeRaster(pred_vals, paste(data_dir,"pred_data_",fire_abbrv[i],"_",out_text,"_s",seed,".tif", sep=''), overwrite=TRUE)
  }
  
  # Calculate final metrics
    raster_range <- global(pred_vals, "range", na.rm = TRUE)
    raster_mean  <- global(pred_vals, "mean", na.rm = TRUE)
    summary_df <- data.frame(
      Metric = c("Min Value", "Max Value", "Average", "Missing (% NA)", "Valid (% Non-NA)"),
      Value  = c(
      raster_range$min,
      raster_range$max,
      raster_mean$mean,
      ((global(is.na(pred_vals), "sum")[, 1])) / ncell(pred_vals) * 100,
      ((ncell(pred_vals) - global(is.na(pred_vals), "sum")[, 1])) / ncell(pred_vals) * 100
    )
  )
  print("#-----------------------------------#")
  print(summary_df)
  print("#-----------------------------------#")
  
  message(paste("Fire done...:",fire_abbrv[i]))
  #if (fire_abbrv[i] == "C") {break}
} # END Running all fires


################################################
#- END
################################################
