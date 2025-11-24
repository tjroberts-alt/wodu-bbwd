###########################################
# Link bird trend population objectives to habitat objectives
# Starting with JV-specific objectives from PIF step-down

###########################################

library(ebirdst) #if using the 2022 version, don't run this here, but after 'remotes::install_version' line
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(terra)
library(arrow) # for read_parquet() function
library(here)
library(mgcv)
library(tidyverse)
library(stringr)

path <- file.path("C:/NAWMP_analysis/NSST_stepdown/ebirdTrends")

### Code to link trends with habitat variables can be found here:
## https://github.com/CLO-Conservation-Science/Habitat-Correlates-Trends-Case-Study/tree/main/Code

# As of the 2023 release, ebird is using a different raster template from 2022. This results in the confusing situation that for the latest version of ebirdst the Status results (version 2023) and Trends results (version 2022) use different templates with different projections. Given this, install the 2022 version of the R package to access the 2022 version of both Status and Trends. This is also ideal as it allows to work with the Status version that the 2012-2022 Trends is based on.
# install the 2022 version of ebirdst to access Status results matching Trends
remotes::install_version("ebirdst", version = "3.2022.3")
library(ebirdst)

# eBird data
eBird_sp_list <- ebirdst_runs  ## identify species codes#

# grid info
srd_covs_all <- read_parquet("H:/NAWMP/NSST/StepDown/srd2022_2012-2022_north-america/srd2022_2012-2022_north-america.parquet")
# just choose the elevation variable we want to use here and rename. Makes it easier later on
srd_covs_all$elev <- srd_covs_all$elevation_30m_median 

# Get individual year data
srd12 <- srd_covs_all[srd_covs_all$year == 2012,]       
srd22 <- srd_covs_all[srd_covs_all$year == 2022,]

#### Species list: whistling ducks and wood ducks
spp_list <- c("bbwduc", "wooduc")


## Multiple species raster prep (This limits the analysis to cells where both spp are able to be modeled. Should consider changing this to either all wodu cells or all wodu cells where bbwd have abundance estimates, ignoring trend)
rasters <- setNames(
  lapply(spp_list, function(sp) 
    load_raster(sp, product = "abundance", period = "weekly",
                metric = "median", resolution = "27km")),
  spp_list
)

# Step 1: collect common_ids for all species
all_ids <- map(spp_list, ~ {
  sp_trend <- load_trends(.x, fold_estimates = TRUE)
  sp_trend$srd_id
})

# Step 2: add in srd12 and srd22 ids
all_ids <- c(all_ids, list(srd12$srd_id, srd22$srd_id))

# Step 3: take intersection
common_ids <- reduce(all_ids, intersect)

# After computing common_ids
srd12 <- srd12 %>% filter(srd_id %in% common_ids) %>% arrange(srd_id)
srd22 <- srd22 %>% filter(srd_id %in% common_ids) %>% arrange(srd_id)

# function to select appropriate cells
get_abd_trend <- function(sp, common_ids) {
  sp_trend <- load_trends(sp, fold_estimates = TRUE)
  
  sp_trend_f <- sp_trend %>%
    filter(srd_id %in% common_ids) %>%
    arrange(srd_id)
  
  trend_wide <- sp_trend_f %>%
    pivot_wider(
      id_cols = c(species_code, fold),
      names_from = srd_id,
      values_from = abd_ppy
    ) %>%
    arrange(fold)
  
  abd_wide <- sp_trend_f %>%
    pivot_wider(
      id_cols = c(species_code, fold),
      names_from = srd_id,
      values_from = abd
    ) %>%
    arrange(fold) %>%
    slice(1) %>%
    select(-fold)
  
  list(abd = abd_wide, trend = trend_wide)
}

# Step 4: run for all species
results <- map(spp_list, get_abd_trend, common_ids = common_ids)

# Step 5: combine
abd_all <- bind_rows(lapply(results, `[[`, "abd"))
abd_all <- abd_all %>%
  column_to_rownames("species_code")

trend_reps_all <- bind_rows(lapply(results, `[[`, "trend"))


# The code below is limited to landscape variables used by eBird team. Consider getting Forest Inventory Analysis data at the scale of interest, or other variables specific to the region of interest

# Before starting the loop below, decide on which habitat/landscape variables to use 
# then can filter srd12, srd22, and srd_agg_pred to those so when using pred_names and creating 'D' below
# Available names are (can create our own variables, but not right now): 
#names(srd_covs_all)

# for wood duck/whistling duck
vars <- c("srd_id","year","latitude","longitude","elev","road_density_c1","road_density_c2","road_density_c3","road_density_c4","road_density_c5","astwbd_c1_ed","astwbd_c1_pland","astwbd_c2_ed","astwbd_c2_pland","astwbd_c3_ed","astwbd_c3_pland","mcd12q1_lccs1_c13_pland","mcd12q1_lccs1_c14_pland","mcd12q1_lccs1_c15_pland","mcd12q1_lccs1_c16_pland","mcd12q1_lccs1_c21_pland","mcd12q1_lccs1_c22_pland","mcd12q1_lccs1_c255_pland","mcd12q1_lccs1_c31_pland","mcd12q1_lccs1_c32_pland","mcd12q1_lccs1_c41_pland","mcd12q1_lccs1_c42_pland","mcd12q1_lccs1_c43_pland","mcd12q1_lccs2_c9_pland","mcd12q1_lccs2_c25_pland","mcd12q1_lccs2_c35_pland",            "mcd12q1_lccs2_c36_pland","mcd12q1_lccs3_c27_pland","mcd12q1_lccs3_c50_pland","mcd12q1_lccs1_c13_ed","mcd12q1_lccs1_c14_ed","mcd12q1_lccs1_c15_ed","mcd12q1_lccs1_c16_ed","mcd12q1_lccs1_c21_ed","mcd12q1_lccs1_c22_ed","mcd12q1_lccs1_c255_ed","mcd12q1_lccs1_c31_ed","mcd12q1_lccs1_c32_ed","mcd12q1_lccs1_c41_ed","mcd12q1_lccs1_c42_ed","mcd12q1_lccs1_c43_ed","mcd12q1_lccs2_c9_ed","mcd12q1_lccs2_c25_ed","mcd12q1_lccs2_c35_ed","mcd12q1_lccs2_c36_ed","mcd12q1_lccs3_c27_ed","mcd12q1_lccs3_c50_ed")

srd12 <- srd12 %>% select(all_of(vars))
srd22 <- srd22 %>% select(all_of(vars))

# For wood duck/whistling duck add both spp abd and trends to srds
# add abundance
id <- "bbwduc"
srd12[[id]] <- abd_all[id, as.character(srd12$srd_id)] %>% as.numeric()
srd22[[id]] <- abd_all[id, as.character(srd22$srd_id)] %>% as.numeric()
id <- "wooduc"
srd12[[id]] <- abd_all[id, as.character(srd12$srd_id)] %>% as.numeric()
srd22[[id]] <- abd_all[id, as.character(srd22$srd_id)] %>% as.numeric()

# add mean trend
bbwd <- load_trends('bbwduc')
bbwd_f <- bbwd %>%
  filter(srd_id %in% common_ids) %>%
  arrange(srd_id)
wodu <- load_trends('wooduc')
wodu_f <- wodu %>%
  filter(srd_id %in% common_ids) %>%
  arrange(srd_id)


srd12 <- srd12 %>%
  left_join(bbwd_f %>% select(srd_id, abd_trend), by = "srd_id") %>%
  rename(bbwd_trend = abd_trend)
srd12 <- srd12 %>%
  left_join(wodu_f %>% select(srd_id, abd_trend), by = "srd_id") %>%
  rename(wodu_trend = abd_trend)

srd22 <- srd22 %>%
  left_join(bbwd_f %>% select(srd_id, abd_trend), by = "srd_id") %>%
  rename(bbwd_trend = abd_trend)
srd22 <- srd22 %>%
  left_join(wodu_f %>% select(srd_id, abd_trend), by = "srd_id") %>%
  rename(wodu_trend = abd_trend)


#### Prep land cover data for modeling and predictions
srd12 <- rename(srd12, ROW_NUM = srd_id)        
srd22 <- rename(srd22, ROW_NUM = srd_id)
srd_agg_pred <- srd12 %>% select(-ROW_NUM)


#### Create folders to store results for each species
for (i in seq_along(spp_list)) {
  dir.create(file.path(path, "Results", spp_list[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "Results", spp_list[i], "GAM_outputs"), recursive = TRUE, showWarnings = FALSE)
}


#### ------------------------- 2. Giant loop over species and folds ----------------------- ####
tic()  # Start timer. Could take up to 1.5 hours on a laptop.

for (iii_spp in 1:length(spp_list)){
  
  #### Organize species-level data
  spp_name <- spp_list[iii_spp] 
  spp_index <- rownames(abd_all) == spp_name
  abd_spp <- abd_all[spp_index, ]
  trend_reps_spp <- trend_reps_all %>% filter(species_code==spp_name) %>%
    select(-c(species_code, fold))
  
  #### Identify Predictors that vary sufficiently for smooths  
  min_ss <- 30
  srd_spp <- srd_agg_pred %>% mutate(abd = as.numeric(abd_spp[1, ])) %>%
    drop_na(abd) %>%
    select(-abd)
  pred_names <- names(srd_spp)
  
  # Finding number of unique values in each column
  pred_unique <- numeric(length(pred_names))
  
  for (iii in seq_along(pred_names)) {
    pred_unique[iii] <- length(unique(srd_agg_pred[[pred_names[iii]]]))
  }
  ##
  
  pred_names <- pred_names[ pred_unique > min_ss ] # min_ss = minimum sample size
  pred_names <-  pred_names[ !is.na(pred_names )]
  pred_names
  print(paste("starting loop through folds for", spp_name))
  
  
  # =============================================================
  # Start loop through 100 folds of trend estimates 
  # =============================================================
  
  #### Loop through 100 folds to propagate trend uncertainty
  r.sq <- NULL
  for(iii_fold in 1:nrow(trend_reps_spp)){
    print(paste("fold",iii_fold, spp_name)) 
    
    #### Prepare the data for modeling   
    resp <- data.frame(
      response = as.numeric(trend_reps_spp[iii_fold, ]),
      ROW_NUM = as.numeric(colnames(trend_reps_spp)), 
      abd = as.numeric(abd_spp[1, ]) )
    D <- cbind(resp, srd_agg_pred)
    
    ## Weight the likelihood by relative abundance (normalized)
    D$weights <- D$abd/mean(D$abd, na.rm = T)
    D <- D[!is.na(D$response),]
    
    ## Set up objects to receive outputs
    if(iii_fold == 1){
      gam.summaries <- as.data.frame(matrix(nrow = 0, ncol = length(pred_names)+1))
      names(gam.summaries) <- c("Fold", "Metric", pred_names[4:length(pred_names)], "Elev", "Lat-Lon")
    }
    
    
    
    # =============================================================
    # Run the GAM
    # =============================================================
    
    #### Construct model Formula
    model_name <- NULL
    for (iii in 4:length(pred_names)){
      model_name[iii-3] <- 
        paste0( "s( ",pred_names[iii],",bs=\"ds\", k=4, m=c(1,0)) + ")
    }      
    model_name <- c("response ~", model_name,
                    "s(elev, bs=\"ds\", k=10, m=c(1,0)) + ",
                    "s(longitude, latitude, bs=\"ds\",k=200, m=c(1,0.5)) +", 
                    "1")
    model_name <- paste(model_name, collapse = " ")
    
    #### Run the model. Note Selection + Extra Gamma penalization
    d.gam <- bam( 
      as.formula(model_name), 
      weights = weights, 
      gamma = 5,
      select = T, 
      discrete = T, 
      data = D )
    #summary(d.gam)
    
    
    # =============================================================
    # Save Model Summaries  
    # =============================================================
    
    save(d.gam, file = file.path(path, "Results", spp_name, "GAM_outputs",
                                 paste0(spp_name, "_gam_fold", iii_fold, ".RData"))
    ) # save model object
    
    #### Save F value, p value, and R2 from each fold
    d.sum <- summary(d.gam)
    d.add <- as.data.frame(d.sum$s.table) %>% select(-c(edf, Ref.df)) %>%
      rownames_to_column("predictor") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("Metric") %>%
      mutate(Fold = rep(iii_fold), .before = Metric) %>%
      slice(-1)
    names(d.add) <- names(gam.summaries)
    gam.summaries <- rbind(gam.summaries, d.add)
    
    r.sq[iii_fold] <- d.sum$r.sq
    
    
    #### Predictor Effect Estimates (type="iterms")
    pred  <- predict.gam(
      d.gam, 
      type = "iterms", se.fit = TRUE) 
    
    ## Save out predictor effect estimates for 1-D marginal effect plots
    d.effect_plotting <- cbind(data.frame(ROW_NUM = D$ROW_NUM), data.frame(pred)) # marginal effects and SE for this fold
    d.effect_plotting <- d.effect_plotting %>%
      mutate(Fold = iii_fold, .before = ROW_NUM)
    
    if(iii_fold == 1){
      effect_plotting <- d.effect_plotting
      rm(d.effect_plotting)
    } else{
      effect_plotting <- rbind(effect_plotting, d.effect_plotting)
    }
    
    
    
    # =============================================================
    # Predict to 2022 land cover surface
    # =============================================================
    
    #### Generate predictions and export landcover data formatted for ROI
    ## Subset srd22 to the same rows used to fit the GAM
    #lc22 <- srd22 %>% select(ROW_NUM, Barren.PLAND:Herbaceous.Wetlands.PLAND) # original. Need to match to names of srd columns used here. Removing year, lat/long, and elev
    lc22 <- srd22 %>% select(c(1, 6:ncol(srd22)))
    D22 <- D %>% select(response:elev, weights) %>% left_join(lc22, by = "ROW_NUM")
    D22 <- D22[names(D)]
    
    ## Predict the GAM to the 2022 land covers
    pred22 <- predict.gam(d.gam, newdata = D22, type = "iterms", se.fit = TRUE) 
    
    ## Export 2012 and 2022 landcover datasets used in the model
    if(iii_fold == 1){
      D22_output <- D22 %>% select(-c(response,weights))
      D_output <- D %>% select(-c(response,weights))
      write.csv(D22_output, file.path(path, "Results", spp_name, paste0(spp_name, "_2022_landcovers.csv")), row.names = FALSE)
      write.csv(D_output, file.path(path, "Results", spp_name, paste0(spp_name, "_2012_landcovers.csv")), row.names = FALSE)
    } 
    
    
    # ===================================================================================
    # Produce datasets for interpretive plotting aids: Geographic effects
    # ===================================================================================
    
    #### Create containers to hold plotting results
    save.geo <- cbind("longitude" = D22$longitude, "latitude" = D22$latitude, "ROW_NUM" = D22$ROW_NUM)
    save.geo.abd <- cbind("longitude" = D22$longitude, "latitude" = D22$latitude, "ROW_NUM" = D22$ROW_NUM)
    
    #### Initiate loop through all predictors
    for(iii_effect in 1:length(colnames(pred22$fit))){
      zzz_name <- names(d.gam$var.summary)[iii_effect]
      
      if (zzz_name != "longitude" & zzz_name != "elev"){
        ## Geographic PPY effect dataset
        zzz <- pred22$fit[, iii_effect]
        save.geo <- cbind(save.geo, zzz)
        colnames(save.geo)[iii_effect+3] <- zzz_name
        
        ## Geographic abundance-weighted PPY effect dataset
        # Zero out areas w/o feature of interest so that importance is not summarized there
        zero_feature <- D22[, names(D22) == zzz_name ] 
        zero_feature[zero_feature!=0] <- 1
        zzz <- pred22$fit[,iii_effect] * D22$abd * zero_feature
        
        zzz_limit <- max(abs(quantile(zzz, probs = c(0.01, 0.99))))    # Plotting aid for color gradients
        zzz[zzz < -1*zzz_limit] <- -1*zzz_limit
        zzz[zzz > zzz_limit] <- zzz_limit
        save.geo.abd <- cbind(save.geo.abd, zzz)
        colnames(save.geo.abd)[iii_effect+3] <- zzz_name
      } 
    } # Close loop through Predictors
    
    
    #### Bind the results from this fold to the main dataset
    save.geo <- as.data.frame(save.geo)
    save.geo.abd <- as.data.frame(save.geo.abd)
    save.geo$fold <- rep(iii_fold, length.out=nrow(save.geo))
    save.geo.abd$fold <- rep(iii_fold, length.out=nrow(save.geo.abd))
    
    if(iii_fold == 1){
      save.geo.folds <- save.geo
      save.geo.abd.folds <- save.geo.abd
    } else{
      save.geo.folds <- rbind(save.geo.folds, save.geo)
      save.geo.abd.folds <- rbind(save.geo.abd.folds, save.geo.abd)
    }
    
    
  } # Close loop through 100 folds
  
  
  #### Export datasets to directory
  write.csv(save.geo.folds, file.path(path, "Results", spp_name, 
                                      paste0(spp_name,"_geo.effect_plotting.csv")), row.names = F)
  write.csv(save.geo.abd.folds, file.path(path, "Results", spp_name, 
                                          paste0(spp_name,"_geo.abd.effect_plotting.csv")), row.names = F)
  write.csv(gam.summaries, file.path(path, "Results", spp_name, 
                                     paste0(spp_name,"_GAM_summaries.csv")), row.names = F)
  write.csv(as.data.frame(r.sq), file.path(path, "Results", spp_name, 
                                           paste0(spp_name,"_R2.csv")), row.names = F)
  write.csv(effect_plotting, file.path(path, "Results", spp_name, 
                                       paste0(spp_name,"_1Deffect_plotting.csv")), row.names = F)
  rm(save.geo.folds, save.geo.abd.folds, save.geo, save.geo.abd, r.sq)
  
} # Close loop through Sagebrush species


#toc()   # end timer


#### ------------------------------------------------------------------------------ ####
#### This script produces plots to visualize results 
##   1. Create maps of species abundance
##   2. Create maps of species trend
##   3. Create maps of % land cover, geographic effect, & effect on abundance for each predictor
##   4. Create 1D marginal effect plot for each predictor
#### ------------------------------------------------------------------------------ ####

spp_names <- c("Black-bellied Whistling-Duck", "Wood Duck")

#### Create Figure directory for each species
for (i in seq_along(spp_list)) {
  dir.create(file.path(path, "Results", "Figures", "Maps", spp_list[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "Results", "Figures", "Marginal_effects_plots", spp_list[i], "GAM_outputs"), recursive = TRUE, showWarnings = FALSE)
}

for(i in 1:length(spp_list)){
  
  #### ------------------------- 1. Create map of species abundance ------------------------ ####
  
  #### Load data
  abd <- read.csv(file.path(path, "Results", spp_list[i], paste0(spp_list[i], "_2022_landcovers.csv")))
  abd <- abd %>% select(abd:longitude)
  
  # Calculate bounding box with 5% buffer for the species to use in all plots
  bbox <- abd %>%
    summarise(
      x_min = min(longitude, na.rm = TRUE),
      x_max = max(longitude, na.rm = TRUE),
      y_min = min(latitude, na.rm = TRUE),
      y_max = max(latitude, na.rm = TRUE)
    ) %>%
    mutate(
      x_range = x_max - x_min,
      y_range = y_max - y_min,
      x_min = x_min - 0.05 * x_range,
      x_max = x_max + 0.05 * x_range,
      y_min = y_min - 0.05 * y_range,
      y_max = y_max + 0.05 * y_range
    )
  
  
  
  #### Create the plot
  ggplot(data = abd) + 
    borders("world", col = "grey50", fill = "grey80", size = 0.075) +
    borders("state", col = "grey50", fill = NA, size = 0.075) +
    geom_point(aes(x = longitude, y = latitude, col = abd, fill = abd), pch = 21) + 
    scale_fill_gradient2(low = "white", mid = "#41b6c4", high = "#023858",
                         midpoint = max(abd$abd)/2, 
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"), 
                         aesthetics = c("color", "fill"),
                         name = "Relative \nabundance") +
    coord_cartesian(xlim = c(bbox$x_min, bbox$x_max), ylim = c(bbox$y_min, bbox$y_max), clip = "on") +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(color = NA, fill = NA), #transparent legend bg
          legend.position = c(0.11, 0.15),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.3, "cm"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = 18)) +
    xlab("Longitude") + ylab("Latitude")
  
  ggsave(file.path(path, "Results", "Figures","Maps",spp_list[i],paste0(spp_list[i],"_Abundance.png")), dpi = 600, width = 5, height = 6, units = "in")
  
  
  
  #### ------------------------- 2. Create map of species trend ------------------------ ####
  
  #### Organize trend data for species i
  spp_trend <- trend_reps_all %>% filter(species_code==spp_list[i]) %>%
    #rename_with(~gsub("X", "", .x)) %>%
    select(-c(species_code, fold)) %>% 
    t() %>%
    as.data.frame() %>% 
    drop_na()
  
  ## Calculate median, upper, and lower 80% trend estimate 
  spp_trend_plot <- data.frame(latitude = abd$latitude, 
                               longitude = abd$longitude,
                               ppy_median = apply(spp_trend, 1, median),
                               ppy_upperCI80 = apply(spp_trend, 1, quantile, probs = 0.9),
                               ppy_lowerCI80 = apply(spp_trend, 1, quantile, probs = 0.1)) 
  
  ## To aid visualization only, max out the trends at +4 and -4
  spp_trend_plot$ppy_median[which(spp_trend_plot$ppy_median > 4)] <- 4
  spp_trend_plot$ppy_median[which(spp_trend_plot$ppy_median < -4)] <- -4
  
  
  #### Create the plot (change "ppy_median" to "pyy_upperCI80" or "ppy_lowerCI80" to plot trend confidence intervals)
  ggplot(data = spp_trend_plot) +
    borders("world", col = "grey50", fill = "grey80", size = 0.075) +
    borders("state", col = "grey50", fill = NA, size = 0.075) +
    geom_point(aes(x = longitude, y = latitude, col = ppy_median, fill = ppy_median, size = abd$abd), pch = 21) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#023858",
                         midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                         aesthetics = c("color", "fill"),
                         name = "PPY \nTrend") +
    scale_radius(range = c(0.0, 1.5), guide = "none") +
    #coord_cartesian(xlim = c(-125, -103), ylim = c(33,52), clip = "on") +
    coord_cartesian(xlim = c(bbox$x_min, bbox$x_max), ylim = c(bbox$y_min, bbox$y_max), clip = "on") +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(color = NA, fill = NA), #transparent legend bg
          legend.position = c(0.1, 0.15),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.3, "cm"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = 18))
  
  ggsave(file.path(path, "Results", "Figures","Maps",spp_list[i],paste0(spp_list[i],"_Trend.png")), dpi = 600, width = 5, height = 6, units = "in")
  
  
  #### --------------- 3. Create map of land cover, geographic effect, & effect on abundance for each predictor -------------- ####
  
  #### Load data
  pland22 <- read.csv(file.path(path, "Results", spp_list[i], paste0(spp_list[i], "_2022_landcovers.csv")))
  pland22 <- pland22 %>% select(-c(ROW_NUM:abd, elev))
  
  pland12 <- read.csv(file.path(path, "Results", spp_list[i], paste0(spp_list[i], "_2012_landcovers.csv")))
  pland12 <- pland12 %>% select(-c(ROW_NUM:abd, elev))
  
  geo <- read.csv(file.path(path, "Results", spp_list[i], paste0(spp_list[i], "_geo.effect_plotting.csv")))
  geo.abd <- read.csv(file.path(path, "Results", spp_list[i], paste0(spp_list[i], "_geo.abd.effect_plotting.csv")))
  
  effect_1D <- read.csv(file.path(path, "Results", spp_list[i], paste0(spp_list[i],"_1Deffect_plotting.csv")))
  
  
  #### Create plotting datasets with mean and sd of effects for each PLAND (% land cover variable)
  
  ## For the geographic trend effect dataset
  geo_plot_mean <- geo %>% group_by(ROW_NUM) %>%
    summarize(across(everything(), mean)) %>%
    select(-c(fold, ROW_NUM))
  geo_plot_sd <- geo %>% group_by(ROW_NUM) %>%
    select(-c(latitude, longitude, fold)) %>%
    summarize(across(everything(), sd)) %>%
    mutate(latitude = geo_plot_mean$latitude,
           longitude = geo_plot_mean$longitude, .before = ) %>%
    select(latitude, longitude, everything(), -ROW_NUM) 
  
  ## For the abundance-weighted effect dataset
  geo.abd_plot_mean <- geo.abd %>% group_by(ROW_NUM) %>%
    summarize(across(everything(), mean)) %>%
    select(-c(fold, ROW_NUM))
  geo.abd_plot_sd <- geo.abd %>% group_by(ROW_NUM) %>%
    select(-c(latitude, longitude, fold)) %>%
    summarize(across(everything(), sd)) %>%
    mutate(latitude = geo_plot_mean$latitude,
           longitude = geo_plot_mean$longitude) %>%
    select(latitude, longitude, everything(), -ROW_NUM) 
  rm(geo, geo.abd)
  
  ## For the 1-D marginal effects dataset
  effect_1D_mean <- effect_1D %>% group_by(ROW_NUM) %>%
    summarize(across(everything(), mean)) %>%
    select(-c(Fold, ROW_NUM))
  # This is equivalent to taking the mean of the mean, upper and lower SE bands for each fold. 
  rm(effect_1D)
  
  #### Start forloop
  for(l in 3:ncol(geo_plot_mean)){
    pland.names <- str_remove(colnames(geo_plot_mean), ".PLAND") %>%
      str_replace_all("[.]"," ") %>%
      tolower()
    
    #### Create the map for PLAND
    pland22 <- pland22[names(geo_plot_mean)]  # make sure column orders match
    pland_sub <- pland22[,c(1:2,l)]
    
    ggplot(data = pland_sub) +
      borders("world", col = "grey50", fill = "grey80", size = 0.075) +
      borders("state", col = "grey50", fill = NA, size = 0.075) +
      geom_point(aes(x = longitude, y = latitude, col = pland_sub[,3], fill = pland_sub[,3]), pch = 21) +
      scale_fill_gradient2(low = "white", mid = "#78c679", high = "#00441B",
                           midpoint = 50,
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                           aesthetics = c("color", "fill"),
                           name = "% in \n2021", limits = c(0, 100)) +
      #coord_cartesian(xlim = c(-125, -103), ylim = c(33,52), clip = "on") +
      coord_cartesian(xlim = c(bbox$x_min, bbox$x_max), ylim = c(bbox$y_min, bbox$y_max), clip = "on") +
      theme(panel.border = element_rect(colour = "black", fill = NA),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(color = NA, fill = NA), #transparent legend bg
            legend.position = c(0.12, 0.14),
            legend.key.width = unit(0.35, "cm"),
            legend.key.height = unit(0.4, "cm"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = 18)) +
      xlab("Longitude") + ylab("Latitude")
    
    ggsave(file.path(path, "Results", "Figures","Maps",spp_list[i],paste0(spp_list[i],"_",pland.names[l],"_2022PLAND.png")), dpi = 600, width = 5, height = 6, units = "in")
    
    
    #### Create the map for the mean geographic effect
    geo_sub <- as.data.frame(geo_plot_mean[,c(1:2,l)])
    # when geo_sub[,3] has only one unique value, breaks vector has length 2 (min and max collapse to the same), but labels vector is still length 3. Fix by handling this dynamically constructing the breaks and labels so they always match in length.
    vals <- geo_sub[,3]
    # unique range
    vals_min <- min(vals, na.rm = TRUE)
    vals_max <- max(vals, na.rm = TRUE)
    
    if (vals_min == vals_max) {
      # Only one value, just use that as the single break
      breaks_use <- vals_min
      labels_use <- "no variation"
    } else {
      # Normal case
      breaks_use <- c(vals_min, 0, vals_max)
      labels_use <- c("more \nnegative", "", "more \npositive")
    }
    
    ggplot(data = geo_sub) +
      borders("world", col = "grey50", fill = "grey80", size = 0.075) +
      borders("state", col = "grey50", fill = NA, size = 0.075) +
      geom_point(aes(x = longitude, y = latitude, col = vals, fill = vals), pch = 21) +
      scale_fill_gradient2(low = "#40004B", mid = "#f7f7f7", high = "#00441B",
                           midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = NA),
                           aesthetics = c("color", "fill"),
                           name = "",
                           breaks = breaks_use,
                           labels = labels_use) +
      #coord_cartesian(xlim = c(-125, -103), ylim = c(33,52), clip = "on") +
      coord_cartesian(xlim = c(bbox$x_min, bbox$x_max), ylim = c(bbox$y_min, bbox$y_max), clip = "on") +
      theme(panel.border = element_rect(colour = "black", fill = NA),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(color = NA, fill = NA), #transparent legend bg
            legend.position = c(0.12, 0.14),
            legend.key.width = unit(0.35, "cm"),
            legend.key.height = unit(0.4, "cm"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = 18),
            legend.spacing.x = unit(0.1, "cm")) +
      xlab("Longitude") + ylab("Latitude")
    
    ggsave(file.path(path, "Results", "Figures","Maps",spp_list[i],paste0(spp_list[i],"_",pland.names[l],"_GeographicEffect_MEAN.png")), dpi = 600, width = 5, height = 6, units = "in")
    
    
    #### Create the map for the SD geographic effect
    geo_sub <- as.data.frame(geo_plot_sd[,c(1:2,l)])
    ggplot(data = geo_sub) +
      borders("world", col = "grey50", fill = "grey80", size = 0.075) +
      borders("state", col = "grey50", fill = NA, size = 0.075) +
      geom_point(aes(x = longitude, y = latitude, col = geo_sub[,3], fill = geo_sub[,3]), pch = 21) +
      scale_fill_gradient2(low = "#FFFFFF", mid = "#CBC2F2", high = "#0800B7",
                           midpoint = max(geo_sub[,3])/2, 
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = NA), 
                           aesthetics = c("color", "fill"),
                           name = "SD") +
      #coord_cartesian(xlim = c(-125, -103), ylim = c(33,52), clip = "on") +
      coord_cartesian(xlim = c(bbox$x_min, bbox$x_max), ylim = c(bbox$y_min, bbox$y_max), clip = "on") +
      theme(panel.border = element_rect(colour = "black", fill = NA),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(color = NA, fill = NA), #transparent legend bg
            legend.position = c(0.12, 0.14),
            legend.key.width = unit(0.35, "cm"),
            legend.key.height = unit(0.4, "cm"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = 18),
            legend.spacing.x = unit(0.1, "cm")) +
      xlab("Longitude") + ylab("Latitude")
    
    ggsave(file.path(path, "Results", "Figures","Maps",spp_list[i],paste0(spp_list[i],"_",pland.names[l],"_GeographicEffect_SD.png")), dpi = 600, width = 5, height = 6, units = "in")
    
    
    #### Create the map for the mean effect on change in abundance
    geo.abd_sub <- as.data.frame(geo.abd_plot_mean[,c(1:2,l)])
    # when geo.abd_sub[,3] has only one unique value, breaks vector has length 2 (min and max collapse to the same), but labels vector is still length 3. Fix by handling this dynamically constructing the breaks and labels so they always match in length.
    vals <- geo.abd_sub[,3]
    # unique range
    vals_min <- min(vals, na.rm = TRUE)
    vals_max <- max(vals, na.rm = TRUE)
    
    if (vals_min == vals_max) {
      # Only one value, just use that as the single break
      breaks_use <- vals_min
      labels_use <- "no variation"
    } else {
      # Normal case
      breaks_use <- c(vals_min, 0, vals_max)
      labels_use <- c("more \nnegative", "", "more \npositive")
    }
    
    ggplot(data = geo.abd_sub) +
      borders("world", col = "grey50", fill = "grey80", size = 0.075) +
      borders("state", col = "grey50", fill = NA, size = 0.075) +
      geom_point(aes(x = longitude, y = latitude, col = vals, fill = vals), pch = 21) +
      scale_fill_gradient2(low = "#40004B", mid = "#f7f7f7", high = "#00441B",
                           midpoint = 0,
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = NA),
                           aesthetics = c("color", "fill"),
                           name = "",
                           breaks = breaks_use,
                           labels = labels_use) +
      #coord_cartesian(xlim = c(-125, -103), ylim = c(33,52), clip = "on") +
      coord_cartesian(xlim = c(bbox$x_min, bbox$x_max), ylim = c(bbox$y_min, bbox$y_max), clip = "on") +
      theme(panel.border = element_rect(colour = "black", fill = NA),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(color = NA, fill = NA), #transparent legend bg
            legend.position = c(0.12, 0.14),
            legend.key.width = unit(0.35, "cm"),
            legend.key.height = unit(0.4, "cm"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = 18)) +
      xlab("Longitude") + ylab("Latitude")
    
    ggsave(file.path(path, "Results", "Figures","Maps",spp_list[i],paste0(spp_list[i],"_",pland.names[l],"_ChangeInAbundance_MEAN.png")), dpi = 600, width = 5, height = 6, units = "in")
    
    
    #### Create the map for the sd effect on change in abundance
    geo.abd_sub <- as.data.frame(geo.abd_plot_sd[,c(1:2,l)])
    ggplot(data = geo.abd_sub) +
      borders("world", col = "grey50", fill = "grey80", size = 0.075) +
      borders("state", col = "grey50", fill = NA, size = 0.075) +
      geom_point(aes(x = longitude, y = latitude, col = geo.abd_sub[,3], fill = geo.abd_sub[,3]), pch = 21) +
      scale_fill_gradient2(low = "#FFFFFF", mid = "#CBC2F2", high = "#0800B7",
                           midpoint = max(geo.abd_sub[,3])/2, 
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = NA), 
                           aesthetics = c("color", "fill"),
                           name = "SD") +
      #coord_cartesian(xlim = c(-125, -103), ylim = c(33,52), clip = "on") +
      coord_cartesian(xlim = c(bbox$x_min, bbox$x_max), ylim = c(bbox$y_min, bbox$y_max), clip = "on") +
      theme(panel.border = element_rect(colour = "black", fill = NA),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(color = NA, fill = NA), #transparent legend bg
            legend.position = c(0.12, 0.14),
            legend.key.width = unit(0.35, "cm"),
            legend.key.height = unit(0.4, "cm"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = 18)) +
      xlab("Longitude") + ylab("Latitude")
    
    ggsave(file.path(path, "Results", "Figures","Maps",spp_list[i],paste0(spp_list[i],"_",pland.names[l],"_ChangeInAbundance_SD.png")), dpi = 600, width = 5, height = 6, units = "in")
    
    
    
    #### ------------------------- 4. Create 1D effect plot for each predictor ------------------------ ####
    
    effect_1D_sub <- effect_1D_mean[,grepl(colnames(geo_plot_mean)[l], colnames(effect_1D_mean))]
    effect_1D_sub[, colnames(geo_plot_mean)[l]] <- pland12[, colnames(geo_plot_mean)[l]]
    effect_1D_sub <- as.data.frame(effect_1D_sub)
    
    
    ggplot(data = effect_1D_sub) +
      geom_line(mapping = aes(x = effect_1D_sub[,3], y = effect_1D_sub[,1])) +
      geom_ribbon(aes(x = effect_1D_sub[,3], ymin = effect_1D_sub[,1]-effect_1D_sub[,2], 
                      ymax = effect_1D_sub[,1]+effect_1D_sub[,2]), alpha = 0.5, fill = "grey80")  +
      coord_cartesian(xlim = c(0, max(effect_1D_sub[,3])+5),
                      ylim = c(-2, 2), expand = FALSE) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill= NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text = element_text(size = 18)) +
      xlab(paste0("Percent ", pland.names[l])) + ylab("Marginal effect size")
    
    ggsave(file.path(path, "Results", "Figures","Marginal_effects_plots",spp_list[i],paste0(spp_list[i],"_",pland.names[l],"_1D_effectplot.png")), dpi = 600, width = 5, height = 6, units = "in")
    
    
  } # l - predictors
  
} # i = species
