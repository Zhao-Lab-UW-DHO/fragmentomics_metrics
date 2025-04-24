library(tidyverse)
library(rsample)
library(tidymodels)
library(parsnip)
library(glmnet)
library(fmsb)
library(pROC)
library(doSNOW)
library(ggplot2)
library(dplyr)
library(themis)
library(patchwork)

cohort <- "GRAIL"

feature_tables <- c(
  "ATAC_entropy",
  "depth",
  "frag_bins",
  "full_gene_depth",
  "mds",
  "se",
  "small_frags",
  "TFBS_entropy"
)

downsample_levels <- c("100M", "50M", "25M", "10M", "5M", "1M")
# downsample_levels <- c("100M")

rep <- "rep1"

starttime <- Sys.time()
GRAIL_comparison_model_results <- tibble()
GRAIL_comparison_ROCs <- tibble()
for (ft in feature_tables) {
  for (dl in downsample_levels) {
    # Set variable for file names
    out_text <- str_c(ft, "_", dl, "_", rep)
    
    # path to feature table
    feature_table_path <- str_c("output/downsample/feature_tables/grail_", ft, "_", dl, ".rds")
    
    # load in data with helper function
    print("Loading Helper Functions and Data")
    
    source("scripts/machine_learning/helper_load_downsample_grail.R")
    
    # set variables
    usegrail <- T
    # trainratio <- 0.7
    # seed <- 1
    
    if(usegrail) {
      df2use <- grail_df
      colors2use <- colorsGRAIL
    } else {
      df2use <- internale1df
      colors2use <- colorsUW
    }
    
    print("All samples table")
    print(table(df2use$phenotype))
    
    
    if (ncol(df2use) > 15000) {
      feat_variance <- sapply(df2use %>% dplyr::select(-c(sample, phenotype, ctdna_fraction)), var) %>% 
        as.data.frame() %>% 
        as_tibble(rownames = "feature") %>% 
        dplyr::rename(variance = 2) %>% 
        arrange(desc(variance)) %>% 
        head(n = 15000) %>% 
        pull(feature)
      
      df2use <- df2use %>% dplyr::select(phenotype, sample, ctdna_fraction, all_of(feat_variance))
    }
    
    print("Feature Table Dimensions:")
    print(dim(df2use))
    
    # TRAINING ----------------------------------------------------------------
    
    for (seed in seq(1,25)) {
      set.seed(seed)
      folds <- vfold_cv(df2use, v=nfold, strata=phenotype, repeats=1)
      
      print('looping through folds to train models')
      print(Sys.time())
      print(out_text)
      print(str_c("Seed: ", seed))
      grail_training_preds <- NULL
      model_coefs <- tibble()
      for(i in 1:nfold) {
        print(i)
        split <- folds$splits[[i]]
        traindatafold <- split %>% analysis() %>% ungroup()
        traindatafold <- traindatafold %>% dplyr::select(-c(ctdna_fraction, sample))
        fit <- trainmodel(traindatafold)
        
        validdatafold <- split %>% assessment()
        pred <- predictmodel(fit,validdatafold)
        grail_training_preds <- bind_rows(grail_training_preds,pred)
        
        split_coefs <- extract_fit_parsnip(fit) %>% broom::tidy() %>%
          mutate(fold = i)
        model_coefs <- bind_rows(model_coefs, split_coefs)
      }
      
      # Data Collecting ---------------------------------------------------------
      
      temp_model_output <- grail_training_preds %>% mutate(feature_table = ft, 
                                                           downsample_level = dl, 
                                                           replicate_number = rep,
                                                           seedused = seed)
      GRAIL_comparison_model_results <- bind_rows(GRAIL_comparison_model_results, temp_model_output)
      
      temp_ROC_results <- extract_ROC_grail_downsample(predictions = grail_training_preds,
                                                       downsample_level = dl, 
                                                       replicate_number = rep, 
                                                       feature = ft, 
                                                       seedused = seed)
      GRAIL_comparison_ROCs <- bind_rows(GRAIL_comparison_ROCs, temp_ROC_results)
      
      # gets timestamp and appends to file output
      file_timestamp <- str_replace_all(string = Sys.time(), pattern = "[-:. ]", replacement = "_") %>% str_sub(1, 19)
      model_output_filename <- str_c("GRAIL_comparison_model_results_downsample_", rep, "_", file_timestamp, ".rds")
      roc_output_filename <- str_c("GRAIL_comparison_ROCs_downsample_", rep, "_", file_timestamp, ".rds")
      
      # saves after every iteration, in case something goes wrong down the line and needs to restart
      if (seed == 25 & dl == "1M") {
        saveRDS(GRAIL_comparison_model_results, str_c("output/machine_learning/", model_output_filename))
        saveRDS(GRAIL_comparison_ROCs, str_c("output/machine_learning/", roc_output_filename))
      }
    }
  }
}

elapsedtime <- Sys.time() - starttime
print(elapsedtime)
