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

usegrail <- TRUE

feature_tables <- c(
  "SE",
  "depth",
  "frag_bins",
  "full_gene_depth",
  "small_frags",
  "mds",
  "ATAC_entropy",
  "TFBS_entropy"
)

usegrail <- TRUE

grail_metadata_simple_path <- "files/grail_ml_metadata.rds"

# starttime <- Sys.time()
grail_mixed_validation_preds <- tibble()
grail_mixed_validation_ROCS <- tibble()
starttime <- Sys.time()
for (ft in feature_tables) {
  # set large if statements pointing to FTs
  grail_mixing_ft_path <- str_c("output/mixing_data/feature_tables/", ft, ".rds")
  
  # load in data set above
  source("scripts/machine_learning/helper_load_mixing_grail.R")
  
  if (usegrail) {
    df2train <- grail_metadata_simple
    mixed_data <- grail_mixing_ft
  }
  
  # loop through ML
  seeds <- read_rds("scripts/machine_learning/paramseed.rds")[1:25]
  # for (SEED in c(1)) { # for testing
  for (SEED in seeds) {
    print(str_c(Sys.time(), "Feature Table:", ft, "Seed:", SEED, sep = " "))
    set.seed(SEED)
    
    # Get samples that will be used for training/validation
    grail_data_split <- initial_split(df2train, prop = 0.7, strata = phenotype)
    grail_training_samples <- training(grail_data_split) %>% pull(sample)
    grail_validation_samples <- testing(grail_data_split) %>% pull(sample)
    
    # extract samples into training/validation
    mixed_training_df <- mixed_data %>% 
      separate(sample, into = c("cancer", "normal", "cancer_ratio", "normal_ratio"), remove = FALSE, sep = "_") %>% 
      filter(cancer %in% grail_training_samples) %>% 
      filter(normal %in% grail_training_samples) %>% 
      dplyr::select(-c(cancer, normal, cancer_ratio, normal_ratio))
    
    mixed_validation_df <- mixed_data %>% 
      separate(sample, into = c("cancer", "normal", "cancer_ratio", "normal_ratio"), remove = FALSE, sep = "_") %>% 
      filter(cancer %in% grail_validation_samples) %>% 
      filter(normal %in% grail_validation_samples) %>% 
      dplyr::select(-c(cancer, normal, cancer_ratio, normal_ratio))
    
    # checking for data leakage
    test_normals_in_training <- mixed_training_df %>% 
      separate(sample, into = c("cancer", "normal", "cancer_ratio", "normal_ratio"), remove = FALSE, sep = "_") %>% 
      pull(normal) %>% 
      unique()
    test_cancer_in_training <- mixed_training_df %>% 
      separate(sample, into = c("cancer", "normal", "cancer_ratio", "normal_ratio"), remove = FALSE, sep = "_") %>% 
      pull(cancer) %>% 
      unique()
    test_samples_in_training <- c(test_normals_in_training, test_cancer_in_training) %>% unique()
    
    test_normals_in_validation <- mixed_validation_df %>% 
      separate(sample, into = c("cancer", "normal", "cancer_ratio", "normal_ratio"), remove = FALSE, sep = "_") %>%  
      pull(normal) %>% 
      unique()
    test_cancer_in_validation <- mixed_validation_df %>% 
      separate(sample, into = c("cancer", "normal", "cancer_ratio", "normal_ratio"), remove = FALSE, sep = "_") %>% 
      pull(cancer) %>% 
      unique()
    test_samples_in_validation <- c(test_normals_in_validation, test_cancer_in_validation) %>% unique()
    
    intersect(test_samples_in_training, test_samples_in_validation) # should yield no overlaps
    
    ### RUN ML ###
    # Large feature tables will take a long time to run
    # starttime <- Sys.time()
    print("Training...")
    mixed_final_fit <- trainmodel(mixed_training_df %>% dplyr::select(-c(ctdna_fraction, sample)), wt=T, seed2use = SEED)
    # traintime <- Sys.time() - starttime
    print("Validating...")
    mixed_validation_preds <- predictmodel(mixed_final_fit, mixed_validation_df) %>% 
      mutate(
        feature_table = ft,
        seedused = SEED
      )
    
    # Extract ROCs by ctdna_fraction bin
    mixed_validation_ROCs <- extract_ROC_by_ctf(predictions = mixed_validation_preds, feature_table = ft, seedused = SEED)

    # combine data
    grail_mixed_validation_preds <- bind_rows(grail_mixed_validation_preds, mixed_validation_preds)
    grail_mixed_validation_ROCS <- bind_rows(grail_mixed_validation_ROCS, mixed_validation_ROCs)

    # save temp file incase it crashed part way through
    saveRDS(grail_mixed_validation_preds, "output/machine_learning/GRAIL_mixing_predictions_MOSTRECENT.rds")
    saveRDS(grail_mixed_validation_ROCS, "output/machine_learning/GRAIL_mixing_AUROCs_MOSTRECENT.rds")
    
    if (SEED == max(seeds)) {
      outpath <- str_c("output/machine_learning/GRAIL_mixing_prediction_", ft, ".rds")
      saveRDS(grail_mixed_validation_preds, outpath)
      if (ft == "TFBS_entropy") {
        saveRDS(grail_mixed_validation_preds, "temp/mixed_metrics_mixedNorm100M_seed16-33_preds_upsample7030_FINAL.rds")
        saveRDS()
      }
    }
    print(str_c(Sys.time(), "Finished!", sep = " "))
  }
}
elapsedtime <- Sys.time() - starttime
print(elapsedtime)
