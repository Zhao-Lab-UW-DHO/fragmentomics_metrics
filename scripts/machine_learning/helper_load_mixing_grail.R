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
library(gghighlight)
library(scales)


options(scipen=999)

nfold <- 10
nthread <- 10

# zhao_panel_v2 <- read_tsv("~/Downloads/panel_v2_zhao_mappable_cbioportal_N542.txt") %>% pull(x)

###############################################################
####################### Model Functions #######################
###############################################################
trainmodel <- function(data2train, seed2use=seed, tunenfold=nfold, tunenthreads=nthread, wt=T) {
  #data2train$ctdna_fraction[is.na(data2train$ctdna_fraction)] <- median(data2train$ctdna_fraction,na.rm=T)
  #data2train$ctdna_fraction <- importance_weights(data2train$ctdna_fraction)
  
  model <- multinom_reg(
    mode = "classification",
    engine = "glmnet",
    penalty = tune(),
    mixture = tune()
  )
  
  grid <- grid_latin_hypercube(
    penalty(),
    mixture(),
    size = nfold
  )
  
  rec <- recipe(phenotype ~ .,data=data2train)
  
  if(wt) {
    rec <- rec %>%
      step_upsample(phenotype)
      #step_smote(phenotype)
      #step_adasyn(phenotype)
      # step_smotenc(phenotype)
  }
  
  wflow <- workflow() %>%
    add_model(model) %>% 
    add_recipe(rec) #%>%
    #add_case_weights(ctdna_fraction)
  
  set.seed(seed2use)
  tunefolds <- vfold_cv(data2train, v=tunenfold, repeats=1)
  
  cl <- makeCluster(nthread)
  registerDoSNOW(cl)
  tunegrid <- tune_grid(
    wflow,
    resamples = tunefolds,
    grid = grid,
    control = control_grid(save_pred = TRUE)
  )
  
  stopCluster(cl)
  
  # best_auc <- select_best(tunegrid, "roc_auc")
  best_auc <- select_best(tunegrid, "accuracy")
  best_wflow <- finalize_workflow(
    wflow,
    best_auc
  )
  
  final_fit <- best_wflow %>% fit(data = data2train)
  
  return(final_fit)
}

predictmodel <- function(modelfit, data2validate) {
  predcont <- predict(modelfit,data2validate,type='prob')
  predcat <- predict(modelfit,data2validate)
  pred <- bind_cols(predcont,predcat)
  pred <- pred %>% add_column(actual=data2validate$phenotype)
  
  if('ctdna_fraction' %in% names(data2validate)) {
    pred <- pred %>%
      add_column(ctfraction=data2validate$ctdna_fraction)
  }
  
  if('sample' %in% names(data2validate)) {
    pred <- pred %>%
      add_column(sample=data2validate$sample)
  }
  
  pred <- pred %>%
    add_column(actualNON=pred$actual=='Normal') %>%
    add_column(actualBLCA=pred$actual=='Bladder') %>%
    add_column(actualLUNG=pred$actual=='Lung') %>%
    add_column(actualBRCA=pred$actual=='Breast') %>%
    add_column(actualPRAD=pred$actual=='Prostate') %>%
    add_column(actualNEPC=pred$actual=='NEPC') %>%
    add_column(actualSCLC=pred$actual=='SCLC') %>%
    add_column(actualPANC=pred$actual=='Pancreas') %>%
    add_column(actualRCC=pred$actual=='RCC')
  
  pred <- rename_with(pred, ~ gsub(".", "", .x, fixed = T),starts_with(".pred"))
  pred <- rename_with(pred, ~ gsub(" ", "_", .x, fixed = T),starts_with("pred"))
  return(pred)
}


############################################################### 
####################### Load Exon 1 Data ######################
###############################################################
namesub <- c('Breast',
             'Bladder',
             'Prostate',
             'Lung',
             'NEPC',
             #'SCLC',
             #'Ovarian',
             #'Pancreas',
             #'RCC',
             'Normal')
names(namesub) <- c('Cancer breast',
                    'Cancer bladder',
                    'Cancer prostate',
                    'Cancer lung',
                    'Cancer NEPC',
                    #'Cancer SCLC',
                    #'Cancer ovarian',
                    #'Cancer pancreas',
                    #'Cancer rcc',
                    'Non cancer')

colorpal <- c('black',
             '#1B9E77',
             '#D95F02',
             '#7570B3',
             '#E7298A',
             '#666666',
             '#66A61E',
             '#A6761D')
names(colorpal) <- c('All',
                     'Bladder',
                     'Breast',
                     'Lung',
                     'NEPC',
                     'Normal',
                     'Prostate',
                     'RCC')
colorsUW <- colorpal[c(-6,-8)]
colorsGRAIL <- colorpal[c(-2,-5,-8)]


grail_mixing_ft <- read_rds(grail_mixing_ft_path)
print(str_c("Loading: ", grail_mixing_ft_path))
grail_metadata_simple <- read_rds(grail_metadata_simple_path)
print(str_c("Loading: ", grail_metadata_simple_path))

print(str_c("Tidying and Normalizing mixed feature table: ", grail_mixing_ft_path))
grail_mixing_ft <- grail_mixing_ft %>% 
  separate(sample, into = c("cancer_sample", "normal_sample", "cancer_pct", "normal_pct"), remove = FALSE, sep = "_") %>% 
  left_join(grail_metadata_simple, by = c("cancer_sample" = "sample")) %>% 
  dplyr::rename(original_ctdna_fraction = ctdna_fraction) %>% 
  mutate(
    cancer_pct = as.double(cancer_pct),
    normal_pct = as.double(normal_pct),
    ctdna_fraction = original_ctdna_fraction * (cancer_pct / 100),
    phenotype = factor(phenotype, levels = c("Normal", "Breast", "Lung", "Prostate"))
    ) %>% 
  relocate(sample, cancer_sample, normal_sample, cancer_pct, normal_pct, phenotype, ctdna_fraction, original_ctdna_fraction) %>% 
  dplyr::select(-c(cancer_sample:normal_pct, original_ctdna_fraction)) %>% 
  dplyr::select(where(~ !any(is.na(.)))) %>% 
  pivot_longer(cols = -c(sample:ctdna_fraction), names_to = "feature", values_to = "value") %>% 
  group_by(sample, phenotype) %>% 
  mutate(value_norm = (value - mean(value)) / sd(value)) %>% 
  dplyr::select(-value) %>% 
  pivot_wider(names_from = feature, values_from = value_norm) %>% 
  ungroup()


############################################################### 
###################### Printing/Plotting ######################
###############################################################

extract_ROC <- function(predictions, feature_table, seedused) {
  rocnon <- roc(predictions, actualNON, `pred_Normal`)
  roclung <- roc(predictions, actualLUNG, `pred_Lung`)
  rocbrca <- roc(predictions, actualBRCA, `pred_Breast`)
  rocprad <- roc(predictions, actualPRAD, `pred_Prostate`)
  
  output <- tibble(
    replicate = seedused,
    feature_table = feature_table,
    phenotype = c("Normal", "Lung", "Breast", "Prostate"),
    roc = c(rocnon$auc[1], roclung$auc[1], rocbrca$auc[1], rocprad$auc[1])
  )
  
  return(output)
}

extract_ROC_by_ctf <- function(predictions, feature_table, seedused) {
  # All Data
  all_rocs <- extract_ROC(predictions = predictions, 
                          feature_table = feature_table, 
                          seedused = seedused) %>% 
    mutate(ctdna_group = "all")
  
  tmp1 <- extract_ROC(predictions = predictions %>% filter(between(ctfraction, 0.1, 1) | actual == "Normal"), 
                      feature_table = feature_table, 
                      seedused = seedused) %>% 
    mutate(ctdna_group = "0.1 - 1")
  
  tmp2 <- extract_ROC(predictions = predictions %>% filter(between(ctfraction, 0.01, 0.1) | actual == "Normal"), 
                      feature_table = feature_table, 
                      seedused = seedused) %>% 
    mutate(ctdna_group = "0.01 - 0.1")
  
  tmp3 <- extract_ROC(predictions = predictions %>% filter(between(ctfraction, 0.001, 0.01) | actual == "Normal"), 
                      feature_table = feature_table, 
                      seedused = seedused) %>% 
    mutate(ctdna_group = "0.001 - 0.01")
  
  tmp4 <- extract_ROC(predictions = predictions %>% filter(between(ctfraction, 0.0001, 0.001) | actual == "Normal"), 
                      feature_table = feature_table, 
                      seedused = seedused) %>% 
    mutate(ctdna_group = "0.0001 - 0.001")
  
  output <- bind_rows(all_rocs, tmp1, tmp2, tmp3, tmp4)
  return(output)
}
