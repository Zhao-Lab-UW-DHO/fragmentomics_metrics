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

###############################################################
####################### Model Functions #######################
###############################################################
trainmodel <- function(data2train, seed2use=seed, tunenfold=nfold, tunenthreads=nthread, wt=T) {
  
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
      #step_upsample(phenotype)
      #step_smote(phenotype)
      #step_adasyn(phenotype)
      step_smotenc(phenotype)
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
    add_column(actualPANC=pred$actual=='Pancreas') %>%
    add_column(actualRCC=pred$actual=='RCC') %>% 
    add_column(actualERpos = pred$actual == "ERpos") %>% 
    add_column(actualERneg = pred$actual == "ERneg") %>% 
    add_column(actualSCLC = pred$actual == "SCLC") %>% 
    add_column(actualNSCLC = pred$actual == "NSCLC")
  
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
             'RCC',
             'Normal')
names(namesub) <- c('Cancer breast',
                    'Cancer bladder',
                    'Cancer prostate',
                    'Cancer lung',
                    'Cancer NEPC',
                    #'Cancer SCLC',
                    #'Cancer ovarian',
                    #'Cancer pancreas',
                    'Cancer rcc',
                    'Non cancer')

colorpal <- c('black',
             '#1B9E77',
             '#D95F02',
             '#7570B3',
             '#E7298A',
             '#666666',
             '#66A61E',
             '#A6761D',
             'red',
             'blue',
             "green",
             "gold")
names(colorpal) <- c('All',
                     'Bladder',
                     'Breast',
                     'Lung',
                     'NEPC',
                     'Normal',
                     'Prostate',
                     'RCC',
                     'ERpos',
                     'ERneg',
                     "SCLC",
                     "NSCLC")
colorsUW <- colorpal[c(-6,-8)]
colorsGRAIL <- colorpal[c(-2,-5,-8)]

# Read in Feature Table and Arrange by sample
file_internal <- feature_table_to_use
print(file_internal)
internal_df <- readRDS(file_internal) %>% arrange(sample)

# Read in Metadata and arrnage by sample
filemetainternal <- 'files/uw_ml_metadata.rds'
print(filemetainternal)
internalmetadf <- readRDS(filemetainternal) %>% arrange(sample)

# Check that Feature Table and Metadata have the same samples
stopifnot(identical(internalmetadf$sample,internal_df$sample))

# Join Feature Table with Metadata
internal_df <- internal_df %>% 
  left_join(internalmetadf, by = "sample") %>% 
  relocate(sample, phenotype, subtype, ctdna_fraction)

## Clean up data
rowsnepc <- internal_df$subtype=='NEPC'
rowsnepc[is.na(rowsnepc)] <- F
internal_df[rowsnepc,'phenotype'] <- 'Cancer NEPC'

internal_df <- internal_df %>% dplyr::select(-subtype)

# Rename Phenotypes to Simpler Names
internal_df$phenotype <- namesub[internal_df$phenotype]
internal_df$phenotype <- factor(internal_df$phenotype)

# Normalization of Data
rowsnum <- which(names(internal_df) %in% c('phenotype','ctdna_fraction', 'sample'))
internal_df[-rowsnum] <- data.frame(t(scale(t(internal_df[-rowsnum]))))


# Replace Breast with HR subtype and Lung with SCLC/NSCLC
ERpos_samples <- internalmetadf %>% filter(phenotype == "Cancer breast") %>% filter(str_detect(subtype, "ER\\+")) %>% pull(sample)
ERneg_samples <- internalmetadf %>% filter(phenotype == "Cancer breast") %>% filter(str_detect(subtype, "ER-")) %>% pull(sample)

SCLC_samples <- internalmetadf %>% filter(phenotype == "Cancer lung") %>% filter(subtype == "sclc") %>% pull(sample)
NSCLC_samples <- internalmetadf %>% filter(phenotype == "Cancer lung") %>% filter(subtype == "nsclc") %>% pull(sample)

pheno_levels <- c("Bladder", "ERpos", "ERneg", "NEPC", "Prostate", "RCC", "Normal", "SCLC", "NSCLC")

internal_df <- internal_df %>% 
  mutate(phenotype = as.character(phenotype)) %>% 
  mutate(phenotype = ifelse(sample %in% ERpos_samples, "ERpos", phenotype)) %>% 
  mutate(phenotype = ifelse(sample %in% ERneg_samples, "ERneg", phenotype)) %>% 
  mutate(phenotype = ifelse(sample %in% SCLC_samples, "SCLC", phenotype)) %>% 
  mutate(phenotype = ifelse(sample %in% NSCLC_samples, "NSCLC", phenotype)) %>% 
  mutate(phenotype = factor(phenotype, levels = pheno_levels))

# Print number of samples in each phenotype
print(table(internal_df$phenotype))


############################################################### 
###################### Printing/Plotting ######################
###############################################################

# Extract UW ROC values
extract_ROC_UW <- function(predictions, feature_table, gene_panel, seedused) {
  
    rocblca <- roc(predictions, actualBLCA, `pred_Bladder`)
    # roclung <- roc(predictions, actualLUNG, `pred_Lung`)
    rocERpos <- roc(predictions, actualERpos, `pred_ERpos`)
    rocERneg <- roc(predictions, actualERneg, `pred_ERneg`)
    rocprad <- roc(predictions, actualPRAD, `pred_Prostate`)
    rocnepc <- roc(predictions, actualNEPC, `pred_NEPC`)
    rocrcc <- roc(predictions, actualRCC, `pred_RCC`)
    rocnon <- roc(predictions, actualNON, `pred_Normal`)
    rocsclc <- roc(predictions, actualSCLC, `pred_SCLC`)
    rocnsclc <- roc(predictions, actualNSCLC, `pred_NSCLC`)
    
    output <- tibble(
      feature_table = feature_table,
      gene_panel = gene_panel,
      seed = seedused,
      phenotype = c("Bladder", "ERpos", "ERneg", "Prostate", "NEPC", "RCC", "Normal", "SCLC", "NSCLC"),
      roc = c(rocblca$auc[1], rocERpos$auc[1], rocERneg$auc[1], rocprad$auc[1], rocnepc$auc[1], rocrcc$auc[1], rocnon$auc[1], rocsclc$auc[1], rocnsclc$auc[1])
    )
  
  return(output)
}

extract_ROC <- function(predictions, feature_table, gene_panel, seedused) {
  rocnon <- roc(predictions, actualNON, `pred_Normal`)
  roclung <- roc(predictions, actualLUNG, `pred_Lung`)
  rocbrca <- roc(predictions, actualBRCA, `pred_Breast`)
  rocprad <- roc(predictions, actualPRAD, `pred_Prostate`)
  
  output <- tibble(
    feature_table = feature_table,
    gene_panel = gene_panel,
    seed = seedused,
    phenotype = c("Normal", "Lung", "Breast", "Prostate"),
    roc = c(rocnon$auc[1], roclung$auc[1], rocbrca$auc[1], rocprad$auc[1])
  )
  
  return(output)
}
