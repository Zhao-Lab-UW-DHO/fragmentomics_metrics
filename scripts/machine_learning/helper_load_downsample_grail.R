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
    add_column(actualSCLC=pred$actual=='SCLC') %>%
    add_column(actualPANC=pred$actual=='Pancreas') %>%
    add_column(actualRCC=pred$actual=='RCC') %>% 
    add_column(actualERpos = pred$actual == "ERpos") %>% 
    add_column(actualERneg = pred$actual == "ERneg")
  
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
             'blue')
names(colorpal) <- c('All',
                     'Bladder',
                     'Breast',
                     'Lung',
                     'NEPC',
                     'Normal',
                     'Prostate',
                     'RCC',
                     'ERpos',
                     'ERneg')
colorsUW <- colorpal[c(-6,-8)]
colorsGRAIL <- colorpal[c(-2,-5,-8)]

file_grail <- feature_table_path
print(file_grail)
grail_df <- readRDS(file_grail)

filemetagrail <- 'files/grail_ml_metadata.rds'
print(filemetagrail)
grailmetadf <- readRDS(filemetagrail)
stopifnot(identical(grailmetadf$sample,grail_df$sample))


grail_df$phenotype <- grailmetadf$phenotype
grail_df$phenotype <- relevel(factor(grail_df$phenotype),ref='Normal')

grail_df$ctdna_fraction <- grailmetadf$ctdna_fraction
grail_df$sample <- grailmetadf$sample
rowsnum <- which(names(grail_df) %in% c('phenotype','ctdna_fraction','sample'))
grail_df[-rowsnum] <- data.frame(t(scale(t(grail_df[-rowsnum]))))
grail_df <- grail_df %>% 
  relocate(sample, phenotype, ctdna_fraction)


############################################################### 
###################### Printing/Plotting ######################
###############################################################

# Get ROC values for GRAIL downsampling
extract_ROC_grail_downsample <- function(predictions, downsample_level, replicate_number, feature, seedused) {
  rocnon <- roc(predictions, actualNON, `pred_Normal`)
  roclung <- roc(predictions, actualLUNG, `pred_Lung`)
  rocbrca <- roc(predictions, actualBRCA, `pred_Breast`)
  rocprad <- roc(predictions, actualPRAD, `pred_Prostate`)
  
  output <- tibble(
    replicate = replicate_number,
    downsample = downsample_level,
    feature_table = feature,
    seed = seedused,
    phenotype = c("Normal", "Lung", "Breast", "Prostate"),
    roc = c(rocnon$auc[1], roclung$auc[1], rocbrca$auc[1], rocprad$auc[1])
  )
  
  return(output)
}
