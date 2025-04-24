# Get ROC values for GRAIL downsampling
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

# Extract UW ROC values
extract_ROC_UW <- function(predictions, feature_table, gene_panel, seedused, RCC = FALSE, ERsub = FALSE, Normals = FALSE, Lungsub = FALSE) {
  if (!RCC & !ERsub & !Normals & !Lungsub) {
    rocblca <- roc(predictions, actualBLCA, `pred_Bladder`)
    roclung <- roc(predictions, actualLUNG, `pred_Lung`)
    rocbrca <- roc(predictions, actualBRCA, `pred_Breast`)
    rocprad <- roc(predictions, actualPRAD, `pred_Prostate`)
    rocnepc <- roc(predictions, actualNEPC, `pred_NEPC`)
    rocrcc <- roc(predictions, actualRCC, `pred_RCC`)
    
    output <- tibble(
      feature_table = feature_table,
      gene_panel = gene_panel,
      seed = seedused,
      phenotype = c("Bladder", "Lung", "Breast", "Prostate", "NEPC", "RCC"),
      roc = c(rocblca$auc[1], roclung$auc[1], rocbrca$auc[1], rocprad$auc[1], rocnepc$auc[1], rocrcc$auc[1])
    )
  }
  
  if (RCC & ERsub & !Normals & !Lungsub) {
    rocblca <- roc(predictions, actualBLCA, `pred_Bladder`)
    roclung <- roc(predictions, actualLUNG, `pred_Lung`)
    rocERpos <- roc(predictions, actualERpos, `pred_ERpos`)
    rocERneg <- roc(predictions, actualERneg, `pred_ERneg`)
    rocprad <- roc(predictions, actualPRAD, `pred_Prostate`)
    rocnepc <- roc(predictions, actualNEPC, `pred_NEPC`)
    rocrcc <- roc(predictions, actualRCC, `pred_RCC`)
    
    output <- tibble(
      feature_table = feature_table,
      gene_panel = gene_panel,
      seed = seedused,
      phenotype = c("Bladder", "Lung", "ERpos", "ERneg", "Prostate", "NEPC", "RCC"),
      roc = c(rocblca$auc[1], roclung$auc[1], rocERpos$auc[1], rocERneg$auc[1], rocprad$auc[1], rocnepc$auc[1], rocrcc$auc[1])
    )
  }
  
  if (RCC & ERsub & Normals & !Lungsub) {
    rocblca <- roc(predictions, actualBLCA, `pred_Bladder`)
    roclung <- roc(predictions, actualLUNG, `pred_Lung`)
    rocERpos <- roc(predictions, actualERpos, `pred_ERpos`)
    rocERneg <- roc(predictions, actualERneg, `pred_ERneg`)
    rocprad <- roc(predictions, actualPRAD, `pred_Prostate`)
    rocnepc <- roc(predictions, actualNEPC, `pred_NEPC`)
    rocrcc <- roc(predictions, actualRCC, `pred_RCC`)
    rocnon <- roc(predictions, actualNON, `pred_Normal`)
    
    output <- tibble(
      feature_table = feature_table,
      gene_panel = gene_panel,
      seed = seedused,
      phenotype = c("Bladder", "Lung", "ERpos", "ERneg", "Prostate", "NEPC", "RCC", "Normal"),
      roc = c(rocblca$auc[1], roclung$auc[1], rocERpos$auc[1], rocERneg$auc[1], rocprad$auc[1], rocnepc$auc[1], rocrcc$auc[1], rocnon$auc[1])
    )
  }
  
  if (RCC & ERsub & Normals & Lungsub) {
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
  }
  
  return(output)
}


extract_ROC_predictCTF <- function(predictions, feature_table, gene_panel, seedused) {
  rochigh <- roc(predictions, actual_high, pred_High)
  rocmid <- roc(predictions, actual_mid, pred_Mid)
  roclow <- roc(predictions, actual_low, pred_Low)
  rochealthy <- roc(predictions, actual_healthy, pred_Healthy)
  
  output <- tibble(
    feature_table = feature_table,
    gene_panel = gene_panel,
    seed = seedused,
    phenotype = c("High", "Mid", "Low", "Healthy"),
    roc = c(rochigh$auc[1], rocmid$auc[1], roclow$auc[1], rochealthy$auc[1])
  )
  
  return(output)
}

extract_ROC_predictCTF_2 <- function(predictions, feature_table, gene_panel, seedused) {
  rochigh <- roc(predictions, actual_high, pred_High)
  # rocmid <- roc(predictions, actual_mid, pred_Mid)
  roclow <- roc(predictions, actual_low, pred_Low)
  rochealthy <- roc(predictions, actual_healthy, pred_Healthy)
  
  output <- tibble(
    feature_table = feature_table,
    gene_panel = gene_panel,
    seed = seedused,
    phenotype = c("High", "Low", "Healthy"),
    roc = c(rochigh$auc[1], roclow$auc[1], rochealthy$auc[1])
  )
  
  return(output)
}

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
              "green2",
              "gold2")
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

plotROC_full <- function(preds, filename='ROC.pdf', usegrail) {
  if(usegrail) {
    rocnon <- roc(preds, actualNON, `pred_Normal`)
    roclung <- roc(preds, actualLUNG, `pred_Lung`)
    rocbrca <- roc(preds, actualBRCA, `pred_Breast`)
    rocprad <- roc(preds, actualPRAD, `pred_Prostate`)
    
    pdf(filename,5,5)
    plot.roc(rocnon,print.auc=T,col=colorpal['Normal'],print.auc.pattern="Normal AUC: %.2f") # Changed "Cancer AUC" to "Normal AUC"
    plot.roc(roclung,add=T,print.auc=T,col=colorpal['Lung'],print.auc.y=0.45,print.auc.pattern="Lung AUC: %.2f")
    plot.roc(rocbrca,add=T,print.auc=T,col=colorpal['Breast'],print.auc.y=0.4,print.auc.pattern="Breast AUC: %.2f")
    plot.roc(rocprad,add=T,print.auc=T,col=colorpal['Prostate'],print.auc.y=0.35,print.auc.pattern="Prostate AUC: %.2f")
  } else {
    rocbladder <- roc(preds, actualBLCA, `pred_Bladder`)
    # roclung <- roc(preds, actualLUNG, `pred_Lung`)
    # rocbrca <- roc(preds, actualBRCA, `pred_Breast`)
    rocprad <- roc(preds, actualPRAD, `pred_Prostate`)
    rocnepc <- roc(preds, actualNEPC, `pred_NEPC`)
    rocrcc <- roc(preds, actualRCC, `pred_RCC`)
    rocnormal <- roc(preds, actualNON, pred_Normal)
    rocERpos <- roc(preds, actualERpos, `pred_ERpos`)
    rocERneg <- roc(preds, actualERneg, `pred_ERneg`)
    rocsclc <- roc(preds, actualSCLC, `pred_SCLC`)
    rocnsclc <- roc(preds, actualNSCLC, `pred_NSCLC`)
    
    pdf(filename,5,5)
    # plot.roc(rocbladder,print.auc=T,col=colorpal['Bladder'],print.auc.pattern="Bladder AUC: %.2f")
    # plot.roc(roclung,add=T,print.auc=T,col=colorpal['Lung'],print.auc.y=0.45,print.auc.pattern="Lung AUC: %.2f")
    # plot.roc(rocbrca,add=T,print.auc=T,col=colorpal['Breast'],print.auc.y=0.4,print.auc.pattern="Breast AUC: %.2f")
    # plot.roc(rocprad,add=T,print.auc=T,col=colorpal['Prostate'],print.auc.y=0.35,print.auc.pattern="Prostate AUC: %.2f")
    # plot.roc(rocnepc,add=T,print.auc=T,col=colorpal['NEPC'],print.auc.y=0.3,print.auc.pattern="NEPC AUC: %.2f")
    # plot.roc(rocrcc,add=T,print.auc=T,col=colorpal['RCC'],print.auc.y=0.25,print.auc.pattern="RCC AUC: %.2f")
    
    
    plot.roc(rocnormal, print.auc=T, col=colorpal['Normal'], print.auc.pattern="Normal AUC: %.2f")
    plot.roc(rocbladder, add=T, print.auc=T, col=colorpal['Bladder'], print.auc.y=0.45, print.auc.pattern="Bladder AUC: %.2f")
    plot.roc(rocprad, add=T, print.auc=T, col=colorpal['Prostate'], print.auc.y=0.4, print.auc.pattern="Prostate AUC: %.2f")
    plot.roc(rocnepc, add=T, print.auc=T, col=colorpal['NEPC'], print.auc.y=0.35, print.auc.pattern="NEPC AUC: %.2f")
    plot.roc(rocrcc, add=T, print.auc=T, col=colorpal['RCC'], print.auc.y=0.3, print.auc.pattern="RCC AUC: %.2f")
    plot.roc(rocERpos, add=T, print.auc=T, col=colorpal['ERpos'], print.auc.y=0.25, print.auc.pattern="ER+ AUC: %.2f")
    plot.roc(rocERneg, add=T, print.auc=T, col=colorpal['ERneg'], print.auc.y=0.2, print.auc.pattern="ER- AUC: %.2f")
    plot.roc(rocsclc, add=T, print.auc=T, col=colorpal['SCLC'], print.auc.y=0.15, print.auc.pattern="SCLC AUC: %.2f")
    plot.roc(rocnsclc, add=T, print.auc=T, col=colorpal['NSCLC'], print.auc.y=0.1, print.auc.pattern="NSCLC AUC: %.2f")
  }
  dev.off()
}

# test_roc <- UW_comparison_results %>% 
#   filter(feature_table == "all_exons_SE_and_depth")
# 
# plotROC_full(preds = test_roc, filename = "temp/UW_new_model_ROC.pdf", usegrail = FALSE)

# read_rds("feature_comparison_scripts/data_output/grail_samples/GRAIL_comparison_ROCs_2023_10_03_18_04_34.rds") %>% 
#   count(feature_table)
# read_rds("feature_comparison_scripts/data_output/UW_all_samples_normal_BCsub_Lungsub/UW_comparison_ROCs_2023_10_09_14_55_37.rds") %>% 
#   count(feature_table, phenotype)
# read_rds("predict_ctf/data_output/grail/GRAIL_comparison_ROCs_predict_ctf_2023_10_08_12_17_46.rds") %>% 
#   count(feature_table)
# read_rds("predict_ctf/data_output/uw/UW_comparison_ROCs_predict_ctf_2023_10_10_19_39_34.rds") %>% 
#   count(feature_table)
# 
# 
# test_roc %>%
#   filter(!is.na(ctfraction)) %>%
#   # filter(between(ctfraction, 0.05, 1)) %>% 
#   count(actual)

extract_roc_by_ctf <- function(predictions, low, high, FEATURE_TABLE, GENE_PANEL, SEED) {
  output <- tibble()
  BIN <- str_c(low, "-", high)
  all_phenos <- predictions %>% count(actual) %>% pull(actual)
  
  filtered_preds <- predictions %>% 
    filter(feature_table == FEATURE_TABLE) %>% 
    filter(gene_panel == GENE_PANEL) %>% 
    filter(seed == SEED) %>% 
    filter(!is.na(ctfraction)) %>% 
    filter(between(ctfraction, low, high))
  
  filtered_phenos <- filtered_preds %>% count(actual) %>% pull(actual)
  
  if ("Bladder" %in% filtered_phenos) {
    roc_bladder <- roc(filtered_preds, actualBLCA, pred_Bladder)
    temp_output <- tibble(phenotype = "Bladder", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_bladder$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  if ("ERpos" %in% filtered_phenos) {
    roc_ERpos <- roc(filtered_preds, actualERpos, pred_ERpos)
    temp_output <- tibble(phenotype = "ERpos", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_ERpos$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  if ("ERneg" %in% filtered_phenos) {
    roc_ERneg <- roc(filtered_preds, actualERneg, pred_ERneg)
    temp_output <- tibble(phenotype = "ERneg", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_ERneg$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  if ("Prostate" %in% filtered_phenos) {
    roc_prostate <- roc(filtered_preds, actualPRAD, pred_Prostate)
    temp_output <- tibble(phenotype = "Prostate", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_prostate$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  if ("RCC" %in% filtered_phenos) {
    roc_RCC <- roc(filtered_preds, actualRCC, pred_RCC)
    temp_output <- tibble(phenotype = "RCC", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_RCC$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  if ("SCLC" %in% filtered_phenos) {
    roc_SCLC <- roc(filtered_preds, actualSCLC, pred_SCLC)
    temp_output <- tibble(phenotype = "SCLC", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_SCLC$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  if ("NSCLC" %in% filtered_phenos) {
    roc_NSCLC <- roc(filtered_preds, actualNSCLC, pred_NSCLC)
    temp_output <- tibble(phenotype = "NSCLC", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_NSCLC$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  if ("NEPC" %in% filtered_phenos) {
    roc_NEPC <- roc(filtered_preds, actualNEPC, pred_NEPC)
    temp_output <- tibble(phenotype = "NEPC", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_NEPC$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  if ("Breast" %in% filtered_phenos) {
    roc_BRCA <- roc(filtered_preds, actualBRCA, pred_Breast)
    temp_output <- tibble(phenotype = "Breast", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_BRCA$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  if ("Lung" %in% filtered_phenos) {
    roc_LUNG <- roc(filtered_preds, actualLUNG, pred_Lung)
    temp_output <- tibble(phenotype = "Lung", feature_table = FEATURE_TABLE, gene_panel = GENE_PANEL, seed = SEED, bin = BIN, roc = roc_LUNG$auc[1])
    output <- bind_rows(output, temp_output)
  }
  
  return(output)
}

# extract_roc_by_ctf(predictions = test_roc, low = 0.5, high = 1, feature_table = "E1SE", gene_panel = "default", seed = 1)
# 
# UW_comparison_results %>% 
#   count(feature_table, gene_panel)

# roc_by_ctf <- tibble()
# for (f in UW_comparison_results$feature_table %>% unique()) {
#   for (g in UW_comparison_results$gene_panel %>% unique()) {
#     for (s in UW_comparison_results$seed %>% unique()) {
#       print(str_c("Extracting", f, g, s, sep = " "))
#       temp_roc_output <- extract_roc_by_ctf(predictions = UW_comparison_results,
#                                             low = 0,
#                                             high = 0.05,
#                                             feature_table = f,
#                                             gene_panel = g,
#                                             seed = s)
#       roc_by_ctf <- bind_rows(roc_by_ctf, temp_roc_output)
#     }
#   }
# }

plot_radar <- function(predictions, grail = usegrail, colors = colorpal, alpha = 1, plot_wrong_grey = FALSE) {
  
  # prefix for outfiles
  if (grail) {
    file_prefix <- "grail"
  } else {
    file_prefix <- "UW"
  }
  
  # get unique cancer types
  cancer_types <- predictions %>% pull(actual) %>% unique()
  
  # keeps colors consistent betweeen datasets
  # color_grid <- tibble(
  #   cancer = c("Normal", "Prostate", "Breast", "Lung", "RCC", "Bladder", "NEPC"),
  #   color_num = seq(1,7)
  # )
  
  # make radar plot for each cancer type
  for (cancer in cancer_types) {
    # plot_color <- color_grid$color_num[color_grid$cancer == cancer]
    plot_color <- colors[cancer][[1]]
    
    radar_data <- predictions %>% 
      filter(actual == cancer) %>% 
      dplyr::select(1:pred_class) %>% 
      dplyr::select(-pred_class) %>% 
      as.data.frame()
    
    color_vector <- predictions %>% 
      filter(actual == cancer) %>% 
      mutate(
        colors = ifelse(pred_class == actual, plot_color, "#808080")
      ) %>%
      pull(colors)
    
    colnames(radar_data) <- gsub("pred_", "", colnames(radar_data))
    
    num_classes <- length(colnames(radar_data))
    
    radar_data <- rbind(rep(1, num_classes), rep(0, num_classes), radar_data)
    
    outfile <- str_c(file_prefix, "_", cancer, "_radar_plot.pdf")
    
    pdf(outfile, width = 5, height = 5)
    if (plot_wrong_grey) {
      radarchart(radar_data, pcol = alpha(color_vector, alpha), plty = 1, title = cancer)
    } else {
      radarchart(radar_data, pcol = alpha(plot_color, alpha), plty = 1, title = cancer)
    }
    dev.off()
  }
}