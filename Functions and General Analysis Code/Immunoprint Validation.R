# Validation 
easypackages::libraries(c("tidyr","dplyr","data.table","cutpointr","readr"))

# insert paths to folders
directory_path <- here::here() # path to directory enclosing code and data sets
folder_path <- paste0(directory_path, "/Supplemental Tables/") # folder containing downloaded data sets
intermediate_folder_path <- paste0(directory_path, "/Intermediate Files Used in Plots/") # path containing intermediate outputs
folder_path_code <- paste0(directory_path, "/Functions and General Analysis Code/") # change this to folder containing downloaded functions and code


source(paste0(folder_path_code, "Immunoprint functions.R"))

# Load Data
combined_metadata_dist <- read_csv(paste0(folder_path, "Data Table S1.csv"))
all_zs_combined <- read_csv(paste0(folder_path, "Data Table S2.csv"))
Viral_Annotations <- read_csv(paste0(folder_path, "Data Table S4.csv"))

# Subset to validation data set
Imm_Val_Zs <- subset_frame_fxn(all_zs_combined, combined_metadata_dist, "imm_val")
Imm_Val_Metadata <- combined_metadata_dist %>% filter(imm_val == TRUE)

# Ensure metadata identifiers exist in matrix (safety)
val_ids <- setdiff(colnames(Imm_Val_Zs), "id")
Imm_Val_Metadata <- Imm_Val_Metadata %>% filter(identifier %in% val_ids)

# Bootstrap settings
n_boot <- 5000

# Preallocate results
Validation_Results <- data.frame(
  Dataset = rep("Validation", n_boot),
  Method  = rep("Final immunoprint", n_boot),
  TP = NA_real_, FP = NA_real_, TN = NA_real_, FN = NA_real_,
  sum_sens_spec = NA_real_,
  accuracy = NA_real_,
  sensitivity = NA_real_,
  specificity = NA_real_,
  AUC = NA_real_,
  nSamples = rep(length(unique(Imm_Val_Metadata$identifier)), n_boot),
  nIndividuals = rep(length(unique(Imm_Val_Metadata$donor)), n_boot)
)

# Helper: compute TP/FP/TN/FN directly from threshold
# Assumes higher immunoprint_match => "Same individual"
compute_cm_from_threshold <- function(df, score_col, label_col, threshold,
                                      pos_label = "Same individual",
                                      neg_label = "Different individuals") {
  score <- df[[score_col]]
  lab   <- df[[label_col]]
  
  # predicted class
  pred <- ifelse(score >= threshold, pos_label, neg_label)
  
  # ensure all four cells exist even if zero
  tab <- table(
    factor(pred, levels = c(pos_label, neg_label)),
    factor(lab,  levels = c(pos_label, neg_label))
  )
  
  tp <- tab[pos_label, pos_label]
  fp <- tab[pos_label, neg_label]
  fn <- tab[neg_label, pos_label]
  tn <- tab[neg_label, neg_label]
  
  c(tp = as.numeric(tp), fp = as.numeric(fp), tn = as.numeric(tn), fn = as.numeric(fn))
}

# Run bootstrap: compare 1/3 (new samples) to 2/3 (database samples)
for (i in seq_len(n_boot)) {
  set.seed(10 + i)
  
  database_samps <- Imm_Val_Metadata %>%
    distinct(identifier) %>%
    sample_frac(2/3) %>%
    pull(identifier)
  
  new_samps <- Imm_Val_Metadata %>%
    distinct(identifier) %>%
    filter(!identifier %in% database_samps) %>%
    pull(identifier)
  
  # Guardrails (skip degenerate splits)
  if (length(database_samps) < 2 || length(new_samps) < 1) next
  
  immunoprint_Val <- immunoprint(
    Imm_Val_Zs %>% select(id, all_of(new_samps)),
    Imm_Val_Zs %>% select(id, all_of(database_samps)),
    Viral_Annotations,
    keepAll = TRUE
  )
  
  immunoprint_Val <- immunoprint_Val %>%
    left_join(
      Imm_Val_Metadata %>% mutate(sample1 = identifier, indiv1 = donor) %>% select(sample1, indiv1),
      by = "sample1"
    ) %>%
    left_join(
      Imm_Val_Metadata %>% mutate(sample2 = identifier, indiv2 = donor) %>% select(sample2, indiv2),
      by = "sample2"
    ) %>%
    mutate(label = ifelse(indiv1 == indiv2, "Same individual", "Different individuals"))
  
  # Evaluate performance
  cp <- cutpointr(
    immunoprint_Val,
    immunoprint_match,
    label,
    method = maximize_metric,
    metric = sum_sens_spec
  )
  
  thr <- cp$optimal_cutpoint
  
  # Version-proof confusion matrix counts (computed directly)
  cm <- compute_cm_from_threshold(
    immunoprint_Val,
    score_col = "immunoprint_match",
    label_col = "label",
    threshold = thr,
    pos_label = "Same individual",
    neg_label = "Different individuals"
  )
  
  Validation_Results[i, c("TP","FP","TN","FN")] <- c(cm["tp"], cm["fp"], cm["tn"], cm["fn"])
  
  # Keep cutpointr metrics (these should now agree with TP/FP/TN/FN)
  Validation_Results[i, c("sum_sens_spec","accuracy","sensitivity","specificity","AUC")] <-
    c(cp$sum_sens_spec, cp$acc, cp$sensitivity, cp$specificity, cp$AUC)
  
  rm(cp)
}

# Summarize (means across bootstrap iterations)
metric_cols <- c("TP","FP","TN","FN","sum_sens_spec","accuracy","sensitivity","specificity","AUC","nSamples","nIndividuals")

Validation_Result_Summary <- data.frame(
  Dataset = "Validation",
  Method  = "Final immunoprint"
)
Validation_Result_Summary[metric_cols] <- lapply(
  metric_cols,
  function(col) mean(Validation_Results[[col]], na.rm = TRUE)
)

write.csv(
  Validation_Results,
  file = paste0(intermediate_folder_path, "Immunoprint_Validation_Results.csv"),
  row.names = FALSE
)

write.csv(
  Validation_Result_Summary,
  file = paste0(intermediate_folder_path, "Immunoprint_Validation_Results_Summary.csv"),
  row.names = FALSE
)

# Also save one immunoprint validation iteration for plot
set.seed(10)

database_samps <- Imm_Val_Metadata %>%
  distinct(identifier) %>%
  sample_frac(2/3) %>%
  pull(identifier)

new_samps <- Imm_Val_Metadata %>%
  distinct(identifier) %>%
  filter(!identifier %in% database_samps) %>%
  pull(identifier)

immunoprint_Val <- immunoprint(
  Imm_Val_Zs %>% select(id, all_of(new_samps)),
  Imm_Val_Zs %>% select(id, all_of(database_samps)),
  Viral_Annotations,
  keepAll = TRUE
)

immunoprint_Val <- immunoprint_Val %>%
  left_join(
    Imm_Val_Metadata %>% mutate(sample1 = identifier, indiv1 = donor) %>% select(sample1, indiv1),
    by = "sample1"
  ) %>%
  left_join(
    Imm_Val_Metadata %>% mutate(sample2 = identifier, indiv2 = donor) %>% select(sample2, indiv2),
    by = "sample2"
  ) %>%
  mutate(label = ifelse(indiv1 == indiv2, "Same individual", "Different individuals"))

write.csv(
  immunoprint_Val,
  file = paste0(intermediate_folder_path, "Immunoprint_Validation_seed10.csv"),
  row.names = FALSE
)
