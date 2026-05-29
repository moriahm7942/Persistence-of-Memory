
# Load packages
easypackages::libraries(c(
  "tidyr","dplyr","data.table","cutpointr","readr", "here"
))

# Paths
# Define paths relative to the root of the cloned GitHub repository
directory_path <- here::here()
folder_path <- paste0(directory_path, "/Supplemental Tables/") # folder containing downloaded data sets
intermediate_folder_path <- paste0(directory_path, "/Intermediate Files Used in Plots/") # path containing intermediate outputs
folder_path_code <- paste0(directory_path, "/Functions and General Analysis Code/") # change this to folder containing downloaded functions and code


source(paste0(folder_path_code,"Immunoprint functions.R"))

# Load data
combined_metadata_dist <- read_csv(paste0(folder_path,"Data Table S1.csv"))
all_zs_combined <- read_csv(paste0(folder_path,"Data Table S2.csv"))
all_hits_combined <- read_csv(paste0(folder_path,"Data Table S3.csv"))
Viral_Annotations <- read_csv(paste0(folder_path,"Data Table S4.csv"))

# Subset to Immunoprint development set
development_Clean_Mean_Zs <- subset_frame_fxn(all_zs_combined, combined_metadata_dist, "imm_dev")
development_metadata <- combined_metadata_dist %>% filter(imm_dev == TRUE)
development_Hits <- subset_frame_fxn(all_hits_combined, combined_metadata_dist, "imm_dev")

# Bounded Z scores
development_boundZs <- bound_Z(development_Clean_Mean_Zs)

# Public epitope list
public_epitopes <- Viral_Annotations$id[Viral_Annotations$public_epitope_2015 == TRUE]

# Correlation calculations
## Hits — Complete
Development_Cor <- stats::cor(development_Hits %>% select(-id), method = "pearson")
Development_Cor[upper.tri(Development_Cor)] <- 100
Development_Cor <- Development_Cor %>% as.data.frame()
Development_Cor$sample1 <- rownames(Development_Cor)

Development_Cor <- Development_Cor %>%
  pivot_longer(-sample1, names_to = "sample2", values_to = "HitCor_Complete") %>%
  filter(HitCor_Complete != 100, sample1 != sample2) %>%
  left_join(development_metadata %>% mutate(sample1 = identifier, indiv1 = donor) %>% select(sample1, indiv1)) %>%
  left_join(development_metadata %>% mutate(sample2 = identifier, indiv2 = donor) %>% select(sample2, indiv2)) %>%
  mutate(label = ifelse(indiv1 == indiv2, "Same individual", "Different individuals")) %>%
  select(sample1, indiv1, sample2, indiv2, label, HitCor_Complete)

## Hits — Public
temp <- stats::cor(
  development_Hits %>% filter(id %in% public_epitopes) %>% select(-id)
) %>% as.data.frame()
temp$sample1 <- rownames(temp)
temp <- temp %>% pivot_longer(-sample1, names_to = "sample2", values_to = "HitCor_Public")
Development_Cor <- left_join(Development_Cor, temp)

## Z scores — Complete
temp <- stats::cor(development_Clean_Mean_Zs %>% select(-id)) %>% as.data.frame()
temp$sample1 <- rownames(temp)
temp <- temp %>% pivot_longer(-sample1, names_to = "sample2", values_to = "ZCor_Complete")
Development_Cor <- left_join(Development_Cor, temp)

## Z scores — Public
temp <- stats::cor(
  development_Clean_Mean_Zs %>% filter(id %in% public_epitopes) %>% select(-id)
) %>% as.data.frame()
temp$sample1 <- rownames(temp)
temp <- temp %>% pivot_longer(-sample1, names_to = "sample2", values_to = "ZCor_Public")
Development_Cor <- left_join(Development_Cor, temp)

## Bounded Z — Complete
temp <- stats::cor(development_boundZs %>% select(-id)) %>% as.data.frame()
temp$sample1 <- rownames(temp)
temp <- temp %>% pivot_longer(-sample1, names_to = "sample2", values_to = "ZCorBound_Complete")
Development_Cor <- left_join(Development_Cor, temp)

## Bounded Z — Public
temp <- stats::cor(
  development_boundZs %>% filter(id %in% public_epitopes) %>% select(-id)
) %>% as.data.frame()
temp$sample1 <- rownames(temp)
temp <- temp %>% pivot_longer(-sample1, names_to = "sample2", values_to = "ZCorBound_Public")
Development_Cor <- left_join(Development_Cor, temp)


# Clean and annotate
Development_Cor <- Development_Cor %>%
  filter(sample1 %in% development_metadata$identifier,
         sample2 %in% development_metadata$identifier) %>%
  left_join(development_metadata %>% mutate(sample1 = identifier, Age_1 = Age_Years) %>% select(sample1, Age_1)) %>%
  left_join(development_metadata %>% mutate(sample2 = identifier, Age_2 = Age_Years) %>% select(sample2, Age_2)) %>%
  mutate(sampling_int = abs(Age_1 - Age_2)) %>%
  filter(sampling_int > 0)


# Performance evaluation
methods_vec <- colnames(
  Development_Cor %>%
    select(-c(sample1,sample2,indiv1,indiv2,label,Age_1,Age_2,sampling_int))
)

Development_Results <- data.frame(
  Method = methods_vec,
  OptimalThreshold = NA_real_,
  TP = NA_real_, FN = NA_real_, FP = NA_real_, TN = NA_real_,
  accuracy = NA_real_, sensitivity = NA_real_, specificity = NA_real_, AUC = NA_real_
)

for(i in seq_along(methods_vec)){
  tmet <- Development_Cor %>%
    select(all_of(methods_vec[i]), label)
  colnames(tmet) <- c("temp", "label")
  
  cp <- cutpointr(
    tmet,
    temp,
    label,
    method = maximize_metric,
    metric = accuracy
  )
  
  Development_Results[i, 2:10] <- c(
    summary(cp)[[7]] %>% unlist() %>% unname(),
    cp$acc, cp$sensitivity, cp$specificity, cp$AUC
  )
}

write.csv(
  Development_Results,
  file = paste0(intermediate_folder_path, "Development_Results.csv"),
  row.names = FALSE
)

# Final Immunoprint evaluation

immunoprint_Dev <- immunoprint(
  development_Clean_Mean_Zs,
  development_Clean_Mean_Zs,
  Viral_Annotations,
  keepAll = TRUE
)

immunoprint_Dev <- immunoprint_Dev %>%
  left_join(development_metadata %>% mutate(sample1 = identifier, indiv1 = donor) %>% select(sample1, indiv1)) %>%
  left_join(development_metadata %>% mutate(sample2 = identifier, indiv2 = donor) %>% select(sample2, indiv2)) %>%
  mutate(label = ifelse(indiv1 == indiv2, "Same individual", "Different individuals"))

cp <- cutpointr(
  immunoprint_Dev,
  immunoprint_match,
  label,
  method = maximize_metric,
  metric = sum_sens_spec
)

Development_Fxn_Results <- data.frame(
  Dataset = "Development",
  Method = "Final Immunoprint",
  TP = cp[[15]][[1]][2,2],
  FP = cp[[15]][[1]][2,4],
  TN = cp[[15]][[1]][2,5],
  FN = cp[[15]][[1]][2,3],
  sum_sens_spec = cp$sum_sens_spec,
  accuracy = cp$acc,
  sensitivity = cp$sensitivity,
  specificity = cp$specificity,
  AUC = cp$AUC,
  nSamples = length(unique(development_metadata$identifier)),
  nIndividuals = length(unique(development_metadata$donor))
)

write.csv(
  Development_Fxn_Results,
  file = paste0(intermediate_folder_path, "Development_Results_Final.csv"),
  row.names = FALSE
)
