# Load Packages
easypackages::libraries(c(
  "readr","readxl","dplyr","tidyr","ggplot2","ggpubr","cowplot","grid","ggplotify",
  "ComplexUpset","purrr","vegan","cluster","tibble"
))

# insert path to folders
directory_path <- "~/HMS Dropbox/Moriah Mitchell/Elledge Lab/AVARS 2025/November_2025_Revision" # CHANGE THIS to path to director enclosing code and data sets
folder_path <- paste0(directory_path, "/Supplemental Tables/") # folder containing downloaded data sets
intermediate_folder_path <- paste0(directory_path, "/Intermediate Files Used in Plots/") # path containing intermediate outputs

# Load Data
combined_metadata_dist <- readr::read_csv(paste0(folder_path,"Data Table S1.csv"))
all_hits_combined      <- readr::read_csv(paste0(folder_path,"Data Table S3.csv"))
Viral_Annotations      <- readr::read_csv(paste0(folder_path,"Data Table S4.csv"))
IgG_Data               <- readr::read_csv(paste0(folder_path,"Data Table S6.csv"))

# Custom functions
subset_frame_fxn <- function(wide_df, sample_metadata, dataset_lab) {
  sample_list <- sample_metadata$identifier[which(sample_metadata[[dataset_lab]] == TRUE)]
  wide_df %>% dplyr::select(c(id, all_of(sample_list)))
}

# Create long DF
adult_child_longzs <- all_hits_combined %>%
  dplyr::select(c(id, all_of(IgG_Data$identifier))) %>%
  tidyr::pivot_longer(-id, names_to = "identifier", values_to = "hit") %>%
  dplyr::left_join(IgG_Data %>% dplyr::select(c(identifier, donor, time_interval, Age_Years, age_group, interval_label)),
                   by = "identifier") %>%
  dplyr::mutate(age_group = ifelse(age_group == ">1yr", "1-3yr", age_group)) %>%
  dplyr::bind_rows(
    subset_frame_fxn(all_hits_combined, combined_metadata_dist, "adult_cross") %>%
      tidyr::pivot_longer(-id, names_to = "identifier", values_to = "hit") %>%
      dplyr::left_join(
        combined_metadata_dist %>%
          dplyr::mutate(time_interval = 40, age_group = "adult", interval_label = "Adult") %>%
          dplyr::select(c(identifier, donor, time_interval, Age_Years, age_group, interval_label)),
        by = "identifier"
      )
  )

# Add sample size to label
temp <- adult_child_longzs %>%
  dplyr::select(identifier, interval_label) %>%
  dplyr::distinct() %>%
  dplyr::group_by(interval_label) %>%
  dplyr::summarise(n_in_group = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(interval_label_updated = paste0(interval_label, " (", n_in_group, ")"))

adult_child_longzs <- adult_child_longzs %>%
  dplyr::left_join(temp %>% dplyr::select(-n_in_group), by = "interval_label")
rm(temp)

# Sample metadata (ONE ROW PER IDENTIFIER; matches main figure)
sample_meta <- adult_child_longzs %>%
  dplyr::group_by(donor, identifier, Age_Years, time_interval, interval_label, interval_label_updated) %>%
  dplyr::summarise(sum_hit = sum(hit), .groups = "drop") %>%
  dplyr::mutate(key = identifier)

# Build sample-by-peptide matrix (rows = identifier; matches main figure)
hit_unfiltered_matrix <- adult_child_longzs %>%
  dplyr::mutate(key = identifier) %>%
  tidyr::pivot_wider(id_cols = key, names_from = id, values_from = hit, values_fill = 0) %>%
  tibble::column_to_rownames(var = "key") %>%
  as.matrix()

# Match newer prevalence filter (>=5 and <= N-5)
cs <- colSums(hit_unfiltered_matrix)
keep <- cs >= 5 & cs <= (nrow(hit_unfiltered_matrix) - 5)
hit_unfiltered_matrix <- hit_unfiltered_matrix[, keep, drop = FALSE]
rm(cs, keep)

# Align sample order for distance + metadata
sample_meta <- sample_meta %>% dplyr::slice(match(rownames(hit_unfiltered_matrix), key))
stopifnot(all(sample_meta$key == rownames(hit_unfiltered_matrix)))

# Jaccard distances on binary matrix
D_jacc <- vegan::vegdist(hit_unfiltered_matrix, method = "jaccard", binary = TRUE)
write.csv(as.matrix(D_jacc),
          file = paste0(intermediate_folder_path, "Jaccard_Distance_Matrix.csv"),
          row.names = TRUE)

# PERMANOVA on Jaccard distances with permutations stratified by donor
adon_age <- vegan::adonis2(D_jacc ~ interval_label, data = sample_meta,
                           permutations = 999, strata = sample_meta$donor)
write.csv(adon_age, file = paste0(intermediate_folder_path, "PERMANOVA_donor_stratified.csv"),
          row.names = FALSE)

# Dispersion test (group spread)
bd <- vegan::betadisper(D_jacc, sample_meta$interval_label)
save(bd, file = paste0(intermediate_folder_path, "betadisper.rdata"))

# Label-permutation null (permute labels; keep distances fixed)
# Use permutations = 0 inside adonis2 since we generate the null externally
set.seed(10)
n_perm <- 500
perm_R2 <- replicate(n_perm, {
  perm <- sample(sample_meta$interval_label)
  a <- vegan::adonis2(D_jacc ~ perm, permutations = 0, strata = sample_meta$donor)
  a$R2[1]
})

write.csv(perm_R2, file = paste0(intermediate_folder_path, "label_permutation_null_R2.csv"),
          row.names = FALSE)