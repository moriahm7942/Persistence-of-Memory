easypackages::libraries(c(
  "readr","readxl","dplyr","tidyr","ggplot2","ggpubr","cowplot","grid","ggplotify",
  "uwot","ComplexUpset","purrr","vegan","cluster","tibble"
))

# insert path to folders
directory_path <- "~/HMS Dropbox/Moriah Mitchell/Elledge Lab/AVARS 2025/November_2025_Revision" # CHANGE THIS to path to director enclosing code and data sets
folder_path <- paste0(directory_path, "/Supplemental Tables/") # folder containing downloaded data sets
intermediate_folder_path <- paste0(directory_path, "/Intermediate Files Used in Plots/") # path containing intermediate outputs

# Load data
combined_metadata_dist <- readr::read_csv(paste0(folder_path, "Data Table S1.csv"))
all_hits_combined      <- readr::read_csv(paste0(folder_path, "Data Table S3.csv"))
Viral_Annotations      <- readr::read_csv(paste0(folder_path, "Data Table S4.csv"))
IgG_Data               <- readr::read_csv(paste0(folder_path, "Data Table S6.csv"))

# Custom functions
subset_frame_fxn <- function(wide_df, sample_metadata, dataset_lab) {
  sample_list <- sample_metadata$identifier[which(sample_metadata[[dataset_lab]] == TRUE)]
  wide_df %>% dplyr::select(c(id, all_of(sample_list)))
}

# Create long DF
adult_child_longzs <- all_hits_combined %>%
  dplyr::select(c(id, all_of(IgG_Data$identifier))) %>%
  tidyr::pivot_longer(-id, names_to = "identifier", values_to = "hit") %>%
  dplyr::left_join(
    IgG_Data %>% dplyr::select(c(identifier, donor, time_interval, Age_Years, age_group, interval_label)),
    by = "identifier"
  ) %>%
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

# Sample metadata (ONE ROW PER IDENTIFIER)
sample_meta <- adult_child_longzs %>%
  dplyr::group_by(donor, identifier, Age_Years, time_interval, interval_label, interval_label_updated) %>%
  dplyr::summarise(sum_hit = sum(hit), .groups = "drop") %>%
  dplyr::mutate(key = identifier)

# Build sample-by-peptide matrix (rows = identifier)
hit_unfiltered_matrix <- adult_child_longzs %>%
  dplyr::mutate(key = identifier) %>%
  tidyr::pivot_wider(id_cols = key, names_from = id, values_from = hit, values_fill = 0) %>%
  tibble::column_to_rownames(var = "key") %>%
  as.matrix()

# Feature filtering by prevalence 
cs <- colSums(hit_unfiltered_matrix)
keep <- cs >= 5 & cs <= (nrow(hit_unfiltered_matrix) - 5)
hit_unfiltered_matrix <- hit_unfiltered_matrix[, keep, drop = FALSE]
rm(cs, keep)

#  DISTANCES 
# Precompute Jaccard distances once and convert to a full matrix for uwot
dist_jaccard <- vegan::vegdist(hit_unfiltered_matrix, method = "jaccard", binary = TRUE)

D_jaccard_mat <- as.matrix(dist_jaccard)
# ensure row/col names exist and match samples
sample_keys <- rownames(hit_unfiltered_matrix)
rownames(D_jaccard_mat) <- sample_keys
colnames(D_jaccard_mat) <- sample_keys

# uwot expects a numeric matrix with finite values
storage.mode(D_jaccard_mat) <- "double"

knn_from_dist <- function(D, k) {
  n <- nrow(D)
  k <- min(k, n - 1)
  
  # For each row, take the k smallest distances excluding self
  idx <- t(vapply(seq_len(n), function(i) {
    ord <- order(D[i, ], na.last = NA)
    ord <- ord[ord != i]
    ord[seq_len(k)]
  }, integer(k)))
  
  dist <- t(vapply(seq_len(n), function(i) {
    D[i, idx[i, ]]
  }, numeric(k)))
  
  list(idx = idx, dist = dist)
}

# Run one UMAP with uwot
run_umap_once <- function(X,
                          n_neighbors = 15,
                          min_dist = 0.1,
                          metric = "jaccard",
                          seed = 1,
                          n_epochs = 200,
                          init = "spectral") {
  
  set.seed(seed)
  n <- nrow(X)
  if (n_neighbors >= n) n_neighbors <- n - 1
  
  if (metric == "jaccard") {
    # Exact Jaccard distances for binary matrix
    D <- as.matrix(vegan::vegdist(X, method = "jaccard", binary = TRUE))
    storage.mode(D) <- "double"
    
    # Precompute nearest neighbors from D (avoid Annoy "precomputed" bug)
    nn <- knn_from_dist(D, k = n_neighbors)
    
    emb <- uwot::umap(
      X = X,
      n_neighbors = n_neighbors,
      n_components = 2,
      nn_method = nn,          # <- precomputed neighbors (idx + dist)
      min_dist = min_dist,
      n_epochs = n_epochs,
      init = init,
      verbose = FALSE
    )
    
  } else {
    # These metrics are supported directly by your uwot build
    emb <- uwot::umap(
      X = X,
      n_neighbors = n_neighbors,
      n_components = 2,
      metric = metric,         # "cosine" or "euclidean"
      min_dist = min_dist,
      n_epochs = n_epochs,
      init = init,
      verbose = FALSE
    )
  }
  
  tibble::tibble(
    key = rownames(X),
    V1  = emb[, 1],
    V2  = emb[, 2],
    n_neighbors = n_neighbors,
    min_dist    = min_dist,
    metric      = metric,
    seed        = seed
  )
}

# Parameter grid
param_grid <- expand.grid(
  n_neighbors = c(5, 15, 30, 50),
  min_dist    = c(0.001, 0.1, 0.5),
  metric      = c("jaccard", "cosine", "euclidean"),
  seed        = 1:10,
  stringsAsFactors = FALSE
)

umap_runs <- purrr::pmap_dfr(
  param_grid,
  function(n_neighbors, min_dist, metric, seed) {
    run_umap_once(hit_unfiltered_matrix, n_neighbors, min_dist, metric, seed)
  }
) %>%
  dplyr::left_join(sample_meta %>% dplyr::select(key, interval_label_updated, interval_label), by = "key")

# Save
write.csv(umap_runs, file = paste0(intermediate_folder_path, "umap_runs.csv"), row.names = FALSE)