# Load packages
easypackages::libraries(c("readr","readxl","dplyr","tidyr","ggplot2","ggpubr","data.table","tibble"))

# insert path to folders
# Define paths relative to the root of the cloned GitHub repository
directory_path <- here::here() sets
folder_path <- paste0(directory_path, "/Supplemental Tables/") # folder containing downloaded data sets
intermediate_folder_path <- paste0(directory_path, "/Intermediate Files Used in Plots/") # path containing intermediate outputs
folder_path_code <- paste0(directory_path, "/Functions and General Analysis Code/") # change this to folder

# Load functions
source(paste0(folder_path_code,"Immunoprint functions.R"))

# Load Supplemental Data
combined_metadata_dist <- read_csv(paste0(folder_path,"Data Table S1.csv"))
all_zs_combined <- read_csv(paste0(folder_path,"Data Table S2.csv"))
all_hits_combined <- read_csv(paste0(folder_path,"Data Table S3.csv"))
Viral_Annotations <- read_csv(paste0(folder_path,"Data Table S4.csv"))

# Adult data
## Filter to adult donors with samples from at least three time points spanning at least four years
temp_donors <- combined_metadata_dist %>% filter(adult_long==TRUE)  %>% group_by(donor) %>% summarise(n_samps=max(row_number()),age_range=max(Age_Years)-min(Age_Years)) %>% filter(n_samps>=3,age_range>=4)%>% select(donor) %>% distinct() %>% unlist %>% unname() 
length(unique(temp_donors))# 96 donors meet criteria
temp_meta<-combined_metadata_dist %>% filter(donor %in% temp_donors)
temp_zs<-all_zs_combined %>% select(c(id,temp_meta$identifier))
## generate trajectory 
Z_fam_cor_long_adult <- generate_trajectory_frame(temp_zs, Viral_Annotations, temp_meta)
Z_fam_cor_long_adult <- Z_fam_cor_long_adult %>% mutate(facet_label="Adults\n(96)")

# One year olds
## Filter to children between 0.8 and 1.5 years at initial sample collection with follow up extending to at least 4.5 years old
tt <- combined_metadata_dist %>% filter(ped_long_2==TRUE,Age_Years>=0.8,Age_Years<=1.5) %>% mutate(dist_from_1=abs(Age_Years-1)) %>% group_by(donor) %>% arrange(dist_from_1) %>% slice_head(n=1)%>% mutate(start_age=Age_Years,start_identifier=identifier) %>%  select(c(donor, start_age,start_identifier)) %>% left_join(combined_metadata_dist) %>% mutate(vec_from_start_age=Age_Years-start_age) %>% filter(vec_from_start_age>=0)
temp_donors<-tt %>% filter(vec_from_start_age>=4) %>% select(donor) %>% distinct() %>% unlist() %>% unname()
temp_meta<-combined_metadata_dist %>% filter(identifier %in% tt$identifier, donor %in% temp_donors)
length(unique(temp_meta$donor)) # 15 donors meet criteria
temp_zs<-all_zs_combined %>% select(c(id,temp_meta$identifier))
## generate trajectory 
Z_fam_cor_long_one <- generate_trajectory_frame(temp_zs, Viral_Annotations, temp_meta)
Z_fam_cor_long_one <- Z_fam_cor_long_one %>% mutate(facet_label="One-year-olds\n(15)")

# Four Year olds
## Filter to children between 3 and 5 years at initial sample collection with follow up extending at least 3 years
tt <- combined_metadata_dist %>% filter(ped_long_2==TRUE,Age_Years>=3,Age_Years<=5) %>% mutate(dist_from_4=abs(Age_Years-4)) %>% group_by(donor) %>% arrange(dist_from_4) %>% slice_head(n=1)%>% mutate(start_age=Age_Years,start_identifier=identifier) %>%  select(c(donor, start_age,start_identifier)) %>% left_join(combined_metadata_dist) %>% mutate(vec_from_start_age=Age_Years-start_age) %>% filter(vec_from_start_age>=0)
temp_donors<-tt %>% filter(vec_from_start_age>=4) %>% select(donor) %>% distinct() %>% unlist() %>% unname()
temp_meta<-combined_metadata_dist %>% filter(identifier %in% tt$identifier, donor %in% temp_donors)
length(unique(temp_meta$donor)) # 14 donors meet criteria
temp_zs<-all_zs_combined %>% select(c(id,temp_meta$identifier))
## generate trajectory 
Z_fam_cor_long_four <- generate_trajectory_frame(temp_zs, Viral_Annotations, temp_meta)
Z_fam_cor_long_four <- Z_fam_cor_long_four %>% mutate(facet_label="Four-year-olds\n(14)")

# Combine and save
Z_fam_cor_combined<-Z_fam_cor_long_four %>% bind_rows(Z_fam_cor_long_adult)%>% bind_rows(Z_fam_cor_long_one)
write.csv(Z_fam_cor_combined,file=paste0(intermediate_folder_path, "Viral_Family_Autocorrelation_Relative_to_initial_samples_by_age.csv"),row.names=FALSE)

# Downsampling analysis
temp_donors <- combined_metadata_dist %>% filter(adult_long==TRUE)  %>% group_by(donor) %>% summarise(n_samps=max(row_number()),age_range=max(Age_Years)-min(Age_Years)) %>% filter(n_samps>=3,age_range>=4)%>% select(donor) %>% distinct() %>% unlist %>% unname() 
length(unique(temp_donors))# 96 donors meet criteria
temp_meta<-combined_metadata_dist %>% filter(donor %in% temp_donors)
temp_zs<-all_zs_combined %>% select(c(id,temp_meta$identifier)) %>% left_join(Viral_Annotations %>% select(c(id, family_label)))

sample_sizes_list <- list(
  Herpesviridae    = c(1000, 2000, 5000, 7500, 10000, 15000, 17000),
  Orthomyxoviridae = c(1000, 2000, 5000, 7500),
  Picornaviridae   = c(1000, 2000, 5000)
)

downsampled_corr <- generate_downsampled_trajectory_frame(
  Zframe            = temp_zs,            # same format you pass to your original function
  viral_anno_frame  = Viral_Annotations,  # with id, family_label
  meta_frame        = temp_meta,        # with identifier, donor, Age_Years
  sample_sizes_list = sample_sizes_list,
  n_iter            = 200,
  seed              = 123
)
downsampled_corr <- downsampled_corr %>% filter(years_between!=0)
base_rows <- downsampled_corr %>%
  select(donor1, donor2, label,family_label,n_peptides_sampled,iter) %>%
  distinct() %>%
  mutate(
    fam_cor       = 1,
    years_between = 0,
    sample1       = NA_character_,
    sample2       = NA_character_,
    Age_Years1    = NA_real_,
    Age_Years2    = NA_real_,
    samp1         = "S1",
    samp2         = "S1"
  )
downsampled_corr <- bind_rows(downsampled_corr, base_rows)
write.csv(downsampled_corr,file=paste0(intermediate_folder_path, "Down_sampled_viral_family_correlation_adult.csv"),row.names=FALSE)
