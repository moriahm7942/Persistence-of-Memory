easypackages::libraries("readr","dplyr","tidyr","data.table")

# insert path to folders
directory_path <- "~/HMS Dropbox/Moriah Mitchell/Elledge Lab/AVARS 2025/November_2025_Revision" # CHANGE THIS to path to director enclosing code and data sets
folder_path <- paste0(directory_path, "/Supplemental Tables/") # folder containing downloaded data sets
intermediate_folder_path <- paste0(directory_path, "/Intermediate Files Used in Plots/") # path containing intermediate outputs
folder_path_code <- paste0(directory_path, "/Functions and General Analysis Code/") # change this to folder containing downloaded functions and code

# Load required data
all_zs_combined <- read_csv(paste0(folder_path,"Data Table S2.csv"))
Viral_Annotations <- read_csv(paste0(folder_path,"Data Table S4.csv"))
source(paste0(folder_path_code,"Immunoprint functions.R"))

# All Zs
all_zs_bound<-bound_Z(all_zs_combined)
cor_complete<-stats::cor(all_zs_bound[,2:ncol(all_zs_bound)])
cor_complete[upper.tri(cor_complete,diag = TRUE)]<-2000
cor_complete<-cor_complete%>% as.data.frame()
cor_complete$sample1<-rownames(cor_complete)
cor_complete<-cor_complete %>% pivot_longer(-sample1,names_to = "sample2",values_to = "cor_complete") %>% filter(cor_complete!=2000)
write.csv(cor_complete,file = paste0(intermediate_folder_path, "complete_repertoire_all_samples_cor.csv"),row.names = FALSE)

# Public epitopes
pub_zs_bound <- all_zs_bound %>% filter(id %in% Viral_Annotations$id[Viral_Annotations$public_epitope_2015==TRUE])
cor_pub<-stats::cor(pub_zs_bound[,2:ncol(pub_zs_bound)])
cor_pub[upper.tri(cor_pub,diag = TRUE)]<-2000
cor_pub<-cor_pub%>% as.data.frame()
cor_pub$sample1<-rownames(cor_pub)
cor_pub<-cor_pub %>% pivot_longer(-sample1,names_to = "sample2",values_to = "cor_pub") %>% filter(cor_pub!=2000)
write.csv(cor_pub,file = paste0(intermediate_folder_path, "pub2015_repertoire_all_samples_cor.csv"),row.names = FALSE)
