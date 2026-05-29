# Load Packages
easypackages::libraries(c("readr","readxl","dplyr","tidyr","ggplot2","ggpubr"))

# insert path to folders
directory_path <- "~/HMS Dropbox/Moriah Mitchell/Elledge Lab/AVARS 2025/November_2025_Revision" # CHANGE THIS to path to director enclosing code and data sets
folder_path <- paste0(directory_path, "/Supplemental Tables/") # folder containing downloaded data sets
intermediate_folder_path <- paste0(directory_path, "/Intermediate Files Used in Plots/") # path containing intermediate outputs
folder_path_code <- paste0(directory_path, "/Functions and General Analysis Code/") # change this to folder containing downloaded functions and code
source(paste0(folder_path_code,"Immunoprint functions.R"))

# Load Supplemental Data
combined_metadata_dist <- read_csv(paste0(folder_path,"Data Table S1.csv"))
all_zs_combined <- read_csv(paste0(folder_path,"Data Table S2.csv"))
all_hits_combined <- read_csv(paste0(folder_path,"Data Table S3.csv"))
Viral_Annotations <- read_csv(paste0(folder_path,"Data Table S4.csv"))

# Pediatric
# Older Cohort
temp_donors<-combined_metadata_dist %>% filter(ped_long_2==TRUE) %>% select(donor) %>% distinct() %>% unlist() %>% unname()
temp_meta<-combined_metadata_dist %>% filter(donor %in% temp_donors)
# Arrange by age to identify consecutive pairs
temp_meta<-temp_meta %>% group_by(donor) %>% arrange(Age_Years) %>% mutate(sample_number=row_number()) %>% ungroup() %>% arrange(donor) %>% select(c(donor,Age_Years,identifier,sample_number))
# Identify consecutive pairs we want
t1<-temp_meta %>% mutate(Age_Years_1=Age_Years,identifier_1=identifier,sample_number_1=sample_number) %>% select(c(donor,Age_Years_1,identifier_1,sample_number_1))
t2<-temp_meta %>% mutate(Age_Years_2=Age_Years,identifier_2=identifier,sample_number_2=sample_number) %>% select(c(donor,Age_Years_2,identifier_2,sample_number_2))
ped_long_pairs<-t1 %>% left_join(t2,relationship = "many-to-many")
ped_long_pairs
# leftjoin Z scores onto this
temp_zs<-all_zs_combined %>% select(c(id,all_of(unique(c(ped_long_pairs$identifier_1,ped_long_pairs$identifier_2)))))
# pivot longer
temp_zs<-temp_zs %>% left_join(Viral_Annotations %>% select(c(id,species,family,family_label))) %>% pivot_longer(-c(id,species,family,family_label),names_to = "identifier",values_to = "ZScore")
ped_long_pairs<-ped_long_pairs %>% left_join(temp_zs %>% mutate(identifier_1=identifier,Z_1=ZScore) %>% select(-c(identifier,ZScore)))
ped_long_pairs<-ped_long_pairs %>% left_join(temp_zs %>% mutate(identifier_2=identifier,Z_2=ZScore) %>% select(-c(identifier,ZScore)))
# Add information
# Recognize is if hit
ped_long_pairs<-ped_long_pairs %>% mutate(Hit_1=ifelse(Z_1>=3.5,1,0),Hit_2=ifelse(Z_2>=3.5,1,0)) %>% mutate(percent_change=100*(Z_2-Z_1)/(Z_1+0.0001))
# Retain is stay recognize
ped_long_pairs<-ped_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==1,"Retained",NA))
# New is newly recognized
ped_long_pairs<-ped_long_pairs %>% mutate(change_cat=ifelse(Hit_1==0 & Hit_2==1,"New",change_cat))
# Lost is if lost
ped_long_pairs<-ped_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==0,"Lost",change_cat))
# Boost is >25% increase in zscore
ped_long_pairs<-ped_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==1 & percent_change>25,"Boosted",change_cat))
# Wane is >25% decrease in Z score
ped_long_pairs<-ped_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==1 & percent_change<(-25),"Waned",change_cat))
# Never recognized
ped_long_pairs<-ped_long_pairs %>% mutate(change_cat=ifelse(Hit_1==0 & Hit_2==0 ,"Not targeted",change_cat))

# Pairs meeting Figure 4 Jaccard criteria
temp_ped_long_pairs <- ped_long_pairs %>% mutate(time_interval=abs(Age_Years_2-Age_Years_1)) %>% filter(time_interval>=4,time_interval<=6) %>% select(-c(sample_number_1,sample_number_2,time_interval))
write.csv(temp_ped_long_pairs,file=paste0(intermediate_folder_path,"Ped_5_apart_Zs.csv"),row.names=FALSE)
ped_long_pairs<-ped_long_pairs %>% filter(sample_number_2==sample_number_1+1) %>% select(-c(sample_number_1,sample_number_2)) # just consecutive
write.csv(ped_long_pairs,file=paste0(intermediate_folder_path,"Ped_Consec_Zs.csv"),row.names=FALSE)
rm(temp_donors,temp_meta,t1,t2,temp_zs)

# Birth
temp_donors<-combined_metadata_dist %>% filter(ped_long_1==TRUE) %>% select(donor) %>% distinct() %>% unlist() %>% unname()
temp_meta<-combined_metadata_dist %>% filter(donor %in% temp_donors)
# Arrange by age to identify consecutive pairs
temp_meta<-temp_meta %>% group_by(donor) %>% arrange(Age_Years) %>% mutate(sample_number=row_number()) %>% ungroup() %>% arrange(donor) %>% select(c(donor,Age_Years,identifier,sample_number))
# Identify consecutive pairs we want
t1<-temp_meta %>% mutate(Age_Years_1=Age_Years,identifier_1=identifier,sample_number_1=sample_number) %>% select(c(donor,Age_Years_1,identifier_1,sample_number_1))
t2<-temp_meta %>% mutate(Age_Years_2=Age_Years,identifier_2=identifier,sample_number_2=sample_number) %>% select(c(donor,Age_Years_2,identifier_2,sample_number_2))
Birth_long_pairs<-t1 %>% left_join(t2,relationship = "many-to-many")
Birth_long_pairs
# leftjoin Z scores onto this
temp_zs<-all_zs_combined %>% select(c(id,all_of(unique(c(Birth_long_pairs$identifier_1,Birth_long_pairs$identifier_2)))))
# pivot longer
temp_zs<-temp_zs %>% left_join(Viral_Annotations %>% select(c(id,species,family,family_label))) %>% pivot_longer(-c(id,species,family,family_label),names_to = "identifier",values_to = "ZScore")
Birth_long_pairs<-Birth_long_pairs %>% left_join(temp_zs %>% mutate(identifier_1=identifier,Z_1=ZScore) %>% select(-c(identifier,ZScore)))
Birth_long_pairs<-Birth_long_pairs %>% left_join(temp_zs %>% mutate(identifier_2=identifier,Z_2=ZScore) %>% select(-c(identifier,ZScore)))
# Add information
# Recognize is if hit
Birth_long_pairs<-Birth_long_pairs %>% mutate(Hit_1=ifelse(Z_1>=3.5,1,0),Hit_2=ifelse(Z_2>=3.5,1,0)) %>% mutate(percent_change=100*(Z_2-Z_1)/(Z_1+0.0001))
# Retain is stay recognize
Birth_long_pairs<-Birth_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==1,"Retained",NA))
# New is newly recognized
Birth_long_pairs<-Birth_long_pairs %>% mutate(change_cat=ifelse(Hit_1==0 & Hit_2==1,"New",change_cat))
# Lost is if lost
Birth_long_pairs<-Birth_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==0,"Lost",change_cat))
# Boost is >25% increase in zscore
Birth_long_pairs<-Birth_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==1 & percent_change>25,"Boosted",change_cat))
# Wane is >25% decrease in Z score
Birth_long_pairs<-Birth_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==1 & percent_change<(-25),"Waned",change_cat))
# Never recognized
Birth_long_pairs<-Birth_long_pairs %>% mutate(change_cat=ifelse(Hit_1==0 & Hit_2==0 ,"Not targeted",change_cat))
# Save
Birth_long_pairs<-Birth_long_pairs %>% filter(sample_number_2==sample_number_1+1) %>% select(-c(sample_number_1,sample_number_2)) # just consecutive
write.csv(Birth_long_pairs,file=paste0(intermediate_folder_path,"Birth_Consec_Zs.csv"),row.names=FALSE)
rm(temp_donors,temp_meta,t1,t2,temp_zs)





# Adult
temp_donors<-combined_metadata_dist %>% filter(adult_long==TRUE)  %>% group_by(donor) %>% summarise(n_samps=max(row_number()),age_range=max(Age_Years)-min(Age_Years)) %>% filter(n_samps>=3,age_range>=5) %>% select(donor) %>% distinct() %>% unlist %>% unname()
temp_meta<-combined_metadata_dist %>% filter(donor %in% temp_donors)
# Arrange by age to identify consecutive pairs
temp_meta<-temp_meta %>% group_by(donor) %>% arrange(Age_Years) %>% mutate(sample_number=row_number()) %>% ungroup() %>% arrange(donor) %>% select(c(donor,Age_Years,identifier,sample_number))
# Identify consecutive pairs we want
t1<-temp_meta %>% mutate(Age_Years_1=Age_Years,identifier_1=identifier,sample_number_1=sample_number) %>% select(c(donor,Age_Years_1,identifier_1,sample_number_1))
t2<-temp_meta %>% mutate(Age_Years_2=Age_Years,identifier_2=identifier,sample_number_2=sample_number) %>% select(c(donor,Age_Years_2,identifier_2,sample_number_2))
adult_long_pairs<-t1 %>% left_join(t2,relationship = "many-to-many")
adult_long_pairs<-adult_long_pairs %>% filter(sample_number_2==sample_number_1+1) %>% select(-c(sample_number_1,sample_number_2))
adult_long_pairs
# left join Z scores onto this
temp_zs<-all_zs_combined %>% select(c(id,all_of(unique(c(adult_long_pairs$identifier_1,adult_long_pairs$identifier_2)))))
# pivot longer
temp_zs<-temp_zs %>% left_join(Viral_Annotations %>% select(c(id,species,family,family_label))) %>% pivot_longer(-c(id,species,family,family_label),names_to = "identifier",values_to = "ZScore")
adult_long_pairs<-adult_long_pairs %>% left_join(temp_zs %>% mutate(identifier_1=identifier,Z_1=ZScore) %>% select(-c(identifier,ZScore)))
adult_long_pairs<-adult_long_pairs %>% left_join(temp_zs %>% mutate(identifier_2=identifier,Z_2=ZScore) %>% select(-c(identifier,ZScore)))
# Add information
# Recognize is if hit
adult_long_pairs<-adult_long_pairs %>% mutate(Hit_1=ifelse(Z_1>=3.5,1,0),Hit_2=ifelse(Z_2>=3.5,1,0)) %>% mutate(percent_change=100*(Z_2-Z_1)/(Z_1+0.0001))
# Retain is stay recognize
adult_long_pairs<-adult_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==1,"Retained",NA))
# New is newly recognized
adult_long_pairs<-adult_long_pairs %>% mutate(change_cat=ifelse(Hit_1==0 & Hit_2==1,"New",change_cat))
# Lost is if lost
adult_long_pairs<-adult_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==0,"Lost",change_cat))
# Boost is >25% increase in zscore
adult_long_pairs<-adult_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==1 & percent_change>25,"Boosted",change_cat))
# Wane is >25% decrease in Z score
adult_long_pairs<-adult_long_pairs %>% mutate(change_cat=ifelse(Hit_1==1 & Hit_2==1 & percent_change<(-25),"Waned",change_cat))
# Never recognized
adult_long_pairs<-adult_long_pairs %>% mutate(change_cat=ifelse(Hit_1==0 & Hit_2==0 ,"Not targeted",change_cat))
write.csv(adult_long_pairs,file=paste0(intermediate_folder_path,"Adult_Consec_Zs.csv"),row.names=FALSE)
rm(temp_donors,temp_meta,t1,t2,temp_zs)
