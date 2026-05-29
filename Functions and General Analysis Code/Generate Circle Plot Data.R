easypackages::libraries(c("readr","readxl","dplyr","tidyr","ggplot2","ggpubr","UpSetR","cowplot","grid","ggplotify"))

# insert path to folders
directory_path <- "~/HMS Dropbox/Moriah Mitchell/Elledge Lab/AVARS 2025/November_2025_Revision" # CHANGE THIS to path to director enclosing code and data sets
folder_path <- paste0(directory_path, "/Supplemental Tables/") # folder containing downloaded data sets
intermediate_folder_path <- paste0(directory_path, "/Intermediate Files Used in Plots/") # path containing intermediate outputs
folder_path_code <- paste0(directory_path, "/Functions and General Analysis Code/") # change this to folder containing downloaded functions and code

# Load Supplemental Data
combined_metadata_dist <- read_csv(paste0(folder_path,"Data Table S1.csv"))
all_zs_combined <- read_csv(paste0(folder_path,"Data Table S2.csv"))
all_hits_combined <- read_csv(paste0(folder_path,"Data Table S3.csv"))
Viral_Annotations <- read_csv(paste0(folder_path,"Data Table S4.csv"))
Birth_Consec_Zs <- read_csv(paste0(intermediate_folder_path,"Birth_Consec_Zs.csv")) # Generated in "Generate Consecutive Z Score Files.R"

# Load custom functions
source(paste0(folder_path_code,"Immunoprint functions.R"))

pal1<-c("#F87E2E", "#FFEC33", "#1FAA79", "#14C4F8", "#0077C8", "#9870F5", "#D53EE1", "#FF36B0", "#D72638", "#6E6E6E")


# Calculate exposures, grouping adenoviruses and Herpes 6A
# If these exposures co-occur, count as one exposure. If they do not, count as a second
# Determine New Exposures
#incorporate public epitopes
n_new_per_age_updated<-Birth_Consec_Zs %>% filter(species!="Hepacivirus C")%>% filter(grepl("C4",identifier_1)) %>% filter(change_cat=="New") %>% mutate(spec_label=ifelse(grepl("mastadenovirus",species),"Human adenovirus",species)) %>% mutate(spec_label=ifelse(species %in% c("Human betaherpesvirus 6A","Human betaherpesvirus 6B"),"Human betaherpesvirus 6",spec_label)) %>% mutate(spec_label=ifelse(species %in% c("Enterovirus A","Enterovirus B"),"Enterovirus A/B",spec_label))%>% mutate(spec_label=ifelse(species %in% c("Rhinovirus A","Rhinovirus B"),"Rhinovirus",spec_label))%>% mutate(spec_label=ifelse(species %in% c("Human alphaherpesvirus 1","Human alphaherpesvirus 2"),"HSV 1/2",spec_label)) %>% mutate(spec_label=ifelse(species %in% c("Influenza A virus","Influenza B virus"),"Influenza A/B",spec_label)) %>% mutate(spec_label=ifelse(family =="Papillomaviridae","Papillomavirus",spec_label))

n_new_per_age_updated<- n_new_per_age_updated%>% group_by(donor,Age_Years_1,Age_Years_2,family,family_label,spec_label,species) %>% summarise(n_new=max(row_number()))

max_n_new_updated<-n_new_per_age_updated %>% group_by(donor,species,family,family_label,spec_label) %>% summarise(max_n_new=max(n_new))
n_new_per_age_updated<-n_new_per_age_updated %>% left_join(max_n_new_updated)

age_at_exposure_updated<-n_new_per_age_updated %>% filter(n_new>=5) %>% filter(n_new>=0.5*max_n_new) # require the number of new peptides targeted to be at least 5 for and at least half the maximum new peptides ever observed for that virus

age_at_exposure_updated<-age_at_exposure_updated %>% left_join(Birth_Consec_Zs %>% filter(Hit_2==1)%>% select(c(species,family,donor,Age_Years_1,Age_Years_2,id))) %>% left_join(Viral_Annotations %>% select(c(id,public_epitope_2015))) %>% group_by(donor,Age_Years_1,Age_Years_2,species,family,family_label,n_new,max_n_new,spec_label) %>% summarise(n_pub=sum(public_epitope_2015==TRUE)) %>% filter(n_pub>0) #incorporate public epitopes

# select youngest age meeeting requirements
age_at_exposure_updated<-age_at_exposure_updated %>% group_by(donor,species,family,family_label,spec_label) %>% arrange(Age_Years_2) %>% slice_head(n=1)

# for the grouped viruses, only count one per timepoint
age_at_exposure_updated<-age_at_exposure_updated %>% group_by(donor,Age_Years_1,Age_Years_2,family,family_label,spec_label) %>% summarise(n_new=sum(n_new),n_pub=sum(n_pub))

# Add shape and color key
virus_key<-age_at_exposure_updated %>% group_by(spec_label,family,family_label) %>% summarise(n_exposed=max(row_number()))
virus_key<-virus_key %>% group_by(family_label) %>%mutate(shape_label=c(0,1,6,8,4,5,18,17,16)[row_number()])
family_cols <- c("Other"=pal1[1], "Herpesviridae"=pal1[9],"Picornaviridae"=pal1[3],"Orthomyxoviridae"=pal1[4], "Pneumoviridae"=pal1[5], "Adenoviridae"=pal1[6], "Coronaviridae"=pal1[7], "Paramyxoviridae"=pal1[8])
virus_key <- virus_key %>%
  mutate(colour = family_cols[family_label]) %>% arrange(family_label)
age_at_exposure_updated<-age_at_exposure_updated %>% left_join(virus_key %>% select(c(spec_label,family_label,shape_label,colour)))
# turn virus_key into two named vectors
species_shapes <- setNames(virus_key$shape_label, virus_key$spec_label)
species_cols   <- setNames(virus_key$colour,      virus_key$spec_label)
write.csv(virus_key, paste0(intermediate_folder_path,"virus_key_figure_5.csv"),row.names = FALSE)

# Create Plot Data
set.seed(10)
temp_donors<-age_at_exposure_updated %>% select(c(donor,Age_Years_2)) %>% distinct() %>% group_by(donor) %>% summarise(max_age=max(Age_Years_2)) %>% arrange(desc(max_age)) %>% filter(max_age>2.5)
# Distribute ages across interval
age_at_exposure_plot<-age_at_exposure_updated %>% filter(donor %in% temp_donors$donor) %>% group_by(Age_Years_1,Age_Years_2,donor) %>% mutate(age_mod=seq(from = Age_Years_1[1]+0.1, to = Age_Years_2[1],length.out=max(row_number())))
# Make measles and rubella occur at the same time
age_at_exposure_plot<-age_at_exposure_plot %>% mutate(age_mod=ifelse(spec_label %in% c("Rubivirus rubellae", "Measles morbillivirus"),(Age_Years_1+Age_Years_2)/2,age_mod)) %>% arrange(Age_Years_2) %>% mutate(ymod=runif(n=max(row_number()), min=3, max=6))

age_at_exposure_plot$spec_label <- factor(age_at_exposure_plot$spec_label,
                                          levels = virus_key$spec_label)

age_at_exposure_plot<-age_at_exposure_plot %>% mutate(age_mod=ifelse(spec_label %in% c("Rubivirus rubellae", "Measles morbillivirus"),(Age_Years_1+Age_Years_2)/2,age_mod)) %>% arrange(Age_Years_2)
counter=1
for (i in seq(age_at_exposure_plot$age_mod)) {
  if (counter==1) {
    age_at_exposure_plot$ymod[i] <- 4
    counter <- 2
  } else if (counter == 2) {
    age_at_exposure_plot$ymod[i] <- 5
    counter <- 3
  } else if (counter == 3) {
    age_at_exposure_plot$ymod[i] <- 6
    counter <- 4
  } else if (counter == 4) {
    age_at_exposure_plot$ymod[i]<- 5
    counter <-1
  }
}
write.csv(age_at_exposure_plot, paste0(intermediate_folder_path,"circle_plot_data_figure_5.csv"),row.names = FALSE)




