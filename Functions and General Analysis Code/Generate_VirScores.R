easypackages::libraries(c("readr","readxl","dplyr","tidyr","stringi","data.table"))
# Returns length of longest common
# substring of X[1..m] and Y[1..n]

LCSubStr<- function(X, Y) {
  m <- as.numeric(stringi::stri_length(X))
  # obtain lengths of strings being compared
  n <- as.numeric(stringi::stri_length(Y))
  LCSuff <- matrix(0, nrow=m+1, ncol=n+1) # create a zero matrix with m rows and n columns to store lengths of longest common suffixies of substrings
  # Note that LCSuff[i][j] contains the
  # length of longest common suffix of
  # X[1...i] and Y[1...j].
  result <- 0 # To store the length of the longest common substring
  # Following steps to build
  # LCSuff[m][n] in bottom up fashion
  for (i in 1:m+1) {
    for (j in 1:n+1) {
      if (i==1 | j==1) {
        LCSuff[i,j] <- 0
      } else if (substring(X,i-1,i-1)==substring(Y,j-1,j-1)) {
        LCSuff[i,j] <- LCSuff[i-1,j-1] + 1
        result <- max(result, LCSuff[i,j])
      }else {
        LCSuff[i,j] <- 0
      }
    }
  }
  return(result)
  # adapted from code by by Soumen Ghosh
}


subset_frame_fxn<-function(wide_df,sample_metadata,dataset_lab) {
  sample_list<-sample_metadata$identifier[which(sample_metadata[[dataset_lab]]==TRUE)]
  filtered_df<-wide_df %>% select(c(id,sample_list))
  return(filtered_df)
}


is_novel_peptide <- function(peptide, assigned_peptides, epitope_len) {
  for (i in seq_along(assigned_peptides)) {
    assigned_peptide <- assigned_peptides[i]
    match_size <- LCSubStr(peptide,assigned_peptide)
    if (match_size >= epitope_len) {
      return(FALSE)
    } 
  }
  return(TRUE)
} 



### Goal
"To bioinformatically remove cross-reactive antibodies, we first sorted the viruses by total
number of hits in descending order. We then iterated through each virus in this order. For
each virus, we iterated through each peptide hit. If the hit shared a subsequence of at least 7
aa with any hit previously observed in any of the viruses from that sample, that hit was
considered to be from a cross-reactive antibody and would be ignored for that virus.
Otherwise, the hit is considered to be specific and the score for that virus is incremented by
one. In this way, we summed only the peptide hits that do not share any linear epitopes. We
compared the final score for each virus to the threshold for that virus to determine whether
the sample is positive for exposure to that virus"
# Set file location
easypackages::libraries(c("dplyr","tidyr","data.table","readr"))
directory_path <- "~/HMS Dropbox/Moriah Mitchell/Elledge Lab/AVARS 2025/November_2025_Revision" # CHANGE THIS to path to director enclosing code and data sets
folder_path <- paste0(directory_path, "/Supplemental Tables/")
# import metadata
combined_metadata<- read_csv(paste0(folder_path,"Data Table S1.csv"))

#import hit file
hits_combined <- read_csv(paste0(folder_path,"Data Table S3.csv"))

hits_combined_temp <-subset_frame_fxn(hits_combined,combined_metadata,"adult_cross")
hits_combined <- subset_frame_fxn(hits_combined,combined_metadata,"ped_cross") %>% left_join(hits_combined_temp)
rm(hits_combined_temp)
#import epitope annotation file
Viral_Annotations <- read_csv(paste0(folder_path,"Data Table S4.csv"))
lib_annotation <- Viral_Annotations ##change this to path for library annotation

#not paralell; this takes a long time to generate
#compute hits per virus for each sample
hits_per_path <- hits_combined%>%left_join(lib_annotation%>%select(c(species,id)))%>%pivot_longer(-c(species,id),names_to = "sample_id",values_to = "hit")%>%group_by(sample_id,species)%>%summarise(nHits=sum(hit))%>%pivot_wider(id_cols = species,names_from = sample_id,values_from = nHits)

#Create data frame to store VirScores
virScores_combined<-hits_per_path%>%pivot_longer(-species,names_to = "sample_id",values_to = "VirScore")%>%mutate(VirScore=0)
#Create data frame to show nonredundant hits
non_redundant_hits<-hits_combined%>%pivot_longer(-id,names_to = "sample_id",values_to = "hit")%>%mutate(hit=0)
hits_combined<-hits_combined%>%left_join(lib_annotation%>%select(c(species, id, peptide))) #Add species and peptide sequence to hits_combined
samples<-colnames(hits_per_path%>%select(-species))
#for each sample, sort pathogens by hits per pathogen
for (i in seq_along(samples)) {
  samp<-samples[i]
  hits_per_path_ord<-hits_per_path%>%select(c(species,samp))
  hits_per_path_ord$nHits<-hits_per_path_ord[,2]
  hits_per_path_ord<-hits_per_path_ord%>%arrange(desc(nHits))%>%filter(nHits>0)
  hits_samp_i<-hits_combined%>%select(c(species,id,peptide,samp))
  hits_samp_i$hit<-hits_samp_i[,4]
  hits_samp_i<-hits_samp_i%>%filter(hit==TRUE)
  
  # We then iterated through each virus in this order. For each virus, we iterated through each peptide hit. 
  # If the hit shared a subsequence of at least 7aa with any hit previously observed in any of the viruses from
  # that sample, that hit was considered to be from a cross-reactive antibody and would be ignored for that virus.
  # Otherwise, the hit is considered to be specific and the score for that virus is incremented by one.
  for (j in seq_along(hits_per_path_ord$species)) {
    pathogen<-hits_per_path_ord$species[j]
    hits_sampi_pathj<-hits_samp_i%>%filter(species==pathogen) #get hits for just the pathogen and sample of interest
    for (k in seq_along(hits_sampi_pathj$id)) {
      if (k==1) {
        score=0
        assigned_peptides_j<-hits_sampi_pathj[1,]
        score=1
      } else if (is_novel_peptide(hits_sampi_pathj$peptide[k],assigned_peptides_j$peptide,7)==TRUE) {
        score=score+1
        assigned_peptides_j<-rbind(assigned_peptides_j,hits_sampi_pathj[k,])
      } 
    }
    virScores_combined<-virScores_combined%>%mutate(VirScore=ifelse(species==pathogen & sample_id==samp,score,VirScore)) #update the scores
    non_redundant_hits<-non_redundant_hits%>%mutate(hit=ifelse(sample_id==samp & id %in% assigned_peptides_j$id,1,hit))
  }
}
rm(assigned_peptides,score,pathogen,hits_per_path,hits_per_path_ord,hits_samp_i,hits_sampi_pathj)

AVARS_non_redundant_hits_wide<-non_redundant_hits %>% pivot_wider(id_cols = "id",names_from = "sample_id",values_from = "hit")
# Change to desired output location
write.csv(virScores_combined,file=paste0(folder_path,"/Data Table S5.csv"),row.names = FALSE)
