easypackages::libraries(c("readr","readxl","dplyr","tidyr","ggplot2","ggpubr","UpSetR","cowplot","grid","ggplotify","data.table", "purrr"))

bound_fxn <- function(x,threshold=1,value=1, below=TRUE) {
  if (below==TRUE){
    out <-ifelse(x<threshold,value,x)
  } else {
    out<- ifelse(x>threshold,value,x)
  }
  return(out)
} ## bounds above if below==FALSE

hit_fxn <- function(x,threshold=3) {
  out <-ifelse(x<threshold,0,1)
  return(out)
} ## default hit if >= 3 

bound_Z <- function(Z_frame,lower=0,upper=15) {
  if (!("id" %in% colnames(Z_frame))) {
    stop("Data frame must contain id column")
  } else {
    Z_DT<-as.data.table(Z_frame)
    .cols <- setdiff(colnames(Z_frame),"id")
    Z_DT[, (.cols) := lapply(.SD,bound_fxn,threshold=lower,value=lower), .SDcols = .cols] ## bound below
    Z_DT[, (.cols) := lapply(.SD,bound_fxn,threshold=upper,value=upper,below=FALSE), .SDcols = .cols] ## bound above
    return(as.data.frame(Z_DT))
  }
}

Z_to_Hit <- function(Z_frame) {
  if (!("id" %in% colnames(Z_frame))) {
    stop("Data frame must contain id column")
  } else {
    Z_DT<-as.data.table(Z_frame)
    .cols <- setdiff(colnames(Z_frame),"id")
    Z_DT[, (.cols) := lapply(.SD,hit_fxn), .SDcols = .cols]
    return(as.data.frame(Z_DT))
  }
}

immunoprint <- function(Z_frame1,Z_frame2,prev,keepAll=FALSE) {
  if (!("id" %in% colnames(Z_frame1))) {
    stop("Data frame 1 must contain id column")} 
    else if (!("id" %in% colnames(Z_frame2))) {
      stop("Data frame 2 must contain id column") 
    } else if (!(("id" %in% colnames(prev))&("public_epitope_2015"%in% colnames(prev)))){
      stop("prev must contain id and public_epitope_2015 columns") 
    } else if(nrow(Z_frame1)!=nrow(Z_frame2)|nrow(Z_frame1)!=nrow(prev))
    {stop("All data frames must contain same number of rows")}
    else {
      zf1<-bound_Z(Z_frame1)
      zf2<-bound_Z(Z_frame2)
      zf_master<-left_join(zf1,zf2)
      # immprint score 1
      pub_eps<-prev$id[prev$public_epitope_2015==TRUE]
      out<-stats::cor(zf1%>%filter(id%in%pub_eps)%>%select(-id),zf2%>%filter(id%in%pub_eps)%>%select(-id),method=c("pearson"))%>%data.frame
      out$sample1<-rownames(out)
      out<-out%>%pivot_longer(-sample1,names_to = "sample2",values_to = "Score")%>%mutate(sorted=sprintf("%s:%s",pmin(sample2,sample2),pmax(sample1,sample2)))%>%distinct(sorted,.keep_all = TRUE)%>%filter(sample1!=sample2)%>%mutate(immunoprint_match=ifelse(Score>0.76,TRUE,FALSE))
      if(keepAll==FALSE){
        out<-out%>%filter(immunoprint_match==TRUE)
      } 
      return(out)
      }
}

subset_frame_fxn<-function(wide_df,sample_metadata,dataset_lab) {
  sample_list<-sample_metadata$identifier[which(sample_metadata[[dataset_lab]]==TRUE)]
  filtered_df<-wide_df %>% select(c(id,sample_list))
  return(filtered_df)
}

# Functions for trajectory analysis
# Function to generate data frame for viral family stratified heatmaps on longitudinal data
generate_heat_frame<-function(Zframe,viral_anno_frame, meta_frame) {
  long_heat_frame<-Zframe %>% mutate(id=as.double(id))%>% left_join(viral_anno_frame %>% select(c(id,species,`Protein names`,family_label))) %>% filter(family_label!="Other")%>% pivot_longer(-c(id,species,`Protein names`,family_label), names_to = "identifier",values_to = "ZScore") %>% mutate(adjZ=ifelse(ZScore<1,0,ZScore)) %>% mutate(adjZ=ifelse(ZScore>15,15,adjZ)) %>% left_join(meta_frame %>% select(c(identifier,donor,Age_Years,samp)) %>% mutate(samp_label=paste0("Age ",Age_Years)))
  return(long_heat_frame)
}

generate_trajectory_frame<-function(Zframe, viral_anno_frame, meta_frame) {
  # Bound Z scores
  Zframe_bound<-bound_Z(Zframe)
  # Alter meta
  meta_frame<-meta_frame %>% group_by(donor) %>% arrange(Age_Years) %>% mutate(samp=paste0("S",row_number()))
  ## Adenovirus
  Z_fam_cor<- data.frame(stats::cor(Zframe_bound %>% filter(id %in% filter(viral_anno_frame,family_label=="Adenoviridae")$id)%>%select(-id)))
  Z_fam_cor[upper.tri(Z_fam_cor)]<-100
  Z_fam_cor<-Z_fam_cor%>%data.frame
  Z_fam_cor$sample1<-rownames(Z_fam_cor) 
  Z_fam_cor<-Z_fam_cor%>%pivot_longer(-sample1,names_to = "sample2",values_to = "Adenoviridae") %>%filter(sample1!=sample2) %>% filter(Adenoviridae!=100)%>%left_join(meta_frame%>%mutate(sample2=identifier,donor2=donor,Age_Years2=Age_Years,samp2=samp)%>% ungroup %>% select(c(sample2,donor2,Age_Years2,samp2)))%>%left_join(meta_frame%>%mutate(sample1=identifier,donor1=donor,Age_Years1=Age_Years,samp1=samp)%>% ungroup %>% select(c(sample1,donor1,Age_Years1,samp1)))%>%mutate(label=ifelse(donor1==donor2,"Same individual","Different individuals"),years_between=abs(Age_Years1-Age_Years2)) # filter out repeat comparisons and join to metadata
  
  Z_fam_cor<-Z_fam_cor %>% filter(label=="Same individual",samp1=="S1"|samp2=="S1") # only interested in within individual and relative to baseline
  #
  ### Herpesviridae
  temp<-data.frame(stats::cor(Zframe_bound%>% filter(id %in% filter(viral_anno_frame,family_label=="Herpesviridae")$id)%>%select(-id)))
  temp$sample1<-rownames(temp)
  temp<-temp%>%pivot_longer(-sample1,names_to = "sample2",values_to = "Herpesviridae")
  Z_fam_cor<-left_join(Z_fam_cor,temp)
  
  ### Orthomyxoviridae
  temp<-data.frame(stats::cor(Zframe_bound%>% filter(id %in% filter(viral_anno_frame,family_label=="Orthomyxoviridae")$id)%>%select(-id)))
  temp$sample1<-rownames(temp)
  temp<-temp%>%pivot_longer(-sample1,names_to = "sample2",values_to = "Orthomyxoviridae")
  Z_fam_cor<-left_join(Z_fam_cor,temp)
  
  ### Picornaviridae
  temp<-data.frame(stats::cor(Zframe_bound%>% filter(id %in% filter(viral_anno_frame,family_label=="Picornaviridae")$id)%>%select(-id)))
  temp$sample1<-rownames(temp)
  temp<-temp%>%pivot_longer(-sample1,names_to = "sample2",values_to = "Picornaviridae")
  Z_fam_cor<-left_join(Z_fam_cor,temp)
  
  ### Pneumoviridae
  temp<-data.frame(stats::cor(Zframe_bound%>% filter(id %in% filter(viral_anno_frame,family_label=="Pneumoviridae")$id)%>%select(-id)))
  temp$sample1<-rownames(temp)
  temp<-temp%>%pivot_longer(-sample1,names_to = "sample2",values_to = "Pneumoviridae")
  Z_fam_cor<-left_join(Z_fam_cor,temp)
  
  ### Coronaviridae
  temp<-data.frame(stats::cor(Zframe_bound%>% filter(id %in% filter(viral_anno_frame,family_label=="Coronaviridae")$id)%>%select(-id)))
  temp$sample1<-rownames(temp)
  temp<-temp%>%pivot_longer(-sample1,names_to = "sample2",values_to = "Coronaviridae")
  Z_fam_cor<-left_join(Z_fam_cor,temp)
  
  ### Paramyxoviridae
  temp<-data.frame(stats::cor(Zframe_bound%>% filter(id %in% filter(viral_anno_frame,family_label=="Paramyxoviridae")$id)%>%select(-id)))
  temp$sample1<-rownames(temp)
  temp<-temp%>%pivot_longer(-sample1,names_to = "sample2",values_to = "Paramyxoviridae")
  Z_fam_cor<-left_join(Z_fam_cor,temp)
  
  # Artificially start at correlations of 1 for sampling interval=0
  temp<-Z_fam_cor %>% select(c(donor1,donor2,label)) %>% distinct() %>% mutate(Adenoviridae=1, Herpesviridae=1, Orthomyxoviridae=1, Picornaviridae=1, Pneumoviridae=1, Coronaviridae=1, Paramyxoviridae=1, years_between=0, Age_Years1=1, Age_Years2=1)
  Z_fam_cor_wide<-Z_fam_cor %>% bind_rows(temp) 
  return(Z_fam_cor_wide)
}


generate_trajectory_frame_single_family <- function(
    Zframe,            # full Zframe with id + sample cols
    viral_anno_frame,  # id -> family_label
    meta_frame,        # sample metadata
    family_name,       # e.g. "Herpesviridae"
    ids_subset = NULL  # optional: vector of peptide ids to restrict to
) {
  # Bound Z scores (your existing function)
  Zframe_bound <- bound_Z(Zframe)
  
  # Prepare metadata with samp = S1, S2, ...
  meta_prep <- meta_frame %>%
    group_by(donor) %>%
    arrange(Age_Years) %>%
    mutate(samp = paste0("S", row_number())) %>%
    ungroup()
  
  # IDs for this family
  fam_ids <- viral_anno_frame %>%
    filter(family_label == family_name) %>%
    pull(id) %>%
    unique()
  
  if (!is.null(ids_subset)) {
    fam_ids <- intersect(fam_ids, ids_subset)
  }
  
  if (length(fam_ids) < 2) {
    # Not enough peptides to compute a correlation
    return(tibble())
  }
  
  # Subset Zframe_bound to those peptides
  mat <- Zframe_bound %>%
    filter(id %in% fam_ids) %>%
    select(-id)
  
  if (nrow(mat) < 2) {
    return(tibble())
  }
  
  # Sample x sample correlation matrix for this family
  cor_mat <- stats::cor(mat)
  
  cor_long <- as.data.frame(cor_mat) %>%
    rownames_to_column("sample1") %>%
    pivot_longer(
      cols = -sample1,
      names_to = "sample2",
      values_to = "fam_cor"
    ) %>%
    filter(sample1 != sample2)
  
  # Join metadata to get donors, ages, etc.
  meta2 <- meta_prep %>%
    mutate(
      sample2   = identifier,
      donor2    = donor,
      Age_Years2 = Age_Years,
      samp2     = samp
    ) %>%
    select(sample2, donor2, Age_Years2, samp2)
  
  meta1 <- meta_prep %>%
    mutate(
      sample1   = identifier,
      donor1    = donor,
      Age_Years1 = Age_Years,
      samp1     = samp
    ) %>%
    select(sample1, donor1, Age_Years1, samp1)
  
  cor_long <- cor_long %>%
    left_join(meta2, by = "sample2") %>%
    left_join(meta1, by = "sample1") %>%
    mutate(
      label         = if_else(donor1 == donor2, "Same individual", "Different individuals"),
      years_between = abs(Age_Years1 - Age_Years2)
    ) %>%
    # Same filters as your original function
    filter(label == "Same individual", samp1 == "S1" | samp2 == "S1")
  
  # Add artificial time-zero rows with correlation = 1, like you do
  base_rows <- cor_long %>%
    select(donor1, donor2, label) %>%
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
  
  out <- bind_rows(cor_long, base_rows) %>%
    mutate(family_label = family_name)
  
  return(out)
}


generate_downsampled_trajectory_frame <- function(
    Zframe,
    viral_anno_frame,
    meta_frame,
    sample_sizes_list,
    n_iter = 200,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # Loop over families and their requested sample sizes
  res <- imap_dfr(sample_sizes_list, function(size_vec, fam_name) {
    
    # All peptide IDs for this family
    fam_ids_all <- viral_anno_frame %>%
      filter(family_label == fam_name) %>%
      pull(id) %>%
      unique()
    
    n_available <- length(fam_ids_all)
    
    if (n_available == 0) {
      warning("No peptides available for family ", fam_name, "; skipping.")
      return(tibble())
    }
    
    # For each requested sample size
    map_dfr(size_vec, function(n_sample) {
      
      if (n_sample > n_available) {
        warning(
          "Requested ", n_sample, " peptides for ", fam_name,
          " but only ", n_available, " available; skipping this size."
        )
        return(tibble())
      }
      
      # Repeat downsampling n_iter times
      map_dfr(seq_len(n_iter), function(iter_i) {
        
        sampled_ids <- sample(fam_ids_all, size = n_sample, replace = FALSE)
        
        traj_df <- generate_trajectory_frame_single_family(
          Zframe          = Zframe,
          viral_anno_frame = viral_anno_frame,
          meta_frame      = meta_frame,
          family_name     = fam_name,
          ids_subset      = sampled_ids
        )
        
        traj_df %>%
          mutate(
            n_peptides_sampled = n_sample,
            iter               = iter_i
          )
      })
    })
  })
  
  res
}


summarize_downsampled_by_model <- function(
    corr_df,
    time_points = c(0, 5, 10),
    df_spline = 3,
    return_errors = FALSE
) {
  # internal helper: always return same columns
  make_empty <- function(tp, err = NA_character_) {
    tibble(
      years_between = tp,
      emmean   = NA_real_,
      lower.CL = NA_real_,
      upper.CL = NA_real_,
      error_msg = err
    )
  }
  
  out <- corr_df %>%
    mutate(
      z_cor = atanh(pmin(pmax(fam_cor, -0.9999), 0.9999))
    ) %>%
    group_by(family_label, n_peptides_sampled, iter) %>%
    group_modify(~{
      dat <- .x
      
      # require enough observations and at least 2 donors for random effect
      if (nrow(dat) < 5 || dplyr::n_distinct(dat$donor1) < 2) {
        return(make_empty(time_points, err = "too_few_rows_or_donors"))
      }
      
      # choose df that is feasible given unique time values
      n_unique_t <- dplyr::n_distinct(dat$years_between)
      df_use <- min(df_spline, max(1, n_unique_t - 1))  # df <= (#unique - 1)
      
      form <- stats::as.formula(
        paste0("z_cor ~ splines::ns(years_between, df = ", df_use, ") + (1 | donor1)")
      )
      
      fit <- try(
        lme4::lmer(
          form,
          data = dat,
          control = lme4::lmerControl(optimizer = "bobyqa")
        ),
        silent = TRUE
      )
      
      if (inherits(fit, "try-error")) {
        msg <- conditionMessage(attr(fit, "condition"))
        return(make_empty(time_points, err = msg))
      }
      
      emm <- try(
        emmeans::emmeans(
          fit,
          ~ years_between,
          at = list(years_between = time_points)
        ),
        silent = TRUE
      )
      
      if (inherits(emm, "try-error")) {
        msg <- conditionMessage(attr(emm, "condition"))
        return(make_empty(time_points, err = msg))
      }
      
      tib <- as_tibble(emm) %>%
        transmute(
          years_between = years_between,
          emmean,
          lower.CL,
          upper.CL,
          error_msg = NA_character_
        )
      
      tib
    }) %>%
    ungroup() %>%
    mutate(
      est_r   = tanh(emmean),
      lower_r = tanh(lower.CL),
      upper_r = tanh(upper.CL)
    )
  
  if (!return_errors) {
    out <- out %>% select(-error_msg)
  }
  
  out
}
