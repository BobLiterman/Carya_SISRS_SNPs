# Function to calculate the modified Z-score for classification
mod_z_test <- function(statMatrix){
  prop_matrix <- prop.test(statMatrix)
  
  median_stat <- median(prop_matrix$estimate)
  mad_stat <- mad(prop_matrix$estimate,constant = 1)
  
  prop_matrix_median <- (prop_matrix$estimate - median_stat)/mad_stat
  
  diff_list <- c()
  p_list <- c()
  
  prop_matrix_df <- data.frame(statMatrix) %>% rownames_to_column("Companion_Sample")
  
  prop_matrix_df <- prop_matrix_df %>%
    mutate(Ratio = 100*(Grouping_Sites/(Grouping_Sites+Non_Grouping_Sites)),Z_Score = prop_matrix_median)
  return(prop_matrix_df)
}

# Pileup handling functions
remove_caret <- function(caret_string){
  if(!str_detect(caret_string,fixed("^"))){
    return(caret_string)
  } else{
    return(str_remove_all(caret_string,"\\^."))
  }
}
remove_dollar <- function(dollar_string){
  if(!str_detect(dollar_string,fixed("$"))){
    return(dollar_string)
  } else{
    return(str_remove_all(dollar_string,"\\$"))
  }
}
create_insertion_regex <- function(regex_insertion){
  digits <- as.integer(str_extract_all(regex_insertion,"([0-9]+)"))
  dots <- paste0(rep(".",digits),collapse = "")
  return(paste0("\\+",as.character(digits),dots,collapse = ''))
}
create_deletion_regex <- function(regex_deletion){
  digits <- as.integer(str_extract_all(regex_deletion,"([0-9]+)"))
  dots <- paste0(rep(".",digits),collapse = "")
  return(paste0("\\-",as.character(digits),dots,collapse = ''))
}
remove_insertions <- function(insertion_string){
  if(!str_detect(insertion_string,fixed("+"))){
    return(insertion_string)
  } else{
    insertion_list <- str_extract_all(insertion_string,"\\+([0-9]+)") %>% unlist()
    
    for(insert in insertion_list){
      insertion_string <- str_remove(insertion_string,create_insertion_regex(insert))
    }
    return(insertion_string)
  }
}
remove_deletions <- function(deletion_string){
  
  if(!str_detect(deletion_string,fixed("-"))){
    return(deletion_string)
  } else{
    deletion_list <- str_extract_all(deletion_string,"\\-([0-9]+)") %>% unlist()
    for(delete in deletion_list){
      deletion_string <- str_remove(deletion_string,create_deletion_regex(delete))
    }
    return(deletion_string)
  }
}
processPileups <- function(path_to_pileup){
  
  valid_pileup_chars <- char_class('A','a','C','c','G','g','T','t',fixed('*'))
  
  # Read in file, strip code column, and add Loc column
  my_pileup <- read_tsv(path_to_pileup,col_types = 'cicicc',col_names = c('Scaffold','Position','Ref_Base','Coverage','Sample_Bases','Code')) %>% 
    select(-Code) %>% 
    mutate(Loc = paste0(Scaffold,"/",Position)) 
  
  pileup_no_cov <- my_pileup %>% filter(Coverage == 0) 
  
  my_pileup <- my_pileup %>%
    filter(Coverage > 0) %>%
    rowwise() %>%
    mutate(Sample_Bases = remove_caret(Sample_Bases)) %>% # Remove ^* from Sample_Bases
    mutate(Sample_Bases = remove_dollar(Sample_Bases)) %>% # Remove $ from Sample_Bases
    mutate(Sample_Bases = remove_insertions(Sample_Bases)) %>% # Remove insertions from Sample_Bases
    mutate(Sample_Bases = remove_deletions(Sample_Bases)) %>% # Remove deletions from Sample_Bases
    mutate(Sample_Bases = str_replace_all(Sample_Bases,fixed("."),Ref_Base)) %>% # Replace periods with ref base 
    mutate(Sample_Bases = str_replace_all(Sample_Bases,fixed(","),Ref_Base)) %>% # Replace commas with ref base
    mutate(Sample_Bases = toupper(Sample_Bases)) %>% # Make bases upper case
    mutate(Allele_Count  = length(unique(unlist(str_extract_all(Sample_Bases,valid_pileup_chars))))) %>%
    ungroup()
  
  return_pileup <- my_pileup %>% filter(Allele_Count == 1) %>% rowwise() %>%
    mutate(Variation = "Homozygous",
           Sample_Base_1 = str_sub(Sample_Bases,start=1,end=1),
           Sample_Base_2 = NA) %>% ungroup() %>%
    rbind(my_pileup %>% filter(Allele_Count == 2) %>% rowwise() %>% 
            mutate(Variation = "Biallelic",
                   Sample_Base_1 = str_sub(paste0(naturalsort(unique(unlist(str_split(Sample_Bases,"")))),collapse = ""),start=1,end=1),
                   Sample_Base_2 = str_sub(paste0(naturalsort(unique(unlist(str_split(Sample_Bases,"")))),collapse = ""),start=2,end=2)) %>% ungroup()) %>%
    rbind(my_pileup %>% filter(Allele_Count > 2) %>% rowwise() %>% 
            mutate(Variation = "Polyallelic",
                   Sample_Base_1 = NA,
                   Sample_Base_2 = NA) %>% ungroup())
  
  if(nrow(pileup_no_cov) > 0){
    return_pileup <- return_pileup %>% rbind(pileup_no_cov %>% mutate(Sample_Bases = NA,Allele_Count = 0,Variation = "No_Coverage",Sample_Base_1 = NA,Sample_Base_2 = NA))
  }
  
  return(return_pileup)
}

# Load libraries
library(Rboretum)
library(rebus)
library(packcircles)
sourceRboretum()

#### Section 00: Set sample names, etc. ####
carya_name_df <- read_csv('data/Carya_Name_Table.csv')
signal_coltypes <- 'iccciciccccccc'

diploid_taxa <- c('CarAqu','CarCat','CarCor','CarIll','CarLac','CarMyr','CarOva','CarPalm')
study_samples <- c('CarAqu_1','CarAqu_2','CarAqu_3','CarCat_1','CarCat_2','CarCor_2','CarCor_3','CarCor_4','CarCor_5','CarCor_6','CarIll_1','CarIll_2','CarIll_3','CarIll_4','CarIll_5','CarLac_1','CarLac_3','CarMyr_1A','CarOva_1','CarOva_2','CarOva_3','CarPalm_1','CarPalm_2')
hybrid_samples <- c('myrxill_1','xbr_1','xbr_2','xbr_3','xbrl_1','xila_1','xila_2','xio_1','xlc_1','xlc_2','xlc_3','xlc_4','xnuss_1')

#### Section 01: Process singleton data ####

# The core set of 8 Carya species were analyzed in three ways. (1) Study samples only (2) Companion samples only (3) Pooled samples
# For each dataset, we extracted the singletons that defined each species. By default, these can include sites with more than one singleton position. 

# Read in signal/split data for each classifier sample

# <SAMPLE>_signal <- getAlignmentSignal(alignment_path = 'alignment_singletons_m0.phylip-relaxed-relaxed',use_gaps = TRUE)

study_singleton_signal <- read_tsv('data/Signal/Study_Singletons_Signal.tsv',col_types = signal_coltypes)
study_singleton_locs <- read_lines('data/Signal/Study_Singletons_LocList')

companion_singleton_signal <- read_tsv('data/Signal/Companion_Singletons_Signal.tsv',col_types = signal_coltypes)
companion_singleton_locs <- read_lines('data/Signal/Companion_Singletons_LocList')

pooled_singleton_signal <- read_tsv('data/Signal/Pooled_Singletons_Signal.tsv',col_types = signal_coltypes)
pooled_singleton_locs <- read_lines('data/Signal/Pooled_Singletons_LocList')

# Add loc information to signal
study_singletons_withLocs <- study_singleton_signal %>% mutate(Loc=study_singleton_locs)
companion_singletons_withLocs <- companion_singleton_signal %>% mutate(Loc=companion_singleton_locs)
pooled_singletons_withLocs <- pooled_singleton_signal %>% mutate(Loc=pooled_singleton_locs)

# Consider only 'true' singletons
study_singletons_withLocs <- study_singletons_withLocs %>% filter(Singleton_Count==1)
companion_singletons_withLocs <- companion_singletons_withLocs %>% filter(Singleton_Count==1)
pooled_singletons_withLocs <- pooled_singletons_withLocs %>% filter(Singleton_Count==1)

# Extract singleton locs
study_singletons <- study_singletons_withLocs$Loc
companion_singletons <- companion_singletons_withLocs$Loc
pooled_singletons <- pooled_singletons_withLocs$Loc

study_singletons_table <- table(study_singletons_withLocs$Singleton_Taxa)
companion_singletons_table <- table(companion_singletons_withLocs$Singleton_Taxa)
pooled_singletons_table <- table(pooled_singletons_withLocs$Singleton_Taxa)

# Get rough number comparison
# SAVED AS 'Carya_Singleton_Data.tsv'
singletons_df <- tibble('Species'=diploid_taxa,'Study'=as.integer(study_singletons_table),'Companion'=as.integer(companion_singletons_table),'Pooled'=as.integer(pooled_singletons_table)) %>%
  mutate(PercErr_Comp_v_Study = 100*(Companion-Study)/Study,
         PercErr_Pooled_v_Study = 100*(Pooled - Study)/Study,
         PercErr_Pooled_v_Comp = 100 *(Pooled-Companion)/Companion)



#### Section 02: Classify study samples with companion data ####

study_species_signal <- read_tsv('data/Signal/Classifier_Study_Samples.tsv',col_types = signal_coltypes)
study_species_locs <- read_lines('data/Signal/Classifier_Study_Samples_LocList')

# Get study sample data for companion samples with only 1 singleton
study_species_signal_withLocs <- study_species_signal %>% mutate(Loc = study_species_locs) %>% filter(Loc %in% companion_singletons)

study_species_classifier_df <- tibble(Study_Sample = character(),Companion_Sample = character(),Called_Sites = integer(),Grouping_Sites=integer())

# For each sample, find out how many sites had coverage, and how many sites match with each reference sample
for(study in study_samples){
  
  # Get relevant signal and get site counts
  study_species_df <- study_species_signal_withLocs %>% filter(Alignment_Name==study)
  study_counts <- nrow(study_species_df)
  
  # For each companion sample, tabulate support (How many sites support the grouping of <STUDY SAMPLE> and <COMPANION SAMPLE>)
  for(companion in diploid_taxa){
    clade <- naturalsort(c(study,companion)) %>% vectorSemi()
    
    # Get support for each grouping
    study_support <- getAlignmentSupport(signal = select(study_species_df,-Loc),clade = clade,dataset_name = study,return_integer = TRUE)
    
    # Add results to df
    study_species_classifier_df <- study_species_classifier_df %>% 
      add_row(Study_Sample = study,Companion_Sample = companion,Called_Sites = study_counts,Grouping_Sites = study_support)
  }
}

# Extract sites that group test samples with companion samples
study_species_grouping_df <- study_species_classifier_df %>%
  group_by(Study_Sample) %>%
  summarise(All_Grouping = sum(Grouping_Sites))

study_species_grouping_table <- study_species_grouping_df %>% 
  pull(All_Grouping) %>%
  `names<-`(study_samples)

study_species_classifier_summary <- study_species_classifier_df %>% 
  rowwise() %>% 
  mutate(Percent = Grouping_Sites/Called_Sites) %>%
  mutate(Non_Grouping_Sites = Called_Sites-Grouping_Sites) %>%
  ungroup() %>%
  select(Study_Sample,Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Percent) %>%
  arrange(Study_Sample,desc(Percent))

study_species_classifier_list <- list()

for(study in study_samples){
  study_matrix <- study_species_classifier_summary %>% filter(Study_Sample == study) %>%
    select(Companion_Sample,Grouping_Sites,Non_Grouping_Sites) %>%
    column_to_rownames(var = "Companion_Sample") %>%
    as.matrix()
  study_species_classifier_list[[study]] <- study_matrix
}

# Each study sample is compared to 8 reference samples, so alpha = 0.05/8 = 6.25E-3
# SAVED AS 'Species_Classifier.tsv'
study_species_stat_df <- purrr::map(.x=1:length(study_species_classifier_list),.f=function(x){mod_z_test(study_species_classifier_list[[x]]) %>% mutate(Study_Sample = names(study_species_classifier_list)[x])}) %>%
  do.call(rbind, .) %>%
  select(Study_Sample, Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Ratio,Z_Score) %>%
  mutate(P_Val = ifelse(Z_Score<0,pnorm(Z_Score),pnorm(-Z_Score))) %>%
  mutate(Significant_Bonferroni = ifelse(P_Val<0.00625 & Z_Score > 0,"Yes",'No')) %>%
  separate(col = Study_Sample,into=c('Species','Sample_Number'),remove = FALSE,sep = "_") %>%
  mutate(Match = ifelse(Species == Companion_Sample,TRUE,FALSE)) %>%
  arrange(Study_Sample,Grouping_Sites) %>%
  select(Study_Sample,Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Ratio,Z_Score,P_Val,Significant_Bonferroni,Match)

# Get top hits
study_species_top_hit_df <- study_species_stat_df %>% 
  group_by(Study_Sample) %>% 
  filter(Ratio == max(Ratio)) %>%
  ungroup() %>%
  arrange(desc(Ratio))

#### Section 03: Classify downsampled study samples with companion data ####

downsample_ids <- c('0.5X','0.25X','0.1X','0.05X','0.025X','0.01X')

# Read in downsampled signal files
ds_05X_signal <- read_tsv('data/Signal/Species_05X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.5X")
ds_025X_signal <- read_tsv('data/Signal/Species_025X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.25X")
ds_01X_signal <- read_tsv('data/Signal/Species_01X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.1X")
ds_005X_signal <- read_tsv('data/Signal/Species_005X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.05X")
ds_0025X_signal <- read_tsv('data/Signal/Species_0025X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.025X")
ds_001X_signal <- read_tsv('data/Signal/Species_001X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.01X")

downsample_study_species_df <- rbind(ds_05X_signal,ds_025X_signal,ds_01X_signal,ds_005X_signal,ds_0025X_signal,ds_001X_signal)
downsample_study_species_classifier_df <- tibble(Study_Sample = character(),Companion_Sample = character(),Called_Sites = integer(),Grouping_Sites=integer(),Downsampling=character())

for(study in study_samples){
  for(ds in downsample_ids){
    
    # Get relevant signal and get site counts
    study_sample_df <- downsample_study_species_df %>% filter(Alignment_Name==study,Downsampling==ds)
    study_counts <- nrow(study_sample_df)
    
    # For each companion sample, tabulate support (How many sites support the grouping of <STUDY SAMPLE> and <COMPANION SAMPLE>)
    for(companion in diploid_taxa){
      clade <- naturalsort(c(study,companion)) %>% vectorSemi()
      
      # Get support for each grouping
      study_support <- getAlignmentSupport(signal = select(study_sample_df,-Downsampling),clade = clade,dataset_name = study,return_integer = TRUE)
      
      # Add results to df
      downsample_study_species_classifier_df <- downsample_study_species_classifier_df %>% 
        add_row(Study_Sample = study,Companion_Sample = companion,Called_Sites = study_counts,Grouping_Sites = study_support,Downsampling = ds)
    }
  }
}

# Extract sites that group test samples with companion samples
downsample_study_species_grouping_df <- downsample_study_species_classifier_df %>%
  group_by(Study_Sample,Downsampling) %>%
  summarise(All_Grouping = sum(Grouping_Sites))

downsample_study_species_classifier_summary <- downsample_study_species_classifier_df %>% 
  rowwise() %>% 
  mutate(Percent = Grouping_Sites/Called_Sites) %>%
  mutate(Non_Grouping_Sites = Called_Sites-Grouping_Sites) %>%
  ungroup() %>%
  select(Study_Sample,Downsampling,Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Percent) %>%
  arrange(Study_Sample,desc(Percent)) %>%
  mutate(Full_Sample_Name = paste0(Study_Sample,"_",Downsampling))

all_downsampled_species <- pull(downsample_study_species_classifier_summary,Full_Sample_Name)

downsample_study_species_classifier_list <- list()

for(study in all_downsampled_species){
  study_matrix <- downsample_study_species_classifier_summary %>% filter(Full_Sample_Name == study) %>%
    select(Companion_Sample,Grouping_Sites,Non_Grouping_Sites) %>%
    column_to_rownames(var = "Companion_Sample") %>%
    as.matrix()
  downsample_study_species_classifier_list[[study]] <- study_matrix
}

# Each study sample is compared to 8 reference samples, so alpha = 0.05/8 = 6.25E-3
# SAVED AS 'Downsampled_Species_Classifier.tsv'
downsample_study_species_stat_df <- purrr::map(.x=1:length(downsample_study_species_classifier_list),.f=function(x){mod_z_test(downsample_study_species_classifier_list[[x]]) %>% mutate(Study_Sample = names(downsample_study_species_classifier_list)[x])}) %>%
  do.call(rbind, .) %>%
  select(Study_Sample,Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Ratio,Z_Score) %>%
  separate(col = Study_Sample,into=c('Species','Sample_Number','Downsampling'),remove = TRUE,sep = "_") %>%
  mutate(Study_Sample = paste0(Species,"_",Sample_Number)) %>%
  mutate(P_Val = ifelse(Z_Score<0,pnorm(Z_Score),pnorm(-Z_Score))) %>%
  mutate(Significant_Bonferroni = ifelse(P_Val<0.00625 & Z_Score > 0,"Yes",'No')) %>%
  mutate(Match = ifelse(Species == Companion_Sample,TRUE,FALSE)) %>%
  arrange(Study_Sample,Grouping_Sites) %>%
  select(Study_Sample,Species,Downsampling,Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Ratio,Z_Score,P_Val,Significant_Bonferroni,Match)

# Get top hits

downsample_base_table <- list('0.5X'=750000000/2,'0.25X'=750000000/4,'0.1X'=750000000/10,'0.05X'=750000000/20,'0.025X'=750000000/40,'0.01X'=750000000/100)

downsample_study_species_top_hit_df <- downsample_study_species_stat_df %>% 
  group_by(Study_Sample,Downsampling) %>% 
  filter(Ratio == max(Ratio)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(Bases = downsample_base_table[[Downsampling]],
         All_Sites = Grouping_Sites + Non_Grouping_Sites) %>%
  arrange(Study_Sample,Bases)

grouping_sites_lm <- lm(All_Sites ~ Bases, data=downsample_study_species_top_hit_df)
ratio_lm <- lm(Ratio ~ Bases, data=downsample_study_species_top_hit_df)
#### Section 04: Classify hybrid samples with pooled data ####

hybrid_pooled_signal <- read_tsv('data/Signal/Classifier_Hybrid_Pooled.tsv',col_types = signal_coltypes)
hybrid_pooled_locs <- read_lines('data/Signal/Classifier_Hybrid_Pooled_LocList')

hybrid_pooled_classifier_df <- tibble(Hybrid_Sample = character(),Companion_Sample = character(),Called_Sites = integer(),Grouping_Sites=integer())

for(hybrid in hybrid_samples){
  
  # Get relevant signal and get site counts
  hybrid_pooled_df <- hybrid_pooled_signal %>% mutate(Loc = hybrid_pooled_locs) %>% filter(Alignment_Name==hybrid,Loc %in% pooled_singletons)
  
  hybrid_pooled_counts <- nrow(hybrid_pooled_df)
  
  # For each companion sample, tabulate support (How many sites support the grouping of <STUDY SAMPLE> and <COMPANION SAMPLE>)
  for(companion in diploid_taxa){
    clade <- naturalsort(c(hybrid,companion)) %>% vectorSemi()
    
    # Get support for each grouping
    hybrid_pooled_support <- getAlignmentSupport(signal = select(hybrid_pooled_df,-Loc),clade = clade,dataset_name = hybrid,return_integer = TRUE)
    
    # Add results to df
    hybrid_pooled_classifier_df <- hybrid_pooled_classifier_df %>% 
      add_row(Hybrid_Sample = hybrid,Companion_Sample = companion,Called_Sites = hybrid_pooled_counts,Grouping_Sites = hybrid_pooled_support)
  }
}

# Extract sites that group test samples with reference samples
hybrid_pooled_grouping_df <- hybrid_pooled_classifier_df %>%
  group_by(Hybrid_Sample) %>%
  summarise(All_Grouping = sum(Grouping_Sites))

hybrid_pooled_grouping_table <- hybrid_pooled_grouping_df %>% 
  pull(All_Grouping) %>%
  `names<-`(hybrid_samples)

hybrid_pooled_classifier_summary <- hybrid_pooled_classifier_df %>% 
  rowwise() %>% 
  mutate(Percent = Grouping_Sites/hybrid_pooled_grouping_table[Hybrid_Sample]) %>%
  mutate(Non_Grouping_Sites = hybrid_pooled_grouping_table[Hybrid_Sample]-Grouping_Sites) %>%
  ungroup() %>%
  select(Hybrid_Sample,Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Percent) %>%
  arrange(Hybrid_Sample,desc(Percent))

hybrid_pooled_classifier_list <- list()

for(hybrid in hybrid_samples){
  
  hybrid_pooled_matrix <- hybrid_pooled_classifier_summary %>% filter(Hybrid_Sample == hybrid) %>%
    select(Companion_Sample,Grouping_Sites,Non_Grouping_Sites) %>%
    column_to_rownames(var = "Companion_Sample") %>%
    as.matrix()
  hybrid_pooled_classifier_list[[hybrid]] <- hybrid_pooled_matrix
}

# Each study sample is compared to 8 reference samples, so alpha = 0.05/8 = 6.25E-3
# SAVED AS 'Hybrid_Classifier.tsv'
hybrid_pooled_stat_df <- purrr::map(.x=1:length(hybrid_pooled_classifier_list),.f=function(x){mod_z_test(hybrid_pooled_classifier_list[[x]]) %>% mutate(Hybrid_Sample = names(hybrid_pooled_classifier_list)[x])}) %>%
  do.call(rbind, .) %>%
  select(Hybrid_Sample, Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Ratio,Z_Score) %>%
  mutate(P_Val = ifelse(Z_Score<0,pnorm(Z_Score),pnorm(-Z_Score))) %>%
  mutate(Significant_Bonferroni = ifelse(P_Val<0.00625 & Z_Score > 0,"Yes",'No')) %>%
  arrange(Hybrid_Sample,desc(Grouping_Sites))

#### Section 05: Classify downsampled hybrid samples with pooled data ####

downsample_ids <- c('0.5X','0.25X','0.1X','0.05X','0.025X','0.01X')

# Read in downsampled signal files
ds_05X_hybrid_signal <- read_tsv('data/Signal/Hybrids_05X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.5X")
ds_025X_hybrid_signal <- read_tsv('data/Signal/Hybrids_025X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.25X")
ds_01X_hybrid_signal <- read_tsv('data/Signal/Hybrids_01X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.1X")
ds_005X_hybrid_signal <- read_tsv('data/Signal/Hybrids_005X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.05X")
ds_0025X_hybrid_signal <- read_tsv('data/Signal/Hybrids_0025X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.025X")
ds_001X_hybrid_signal <- read_tsv('data/Signal/Hybrids_001X_Classifier.tsv',col_types = signal_coltypes) %>% mutate(Downsampling = "0.01X")

downsample_study_hybrid_df <- rbind(ds_05X_hybrid_signal,ds_025X_hybrid_signal,ds_01X_hybrid_signal,ds_005X_hybrid_signal,ds_0025X_hybrid_signal,ds_001X_hybrid_signal)
downsample_study_hybrid_classifier_df <- tibble(Study_Sample = character(),Companion_Sample = character(),Called_Sites = integer(),Grouping_Sites=integer(),Downsampling=character())

for(study in hybrid_samples){
  for(ds in downsample_ids){
    
    # Get relevant signal and get site counts
    study_sample_df <- downsample_study_hybrid_df %>% filter(Alignment_Name==study,Downsampling==ds)
    study_counts <- nrow(study_sample_df)
    
    # For each companion sample, tabulate support (How many sites support the grouping of <STUDY SAMPLE> and <COMPANION SAMPLE>)
    for(companion in diploid_taxa){
      clade <- naturalsort(c(study,companion)) %>% vectorSemi()
      
      # Get support for each grouping
      study_support <- getAlignmentSupport(signal = select(study_sample_df,-Downsampling),clade = clade,dataset_name = study,return_integer = TRUE)
      
      # Add results to df
      downsample_study_hybrid_classifier_df <- downsample_study_hybrid_classifier_df %>% 
        add_row(Study_Sample = study,Companion_Sample = companion,Called_Sites = study_counts,Grouping_Sites = study_support,Downsampling = ds)
    }
  }
}

# Extract sites that group test samples with companion samples
downsample_study_hybrid_grouping_df <- downsample_study_hybrid_classifier_df %>%
  group_by(Study_Sample,Downsampling) %>%
  summarise(All_Grouping = sum(Grouping_Sites))

downsample_study_hybrid_classifier_summary <- downsample_study_hybrid_classifier_df %>% 
  rowwise() %>% 
  mutate(Percent = Grouping_Sites/Called_Sites) %>%
  mutate(Non_Grouping_Sites = Called_Sites-Grouping_Sites) %>%
  ungroup() %>%
  select(Study_Sample,Downsampling,Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Percent) %>%
  arrange(Study_Sample,desc(Percent)) %>%
  mutate(Full_Sample_Name = paste0(Study_Sample,"_",Downsampling))

all_downsampled_hybrid <- pull(downsample_study_hybrid_classifier_summary,Full_Sample_Name)

downsample_study_hybrid_classifier_list <- list()

for(study in all_downsampled_hybrid){
  study_matrix <- downsample_study_hybrid_classifier_summary %>% filter(Full_Sample_Name == study) %>%
    select(Companion_Sample,Grouping_Sites,Non_Grouping_Sites) %>%
    column_to_rownames(var = "Companion_Sample") %>%
    as.matrix()
  downsample_study_hybrid_classifier_list[[study]] <- study_matrix
}

# Each study sample is compared to 8 reference samples, so alpha = 0.05/8 = 6.25E-3
# SAVED AS 'Downsampled_Hybrid_Classifier.tsv'
downsample_study_hybrid_stat_df <- purrr::map(.x=1:length(downsample_study_hybrid_classifier_list),.f=function(x){mod_z_test(downsample_study_hybrid_classifier_list[[x]]) %>% mutate(Study_Sample = names(downsample_study_hybrid_classifier_list)[x])}) %>%
  do.call(rbind, .) %>%
  select(Study_Sample,Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Ratio,Z_Score) %>%
  separate(col = Study_Sample,into=c('Species','Sample_Number','Downsampling'),remove = TRUE,sep = "_") %>%
  mutate(Study_Sample = paste0(Species,"_",Sample_Number)) %>%
  mutate(P_Val = ifelse(Z_Score<0,pnorm(Z_Score),pnorm(-Z_Score))) %>%
  mutate(Significant_Bonferroni = ifelse(P_Val<0.00625 & Z_Score > 0,"Yes",'No')) %>%
  arrange(Study_Sample,Grouping_Sites) %>%
  select(Study_Sample,Species,Downsampling,Companion_Sample,Grouping_Sites,Non_Grouping_Sites,Ratio,Z_Score,P_Val,Significant_Bonferroni)


#### Section 06: Alleles/Coverage ####

# Get singleton positions
caraqu_singletons <- pooled_singletons_withLocs %>% filter(Singleton_Taxa == "CarAqu") %>% pull(Loc)
carcor_singletons <- pooled_singletons_withLocs %>% filter(Singleton_Taxa == "CarCor") %>% pull(Loc)
carill_singletons <- pooled_singletons_withLocs %>% filter(Singleton_Taxa == "CarIll") %>% pull(Loc)
carmyr_singletons <- pooled_singletons_withLocs %>% filter(Singleton_Taxa == "CarMyr") %>% pull(Loc)
carova_singletons <- pooled_singletons_withLocs %>% filter(Singleton_Taxa == "CarOva") %>% pull(Loc)
carlac_singletons <- pooled_singletons_withLocs %>% filter(Singleton_Taxa == "CarLac") %>% pull(Loc)

# For screening xlc
ill_aqu_singletons <- c(carill_singletons,caraqu_singletons)

# For screening xbr
ill_cor_singletons <- c(carill_singletons,carcor_singletons)

# For screening xnuss
ill_lac_singletons <- c(carill_singletons,carlac_singletons)

# For screening myrxill
ill_myr_singletons <- c(carill_singletons,carmyr_singletons)

# For screening xio
ill_ova_singletons <- c(carill_singletons,carova_singletons)

# write_lines(ill_cor_singletons,'Manuscript_Output/CarIll_CarCor_Singletons.txt')
# write_lines(ill_aqu_singletons,'Manuscript_Output/CarIll_CarAqu_Singletons.txt')
# write_lines(ill_lac_singletons,'Manuscript_Output/CarIll_CarLac_Singletons.txt')
# write_lines(ill_myr_singletons,'Manuscript_Output/CarIll_CarMyr_Singletons.txt')
# write_lines(ill_ova_singletons,'Manuscript_Output/CarIll_CarOva_Singletons.txt')

# Use these singletons to extract all SNPs from both parents, and each hybrid pileup file

xlc_caraqu_pileup <- processPileups('data/Focal/xlc/Focal_xlc_CarAqu.pileups') %>% mutate(Sample = "xlc_CarAqu")
xlc_carill_pileup <- processPileups('data/Focal/xlc/Focal_xlc_CarIll.pileups') %>% mutate(Sample = "xlc_CarIll")

xlc1_pileup <- processPileups('data/Focal/xlc/Focal_xlc_1.pileups') %>% mutate(Sample = "xlc_1")
xlc2_pileup <- processPileups('data/Focal/xlc/Focal_xlc_2.pileups') %>% mutate(Sample = "xlc_2")
xlc3_pileup <- processPileups('data/Focal/xlc/Focal_xlc_3.pileups') %>% mutate(Sample = "xlc_3")
xlc4_pileup <- processPileups('data/Focal/xlc/Focal_xlc_4.pileups') %>% mutate(Sample = "xlc_4")

xbr_carcor_pileup <- processPileups('data/Focal/xbr/Focal_xbr_CarCor.pileups') %>% mutate(Sample = "xbr_CarCor")
xbr_carill_pileup <- processPileups('data/Focal/xbr/Focal_xbr_CarIll.pileups') %>% mutate(Sample = "xbr_CarIll")

xbr1_pileup <- processPileups('data/Focal/xbr/Focal_xbr_1.pileups') %>% mutate(Sample = "xbr_1")
xbr2_pileup <- processPileups('data/Focal/xbr/Focal_xbr_2.pileups') %>% mutate(Sample = "xbr_2")
xbr3_pileup <- processPileups('data/Focal/xbr/Focal_xbr_3.pileups') %>% mutate(Sample = "xbr_3")

xio_carova_pileup <- processPileups('data/Focal/xio/Focal_xio_CarOva.pileups') %>% mutate(Sample = "xio_CarOva")
xio_carill_pileup <- processPileups('data/Focal/xio/Focal_xio_CarIll.pileups') %>% mutate(Sample = "xio_CarIll")

xio_pileup <- processPileups('data/Focal/xio/Focal_xio.pileups') %>% mutate(Sample = "xio")

xnuss_carlac_pileup <- processPileups('data/Focal/xnuss/Focal_xnuss_CarLac.pileups') %>% mutate(Sample = "xnuss_CarLac")
xnuss_carill_pileup <- processPileups('data/Focal/xnuss/Focal_xnuss_CarIll.pileups') %>% mutate(Sample = "xnuss_CarIll")

xnuss_pileup <- processPileups('data/Focal/xnuss/Focal_xnuss.pileups') %>% mutate(Sample = "xnuss")

myrxill_carmyr_pileup <- processPileups('data/Focal/myrxill/Focal_myrxill_CarMyr.pileups') %>% mutate(Sample = "myrxill_CarMyr")
myrxill_carill_pileup <- processPileups('data/Focal/myrxill/Focal_myrxill_CarIll.pileups') %>% mutate(Sample = "myrxill_CarIll")

myrxill_pileup <- processPileups('data/Focal/myrxill/Focal_myrxill.pileups') %>% mutate(Sample = "myrxill")

# For each hybrid, add in rows where no coverage was present in the pileup

xlc1_pileup <- xlc1_pileup %>% rbind(xlc_carill_pileup %>% filter(!Loc %in% xlc1_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                        Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "xlc_1"))
xlc2_pileup <- xlc2_pileup %>% rbind(xlc_carill_pileup %>% filter(!Loc %in% xlc2_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                        Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "xlc_2"))
xlc3_pileup <- xlc3_pileup %>% rbind(xlc_carill_pileup %>% filter(!Loc %in% xlc3_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                        Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "xlc_3"))
xlc4_pileup <- xlc4_pileup %>% rbind(xlc_carill_pileup %>% filter(!Loc %in% xlc4_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                        Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "xlc_4"))

xbr1_pileup <- xbr1_pileup %>% rbind(xbr_carill_pileup %>% filter(!Loc %in% xbr1_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                        Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "xbr_1"))
xbr2_pileup <- xbr2_pileup %>% rbind(xbr_carill_pileup %>% filter(!Loc %in% xbr2_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                        Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "xbr_2"))
xbr3_pileup <- xbr3_pileup %>% rbind(xbr_carill_pileup %>% filter(!Loc %in% xbr3_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                        Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "xbr_3"))

xio_pileup <- xio_pileup %>% rbind(xio_carill_pileup %>% filter(!Loc %in% xio_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                     Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "xio"))

xnuss_pileup <- xnuss_pileup %>% rbind(xnuss_carill_pileup %>% filter(!Loc %in% xnuss_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                             Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "xnuss"))

myrxill_pileup <- myrxill_pileup %>% rbind(myrxill_carill_pileup %>% filter(!Loc %in% myrxill_pileup$Loc) %>% mutate(Coverage = 0,Allele_Count = 0,Variation = "No_Coverage",
                                                                                                                     Sample_Base_1 = NA,Sample_Base_2 = NA,Sample = "myrxill"))

# Make dataframes for each cross
xbr_df <- xbr_carcor_pileup %>% select(Loc,Sample_Bases) %>%
  rowwise() %>%
  mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
  ungroup() %>%
  rename(CarCor_Base = 'Sample_Bases') %>%
  left_join(xbr_carill_pileup %>% select(Loc,Sample_Bases) %>%
              rowwise() %>%
              mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
              ungroup() %>%
              rename(CarIll_Base = 'Sample_Bases'))
xbr1_df <- xbr_df %>%
  left_join(xbr1_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc')

xbr2_df <- xbr_df %>%
  left_join(xbr2_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc')

xbr3_df <- xbr_df %>%
  left_join(xbr3_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc')

xbr_df <- rbind(xbr1_df,xbr2_df,xbr3_df) %>% rowwise() %>%
  mutate(Base_1_Match = ifelse(is.na(Sample_Base_1),'No_Coverage',ifelse(Sample_Base_1 == CarCor_Base,'CarCor',ifelse(Sample_Base_1 == CarIll_Base,'CarIll','Other'))),
         Base_2_Match = ifelse(is.na(Sample_Base_2),'No_Coverage',ifelse(Sample_Base_2 == CarCor_Base,'CarCor',ifelse(Sample_Base_2 == CarIll_Base,'CarIll','Other')))) %>%
  select(Sample,Loc,Variation,Coverage,Base_1_Match,Base_2_Match)

xlc_df <- xlc_caraqu_pileup %>% select(Loc,Sample_Bases) %>%
  rowwise() %>%
  mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
  ungroup() %>%
  rename(CarAqu_Base = 'Sample_Bases') %>%
  left_join(xlc_carill_pileup %>% select(Loc,Sample_Bases) %>%
              rowwise() %>%
              mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
              ungroup() %>%
              rename(CarIll_Base = 'Sample_Bases'))

xlc1_df <- xlc_df %>%
  left_join(xlc1_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc')

xlc2_df <- xlc_df %>%
  left_join(xlc2_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc')

xlc3_df <- xlc_df %>%
  left_join(xlc3_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc')

xlc4_df <- xlc_df %>%
  left_join(xlc4_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc')

xlc_df <- rbind(xlc1_df,xlc2_df,xlc3_df,xlc4_df)  %>%
  mutate(Base_1_Match = ifelse(is.na(Sample_Base_1),'No_Coverage',ifelse(Sample_Base_1 == CarAqu_Base,'CarAqu',ifelse(Sample_Base_1 == CarIll_Base,'CarIll','Other'))),
         Base_2_Match = ifelse(is.na(Sample_Base_2),'No_Coverage',ifelse(Sample_Base_2 == CarAqu_Base,'CarAqu',ifelse(Sample_Base_2 == CarIll_Base,'CarIll','Other')))) %>%
  select(Sample,Loc,Variation,Coverage,Base_1_Match,Base_2_Match)

xio_df <- xio_carova_pileup %>% select(Loc,Sample_Bases) %>%
  rowwise() %>%
  mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
  ungroup() %>%
  rename(CarOva_Base = 'Sample_Bases') %>%
  left_join(xio_carill_pileup %>% select(Loc,Sample_Bases) %>%
              rowwise() %>%
              mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
              ungroup() %>%
              rename(CarIll_Base = 'Sample_Bases')) %>%
  left_join(xio_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc') %>%
  mutate(Base_1_Match = ifelse(is.na(Sample_Base_1),'No_Coverage',ifelse(Sample_Base_1 == CarOva_Base,'CarOva',ifelse(Sample_Base_1 == CarIll_Base,'CarIll','Other'))),
         Base_2_Match = ifelse(is.na(Sample_Base_2),'No_Coverage',ifelse(Sample_Base_2 == CarOva_Base,'CarOva',ifelse(Sample_Base_2 == CarIll_Base,'CarIll','Other')))) %>%
  select(Sample,Loc,Variation,Coverage,Base_1_Match,Base_2_Match)


xnuss_df <- xnuss_carlac_pileup %>% select(Loc,Sample_Bases) %>%
  rowwise() %>%
  mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
  ungroup() %>%
  rename(CarLac_Base = 'Sample_Bases') %>%
  left_join(xnuss_carill_pileup %>% select(Loc,Sample_Bases) %>%
              rowwise() %>%
              mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
              ungroup() %>%
              rename(CarIll_Base = 'Sample_Bases')) %>%
  left_join(xnuss_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc') %>%
  mutate(Base_1_Match = ifelse(is.na(Sample_Base_1),'No_Coverage',ifelse(Sample_Base_1 == CarLac_Base,'CarLac',ifelse(Sample_Base_1 == CarIll_Base,'CarIll','Other'))),
         Base_2_Match = ifelse(is.na(Sample_Base_2),'No_Coverage',ifelse(Sample_Base_2 == CarLac_Base,'CarLac',ifelse(Sample_Base_2 == CarIll_Base,'CarIll','Other')))) %>%
  select(Sample,Loc,Variation,Coverage,Base_1_Match,Base_2_Match)


myrxill_df <- myrxill_carmyr_pileup %>% select(Loc,Sample_Bases) %>%
  rowwise() %>%
  mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
  ungroup() %>%
  rename(CarMyr_Base = 'Sample_Bases') %>%
  left_join(myrxill_carill_pileup %>% select(Loc,Sample_Bases) %>%
              rowwise() %>%
              mutate(Sample_Bases = str_sub(Sample_Bases,start=1,end=1)) %>%
              ungroup() %>%
              rename(CarIll_Base = 'Sample_Bases')) %>%
  left_join(myrxill_pileup %>% select(Loc,Sample,Variation,Coverage,Allele_Count,Sample_Base_1,Sample_Base_2),by='Loc') %>%
  mutate(Base_1_Match = ifelse(is.na(Sample_Base_1),'No_Coverage',ifelse(Sample_Base_1 == CarMyr_Base,'CarMyr',ifelse(Sample_Base_1 == CarIll_Base,'CarIll','Other'))),
         Base_2_Match = ifelse(is.na(Sample_Base_2),'No_Coverage',ifelse(Sample_Base_2 == CarMyr_Base,'CarMyr',ifelse(Sample_Base_2 == CarIll_Base,'CarIll','Other')))) %>%
  select(Sample,Loc,Variation,Coverage,Base_1_Match,Base_2_Match)



allele_df <- rbind(xbr_df,xlc_df,xio_df,xnuss_df,myrxill_df) %>% ungroup()

allele_nocov_df <- allele_df %>% filter(Coverage == 0)

allele_cov_df <- allele_df %>% filter(Coverage > 0)

# Get total sites with coverage
coverage_site_df <- allele_cov_df %>% group_by(Sample,Coverage) %>% 
  summarize(Coverage_Sites=n()) %>% 
  ungroup() %>% rowwise() %>%
  mutate(Name = paste(c(Sample,Coverage),collapse = "_")) %>%
  ungroup() %>%
  select(Name,Coverage_Sites)

coverage_site_table <- as.integer(coverage_site_df$Coverage_Sites) %>% `names<-`(coverage_site_df$Name)

# Get total sites by sample
total_site_df <- allele_cov_df %>% group_by(Sample) %>% summarize(Total_Sites = n())
total_site_table <-  as.integer(total_site_df$Total_Sites) %>% `names<-`(total_site_df$Sample)

# Get total sites by pattern
variation_site_df <- allele_cov_df %>% group_by(Sample,Variation) %>% summarize(Variation_Sites = n()) %>% rowwise() %>%
  mutate(variation_Name = paste(c(Sample,Variation),collapse = "_")) %>% ungroup()
variation_site_table <-  as.integer(variation_site_df$Variation_Sites) %>% `names<-`(variation_site_df$variation_Name)

# Get coverage statistics
coverage_df <- allele_cov_df %>% group_by(Sample,Coverage,Variation) %>% 
  summarize(Pattern_Coverage_Sites = n()) %>% 
  ungroup() %>% rowwise() %>%
  mutate(Name = paste(c(Sample,Coverage),collapse = "_"),
         Variation_Name = paste(c(Sample,Variation),collapse = "_")) %>%
  ungroup() %>% 
  left_join(coverage_site_df,by='Name') %>% select(-Name) %>%
  rowwise() %>%
  mutate(Total_Sites = total_site_table[[Sample]],
         Variation_Sites = variation_site_table[[Variation_Name]]) %>% ungroup() %>%
  select(Sample,Total_Sites,Variation,Variation_Sites,Coverage,Coverage_Sites,Pattern_Coverage_Sites)

coverage_breakdown <- coverage_df %>% mutate(Percent_Coverage = Coverage_Sites/Total_Sites) %>%
  group_by(Sample) %>% filter(!duplicated(Coverage)) %>% ungroup() %>%
  select(Sample,Coverage,Coverage_Sites,Percent_Coverage)

coverage_breakdown[is.na(coverage_breakdown)] <- 0

# SAVED AS 'Hybrid_Coverage_Breakdown.tsv'
pattern_breakdown <- coverage_df %>% mutate(Percent_Variation = Variation_Sites/Total_Sites) %>%
  group_by(Sample) %>% filter(!duplicated(Variation)) %>% ungroup() %>%
  select(Sample,Total_Sites,Variation,Variation_Sites,Percent_Variation) %>%
  pivot_wider(names_from = 'Variation',values_from = c('Variation_Sites','Percent_Variation')) %>%
  `names<-`(c('Sample','Total_Sites',
              'Homozygous_Sites','Biallelic_Sites','Polyallelic_Sites',
              'Percent_Homozygous','Percent_Biallelic','Percent_Polyallelic'))

# Get odds of species match at homozygous sites
# SAVED AS 'Hybrid_Homozygous_Breakdown.tsv'
homozygous_df <- allele_cov_df %>% 
  filter(Variation == 'Homozygous') %>% 
  mutate(Base_1_Match = ifelse(Base_1_Match %in% c('CarCor','CarAqu','CarLac','CarOva','CarMyr'),'Other_Parent',Base_1_Match)) %>% 
  group_by(Sample,Base_1_Match) %>% 
  summarize(Sites = n()) %>% pivot_wider(names_from = 'Base_1_Match',values_from='Sites') %>%
  rowwise() %>%
  mutate(Proportion_CarIll = CarIll/(CarIll+Other+Other_Parent),
         Proportion_Other = Other/(CarIll+Other+Other_Parent)) %>%
  select(Sample,CarIll,Other_Parent,Proportion_CarIll,Other,Proportion_Other) %>% ungroup()

# SAVED AS 'Hybrid_Biallelic_Breakdown.tsv'
biallelic_df <- allele_cov_df %>% 
  filter(Variation == 'Biallelic') %>%
  mutate(Base_1_Match = ifelse(Base_1_Match %in% c('CarCor','CarAqu','CarLac','CarOva','CarMyr'),'Other_Parent',Base_1_Match),
         Base_2_Match = ifelse(Base_2_Match %in% c('CarCor','CarAqu','CarLac','CarOva','CarMyr'),'Other_Parent',Base_2_Match)) %>%
  mutate(Has_Other = ifelse(Base_1_Match == "Other" & Base_2_Match == "Other","Both",ifelse(Base_1_Match == "Other" | Base_2_Match == "Other","Mixed","Pure"))) %>%
  group_by(Sample,Has_Other) %>% summarize(Sites = n()) %>%
  rowwise() %>%
  pivot_wider(names_from='Has_Other',values_from='Sites') %>%
  mutate(Both = ifelse(is.na(Both),0,Both)) %>%
  mutate(Percent_Mixed = Mixed/(Mixed+Pure+Both),
         Percent_Pure = Pure/(Mixed+Pure+Both),
         Percent_Both = Both/(Mixed+Pure+Both))

#### Section 07: Figures ####

# Set colors for plots

species_list <- list('CarAqu'="CarAqu",
                     'CarCat'="CarCat",
                     'CarCor'="CarCor",
                     'CarIll'="CarIll",
                     'CarLac'="CarLac",
                     'CarMyr'="CarMyr",
                     'CarOva'="CarOva",
                     'CarPalm'="CarPalm")

species_list2 <- list('Carya aquatica'="Carya aquatica",
                      'Carya cathayensis'="Carya cathayensis",
                      'Carya cordiformis'="Carya cordiformis",
                      'Carya illioinensis'="Carya illioinensis",
                      'Carya laciniosa'="Carya laciniosa",
                      'Carya myrisiticiformis'="Carya myrisiticiformis",
                      'Carya ovata'="Carya ovata",
                      'Carya palmeri'="Carya palmeri")

species_color_list <- list('CarAqu'="#AA4499", # Purple
                           'CarCat'="#882255", # Maroon
                           'CarCor'="#332288", # Dark blue
                           'CarIll'="#DDCC77", # Yellow
                           'CarLac'="#88CCEE", # Light blue
                           'CarMyr'="#44AA99", # Light green
                           'CarOva'="#117733", # Dark green
                           'CarPalm'="#CC6677", # Salmon
                           'N.S.'="#DDDDDD")

species_color_list2 <- list('Carya aquatica'="#AA4499", # Purple
                            'Carya cathayensis'="#882255", # Maroon
                            'Carya cordiformis'="#332288", # Dark blue
                            'Carya illioinensis'="#DDCC77", # Yellow
                            'Carya laciniosa'="#88CCEE", # Light blue
                            'Carya myrisiticiformis'="#44AA99", # Light green
                            'Carya ovata'="#117733", # Dark green
                            'Carya palmeri'="#CC6677", # Salmon
                            'N.S.'="#DDDDDD")

species_colors <- c("#AA4499","#882255","#332288","#DDCC77","#88CCEE","#44AA99","#117733","#CC6677")

study_bar_plot_df <- study_species_stat_df %>% select(Study_Sample,Companion_Sample,Grouping_Sites,Significant_Bonferroni) %>%
  mutate(Plot_Color = species_color_list[Companion_Sample]) %>%
  separate(Study_Sample,into = c("Species","Sample"),remove = FALSE) %>%
  select(-Sample)

ms_species_plot_df <- study_bar_plot_df %>% filter(Study_Sample %in% c('CarAqu_1','CarCat_1','CarCor_2','CarIll_1','CarLac_1','CarMyr_1A','CarOva_1','CarPalm_2')) %>%
  mutate(Species = factor(Species,levels = c('CarLac','CarOva','CarMyr','CarPalm','CarCor','CarAqu','CarIll','CarCat')))

ms_species_bar_plot <- ggplot(ms_species_plot_df,aes(x=factor(Study_Sample,levels = rev(c('CarLac_1','CarOva_1','CarMyr_1A','CarPalm_2','CarCor_2','CarAqu_1','CarIll_1','CarCat_1'))),y=Grouping_Sites,fill=Plot_Color,pattern_fill = Plot_Color)) +
  geom_bar(stat="identity",position = "fill",width = 0.5,color="black")+
  scale_fill_identity() +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",plot.background=element_blank()) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),expand = expansion(mult=c(0,0.01)), limits = c(0, NA)) +
  xlab("Study Sample") +
  ylab("Percent of SNPs (%)") +
  coord_flip()


#### Figure S1: SNP Bar graphs for all samples ####

# NOTE: HAVE TO ADD LEGEND
fig_s1 <- ggplot(study_bar_plot_df,aes(x=Study_Sample,y=Grouping_Sites,fill=Plot_Color,pattern_fill = Plot_Color)) +
  geom_bar(stat="identity",position = "fill",width = 0.5)+
  scale_fill_identity() +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),expand = expansion(mult=c(0,0.01)), limits = c(0, NA)) +
  ylab("Proportion of SNPs Matching Reference Species") +
  facet_wrap(~Species,scales="free_x",nrow=2)

#### Figure S2: Downsampling grouping sites ####
fig_s2 <- ggplot(downsample_study_species_top_hit_df,aes(x=Bases,y=Grouping_Sites)) +
  geom_point(aes(color=Study_Sample)) +
  geom_smooth(method="lm",aes(fill=Study_Sample,color=Study_Sample)) +
  theme_bw() +
  facet_wrap(~Species,scales = "free_y") +
  xlab("Sample Bases") +
  ylab("SNPs Supporting Top Hit")

#### Figure S3: Downsampling grouping ratios ####
fig_s3 <- ggplot(downsample_study_species_top_hit_df,aes(x=Bases,y=Ratio)) +
  geom_point(aes(color=Study_Sample)) +
  geom_smooth(method="lm",aes(fill=Study_Sample,color=Study_Sample)) +
  theme_bw() +
  facet_wrap(~Species) +
  xlab("Sample Bases") +
  ylab("Proportion of SNPs Supporting Top Hit")

#### Figure 1B-D: Bubble plots for hybrids ####

hybrid_pooled_circle_plot_list <- list()

for(hybrid in hybrid_samples){
  
  pooled_plot_df <- hybrid_pooled_stat_df %>%
    filter(Grouping_Sites > 0,Hybrid_Sample == hybrid) %>% 
    rowwise() %>%
    mutate(Plot_Species = ifelse(Significant_Bonferroni == "Yes",Companion_Sample,"N.S."),
           Plot_Color = species_color_list[Plot_Species],
           Ratio_Label = ifelse(Significant_Bonferroni == "Yes",paste(c(Companion_Sample,"\n",round(Ratio,2),"%"),collapse = ""),"")) %>%
    select(Hybrid_Sample,Companion_Sample,Grouping_Sites,Ratio,Plot_Color,Ratio_Label)
  
  pooled_hybrid_grouping_sites <- sum(pooled_plot_df$Grouping_Sites)
  
  pooled_plot_color_list <- pull(pooled_plot_df,Plot_Color) %>% `names<-`(pooled_plot_df$Companion_Sample)
  
  pooled_packing <- circleProgressiveLayout(pooled_plot_df$Ratio, sizetype='area')
  
  pooled_plot_df <- cbind(pooled_plot_df, pooled_packing)
  
  pooled_ranked_matches <- pooled_plot_df$Companion_Sample
  names(pooled_ranked_matches)  <- c(1:length(pooled_ranked_matches))
  
  pooled_dat.gg <- circleLayoutVertices(pooled_packing,npoints = 100)
  pooled_dat.gg$Label <- pooled_ranked_matches[pooled_dat.gg$id]
  pooled_dat.gg <- pooled_dat.gg %>% rowwise() %>% mutate(Plot_Color = pooled_plot_color_list[Label]) %>% ungroup()
  
  hybrid_pooled_circle_plot_list[[hybrid]] <- ggplot() + 
    
    # Make the bubbles
    geom_polygon(data = pooled_dat.gg, colour = "black", alpha = 1,aes(x, y, group = as.factor(id),fill=Plot_Color)) +
    scale_color_identity() +
    # Add text in the center of each bubble + control its size
    geom_text(data = pooled_plot_df, aes(x, y, size=Ratio, label = Ratio_Label)) +
    scale_size_continuous(range = c(1,5)) +
    # General theme:
    theme_void() + 
    ggtitle(paste0(hybrid,"\nPooled: ",format(pooled_hybrid_grouping_sites,big.mark = ",")," sites")) +
    theme(legend.position="none",plot.title = element_text(hjust = 0.5)) +
    coord_equal()
}

# Figure S4: All hybrid bubble plots

# Pecan crosses
plot_myrxill <- hybrid_pooled_circle_plot_list[["myrxill_1"]]
plot_xbr <- tandemPlotter(hybrid_pooled_circle_plot_list[["xbr_1"]],hybrid_pooled_circle_plot_list[["xbr_2"]],hybrid_pooled_circle_plot_list[["xbr_3"]])
plot_xio <- hybrid_pooled_circle_plot_list[["xio_1"]]
plot_xnuss <- hybrid_pooled_circle_plot_list[["xnuss_1"]]
plot_xlc <- tandemPlotter(hybrid_pooled_circle_plot_list[["xlc_1"]],hybrid_pooled_circle_plot_list[["xlc_3"]],hybrid_pooled_circle_plot_list[["xlc_2"]],hybrid_pooled_circle_plot_list[["xlc_4"]])

all_pecan <- tandemPlotter(plot_myrxill,plot_xbr,plot_xio,plot_xnuss,plot_xlc,vertical = TRUE)
figure_pecan <- tandemPlotter(plot_myrxill,hybrid_pooled_circle_plot_list[["xbr_1"]],plot_xio,plot_xnuss,hybrid_pooled_circle_plot_list[["xlc_2"]])

# Other Crosses
plot_xbrl <- hybrid_pooled_circle_plot_list[["xbrl_1"]]
plot_xila <- tandemPlotter(hybrid_pooled_circle_plot_list[["xila_1"]],hybrid_pooled_circle_plot_list[["xila_2"]])

all_other <- tandemPlotter(plot_xbrl,plot_xila,vertical=TRUE)