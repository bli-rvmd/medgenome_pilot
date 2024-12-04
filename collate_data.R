library(dplyr)
library(openxlsx)
library(stringr)

####
# load data that was queried from cbioportal (also previously sent to MedGenome) for samples having both WES and RNA-Seq available
#### 
df_c <- read.csv("./MedGenome_metadata_sent/20241121_medgenome_samples_list.csv", header = T)


####
# load and check data sent from JJ
####
## 6236 dataset (N = 88) - 
df_6236 <- read.xlsx("./MedGenome_metadata_from_JJ/RMC-6236 MCT TV% change (cancer discovery).xlsx") %>%
  mutate(Sample_ID = str_remove_all(Model.Name, "-")) %>% # remove '-' in sample_id
  mutate(Sample_ID = str_to_upper(Sample_ID)) %>% # convert to all upper case
  mutate(Sample_ID = str_trim(Sample_ID)) %>% # trim off trailing whitespaces
  mutate(Sample_ID = gsub("LUC$", "", Sample_ID)) # rename LU99-luc to LU99

# JJ's 6236 samples not found in df_c
df_6236$Sample_ID[!df_6236$Sample_ID %in% df_c$Sample_ID]
# [1] "CTG0860"  "LXFA2889" "OVCAR5"  # likely due to not having WES or RNA-seq data or both?
  
  
## 6291 dataset (N = 40) - 
df_6291 <- read.xlsx("./MedGenome_metadata_from_JJ/RMC-6291.xlsx") %>%
  mutate(Sample_ID = str_remove_all(Models, "-")) %>% # remove '-' in sample_id
  mutate(Sample_ID = str_remove_all(Sample_ID, "\\*")) %>% # remove "*" in sample_id
  mutate(Sample_ID = str_to_upper(Sample_ID)) %>% # convert to all upper case
  mutate(Sample_ID = str_trim(Sample_ID)) # trim off trailing whitespaces
  
# JJ"s 6291 samples not found in df_c
df_6291$Sample_ID[!df_6291$Sample_ID %in% df_c$Sample_ID]
# [1] "ST1384" #likely to due typo of either ST1384? or ST1348? 


## 9805 dataset (N = 30)
df_9805 <- read.xlsx("./MedGenome_metadata_from_JJ/RMC-9805.xlsx") %>%
  mutate(Sample_ID = str_remove_all(Model, "-")) %>% # remove '-' in sample_id
  mutate(Sample_ID = str_to_upper(Sample_ID)) %>% # convert to all upper case
  mutate(Sample_ID = str_trim(Sample_ID)) # trim off trailing whitespaces
  
# JJ's 9805 samples not found in df_c
df_9805$Sample_ID[!df_9805$Sample_ID %in% df_c$Sample_ID]
# [1] "LXFA2889"


####
# backtrack samples contained in df_c but not found in JJ's data
####
## 6236
df_c_6236 <- df_c %>% 
  filter(!is.na(M_Recist_6236))

table(df_c_6236$Sample_ID %in% df_6236$Sample_ID)
# 
# FALSE  TRUE 
# 35    85 


## 6291
df_c_6291 <- df_c %>% 
  filter(!is.na(M_Recist_6291))

table(df_c_6291$Sample_ID %in% df_6291$Sample_ID)
# 
# FALSE  TRUE 
# 1    39 


## 9805
df_c_9805 <- df_c %>% 
  filter(!is.na(M_Recist_9805))

table(df_c_9805$Sample_ID %in% df_9805$Sample_ID)
# 
# TRUE 
# 28 



####
# collate JJ samples into df_c format 
# ST1384 -> ST1348
# LU99-luc -> LU99
# remove LXFA2889, CTG0860, and OVCAR5
####

df_6236_u <- df_6236 %>%
  filter(Sample_ID %in% df_c$Sample_ID) %>%
  arrange(Sample_ID) %>%
  select(Sample_ID, `RMC-6236.%.Mean.Tumor.Volume.change`)
colnames(df_6236_u) <- c("Sample_ID", "TV_Change_6236")
df_6236_u <- df_6236_u %>%
  mutate(TV_Change_6236 = round(TV_Change_6236))

# df_c_6236_u <- df_c_6236 %>%
#   filter(Sample_ID %in% df_6236_u$Sample_ID) %>%
#   arrange(Sample_ID)

df_6291_u <- df_6291 %>%
  mutate(Sample_ID = gsub("ST1384", "ST1348", Sample_ID)) %>%
  arrange(Sample_ID) %>%
  mutate(TV_Change_6291 = str_extract(`Tumor.volume,.%.change.(mean.±.sem)†`, "-?\\d+")) %>%
  select(Sample_ID, TV_Change_6291)


df_9805_u <- df_9805 %>%
  filter(Sample_ID %in% df_c$Sample_ID) %>%
  arrange(Sample_ID) %>%
  mutate(TV_Change_9805 = round(`RMC-9805.%.Mean.TV.change`)) %>%
  select(Sample_ID, TV_Change_9805)

## merge into df_c format
sample_vec <- unique(c(df_6236_u$Sample_ID, df_6291_u$Sample_ID, df_9805_u$Sample_ID))

df_c_u <- df_c %>%
  filter(Sample_ID %in% sample_vec) %>%
  arrange(Sample_ID) %>%
  select(Sample_ID, Cancer_Type, Ras_Genotype, Tumor_Volume_Change_6236, Tumor_Volume_Change_6291, Tumor_Volume_Change_9805)

## create output df
df_o <- df_c_u %>% # join TV change values from JJ"s data
  left_join(df_6236_u, by = "Sample_ID") %>%
  left_join(df_6291_u, by = "Sample_ID") %>%
  left_join(df_9805_u, by = "Sample_ID") %>%
  mutate( # use JJ's TV Change values
    Tumor_Volume_Change_6236 = TV_Change_6236, 
    Tumor_Volume_Change_6291 = TV_Change_6291, 
    Tumor_Volume_Change_9805 = TV_Change_9805, 
  ) %>%
  select(Sample_ID, Cancer_Type, Ras_Genotype, Tumor_Volume_Change_6236, Tumor_Volume_Change_6291, Tumor_Volume_Change_9805)

write.csv(df_o, "./MedGenome_metadata_sent/20241203_medgenome_samples_list.csv", sep = "\t", row.names = F)
