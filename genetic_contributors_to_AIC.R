

# 038 Papier CAI ----------------------------------------------------








# ___UKB Table variant LOF not HC cols -------------------------------------------

# annotate cols with LOF variants without HC, age, binaries, survival data

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)


# main df
df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260220_Summary_age_delays_diag_CH_PRS_REVEL_LOF_alphaMS.txt")

# Import new df with LOF count column
new_df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/001_RareVariants_allscores_bylist_20260414.txt") %>% 
  select(`Participant.ID`, WES_500k_LoF_MAF01__Hauck_ITPnHAI, WES_500k_LoF_MAF1__Hauck_ITPnHAI)

# Import PRS table to get indivs ages
prs_table <- read_xlsx("/home/stn/Documents/GitHub/BIN6007_Stage/project_icon_invitae/filtered_data/PRS/PRS_SLE_500k.xlsx") %>% 
  select(`Participant.ID`, `Age.at.recruitment` )



df_survival <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/004_20240421_SurvivalAnalysis_Censoring_final.txt", select = c("participant_id", "death_date", "death_censor_date", "hesin_last_date")) %>% 
  rename(`Participant.ID` = participant_id)


df_date_blood_drawing <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/4_500k_bloodcount.txt", select = c("Participant.ID", "Time.blood.sample.collected...Instance.0...Array.0"))



# Annotate LOF not CH cols
annotated_table <- left_join(df, new_df, by = c("Participant.ID")) %>% # Annotates LOF variants MAF < 0.1% and MAF < 1%
  
  # Replace NAs with 0 in LoF variants column
  mutate(WES_500k_LoF_MAF01__Hauck_ITPnHAI = ifelse(is.na(WES_500k_LoF_MAF01__Hauck_ITPnHAI), 0, WES_500k_LoF_MAF01__Hauck_ITPnHAI)) %>% 
  mutate(WES_500k_LoF_MAF1__Hauck_ITPnHAI = ifelse(is.na(WES_500k_LoF_MAF1__Hauck_ITPnHAI), 0, WES_500k_LoF_MAF1__Hauck_ITPnHAI)) %>% 
  
  # MAF < 0.1 %
  mutate("MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI" = WES_500k_LoF_MAF01__Hauck_ITPnHAI + `list_ITPnHAI_maf_MAF01_revel_0.9_count`) %>% 
  
  # MAF < 1 %
  mutate("MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI" = WES_500k_LoF_MAF1__Hauck_ITPnHAI + `list_ITPnHAI_maf_MAF1_revel_0.9_count` )


# Annotate age column
age_annotation <- left_join(annotated_table, prs_table, by = "Participant.ID") 



### Column annotations
# Annotate CH indivs
df_binary_CH_Age_sample <- age_annotation %>% 
  mutate(binary_CH_Age_sample = ifelse(!is.na(CH_Age_sample), 1, 0))

# Annotate ITP or AIHA binary
df_ITP_or_AIHA_binary <- df_binary_CH_Age_sample %>%
  mutate("binary_ITP_or_AIHA_indivs" = case_when(
    !is.na(TPI) ~ "1",
    !is.na(HAI) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 


# AIHA binary col
df_AIHA_binary <- df_ITP_or_AIHA_binary %>%
  mutate("binary_AIHA" = case_when(
    !is.na(HAI) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 



# Annotate binary ITP
df_ITP_binary <- df_AIHA_binary %>%
  mutate("binary_ITP" = case_when(
    !is.na(TPI) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 


# SLE binary 
df_binary_SLE <- df_ITP_binary %>% 
  mutate("binary_SLE" = case_when(
    !is.na(SLE) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 


# binary LYMPHOME
df_binary_lymphome <- df_binary_SLE %>%
  mutate("binary_LYMPHOME" = case_when(
    !is.na(LYMPHOME) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 

# binary Other_LEUK
df_binary_Other_LEUK <- df_binary_lymphome %>%
  mutate("binary_Other_LEUK" = case_when(
    !is.na(Other_LEUK) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 

# binary LLC
df_binary_LLC <- df_binary_Other_LEUK %>%
  mutate("binary_LLC" = case_when(
    !is.na(LLC) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 



# Create col with lower age value for cells containing multiple ages
df_lower_age_AIHA <- df_binary_LLC %>%
  mutate(AIHA_lower_age = sapply(strsplit(as.character(HAI), ";"), function(x) min(as.numeric(x))))

df_lower_age_ITP <- df_lower_age_AIHA %>%
  mutate(TPI_lower_age = sapply(strsplit(as.character(TPI), ";"), function(x) min(as.numeric(x))))

df_lower_age_CAI <- df_lower_age_ITP %>% 
  mutate(CAI_lower_age = pmin(TPI,HAI, na.rm = T))

df_lower_age_SLE <- df_lower_age_CAI %>%
  mutate(SLE_lower_age = sapply(strsplit(as.character(SLE), ";"), function(x) min(as.numeric(x))))

df_lower_age_lymphoma <- df_lower_age_SLE %>%
  mutate(LYMPHOME_lower_age = sapply(strsplit(as.character(LYMPHOME), ";"), function(x) min(as.numeric(x))))

df_lower_age_lymphoma_Other_LEUK <- df_lower_age_lymphoma %>% 
  mutate(Other_LEUK_lower_age = sapply(strsplit(as.character(Other_LEUK), ";"), function(x) min(as.numeric(x))))

df_lower_age_lymphoma_Other_LEUK_LLC <- df_lower_age_lymphoma_Other_LEUK %>% 
  mutate(LLC_lower_age = sapply(strsplit(as.character(LLC), ";"), function(x) min(as.numeric(x))))

df_lower_age_lymphoma_Other_LEUK_LLC_HC <- df_lower_age_lymphoma_Other_LEUK_LLC %>% 
  mutate(CH_Age_sample_lower_age = sapply(strsplit(as.character(CH_Age_sample), ";"), function(x) min(as.numeric(x))))




### Annotate survival analysis data 
# Annotate time of blood sample collection (age_inclusion)
last_annotated_df <- df_lower_age_lymphoma_Other_LEUK_LLC_HC 

df_age_inclusion <- left_join(last_annotated_df, df_date_blood_drawing, by = "Participant.ID") %>% # Annotate date of blood drawing
  mutate(blood_draw_date = sub("T.*$", "", Time.blood.sample.collected...Instance.0...Array.0)) %>% # Formats date of blood drawing by removing the hour of drawing
  mutate(formatted_DATE_BIRTH = parse_date_time(DATE_BIRTH, orders = c("dmy", "ymd"))) %>% # Formats the date of birth to YYYY-MM-DD
  mutate(formatted_DATE_BIRTH = format(formatted_DATE_BIRTH, "%Y-%m-%d")) %>% # Formats the date of birth to YYYY-MM-DD
  mutate(age_inclusion = as.numeric(as.Date(blood_draw_date) - as.Date(formatted_DATE_BIRTH))) # Calculates the age at the inclusion

# Annotate date and age end of follow-up 
df_end_follow_up <- left_join(df_age_inclusion, df_survival, by = "Participant.ID") %>% 
  mutate(date_end_follow_up = (pmax(death_date, death_censor_date, hesin_last_date, na.rm = T))) %>% 
  select(-death_date, -death_censor_date, -hesin_last_date) %>% 
  mutate(age_end_follow_up = as.numeric(as.Date(date_end_follow_up) - as.Date(formatted_DATE_BIRTH)))

# Indivs to be used as control
# Annotate span age end of follow-up - age inclusion
df_span_age_end_follow_up_x_age_inclusion <- df_end_follow_up %>% 
  mutate(span_age_end_follow_up_age_inclusion = age_end_follow_up - age_inclusion)






# Write output

output <- df_span_age_end_follow_up_x_age_inclusion

print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# ___AoU master table --------------------------------------------------------


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




### Import files
df_main <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/data/cohort/DateEndFollow_reliable.tsv") %>% 
  rename(date_end_follow_up = followup_end_conservative,
         age_end_follow_up = age_at_followup)

df_sequenced_indivs <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/results/1_diagnostics/autoimmune_diseases_CDRv8/wgs_samples.tsv")

df_ITP_693_indivs <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/results/1_diagnostics/autoimmune_diseases_CDRv8/TPI_693_only.tsv")

df_AIHA_indivs <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/results/1_diagnostics/autoimmune_diseases_CDRv8/AHAI_only.tsv")

df_general_info <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/data/cohort/general_info_allparticipants_CDRv8.tsv", 
                         select = c("person_id","sex_genetic", "year_of_birth"))

df_sle_indivs <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/results/1_diagnostics/autoimmune_diseases_CDRv8/summary_by_diag/SLE_participant.txt")

df_llc_indivs <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/results/1_diagnostics/autoimmune_diseases_CDRv8/summary_by_diag/LLC_participant.txt")

df_lymphome_indivs <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/results/1_diagnostics/autoimmune_diseases_CDRv8/summary_by_diag/LYMPHOME_participant.txt")

df_other_leuk_indivs <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/results/1_diagnostics/autoimmune_diseases_CDRv8/summary_by_diag/Other_LEUK_participant.txt")

df_LoF_MAF01__Hauck_ITPnHAI <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/002_AoU_LoF_MAF01__Hauck_ITPnHAI.tsv", 
                                     select = c("s", "AoU_LoF_MAF01__Hauck_ITPnHAI")) %>% 
  distinct(s, .keep_all = T)

df_LoF_MAF1__Hauck_ITPnHAI <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/005_AoU_LoF_MAF1__Hauck_ITPnHAI.tsv", 
                                    select = c("s", "AoU_LoF_MAF1__Hauck_ITPnHAI")) %>% 
  distinct(s, .keep_all = T)

df_CHIP <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/003_AoU_CHIP.tsv")

df_PGS_196 <- fread("/home/rstudio/workspace/workspace-bucket/data/from_workspace_PRS_calculation/PGS000196_scores.csv") %>% 
  select(-N_variants) %>% 
  rename(PGS000196 = sum_weights)

df_PGS_4917 <- fread("/home/rstudio/workspace/workspace-bucket/data/from_workspace_PRS_calculation/PGS004917_scores.csv") %>% 
  select(-N_variants) %>% 
  rename(PGS004917 = sum_weights)

df_all_covariables <- fread("/home/rstudio/workspace/workspace-bucket/data/cohort/20260618_all_covariables.tsv") %>% 
  select(person_id, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, biosample_collection_age)

df_first_condition <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/data/cohort/archive/lastcondition_allparticipants.tsv") %>% 
  select(person_id, first_condition_date)






### Code




# Keep indivs with sequencing only
df_main_sequenced <- inner_join(df_main, df_sequenced_indivs, by = c("person_id" = "V1")) # 270634 indivs





### Annotate cols


# General data
# Sex, year of birth
df_sex_indivs <- left_join(df_main_sequenced, df_general_info, by = c("person_id"))

# First condition date
df_first_condition_date <- left_join(df_sex_indivs, 
                                     df_first_condition, 
                                     by = c("person_id")) # 10687 indivs doesn't have first condition annotation

# First condition year
df_first_condition_year <- df_first_condition_date %>%
  mutate(first_condition_year = year(first_condition_date))

# First condition age
df_first_condition_age <- df_first_condition_year %>% 
  mutate(first_condition_age = (first_condition_year - year_of_birth) * 365.25)

# age_end_follow_up in days
df_age_end_follow_up <- df_first_condition_age %>% 
  mutate(age_end_follow_up = age_end_follow_up*365.25)




# Covariables
# LoF variants
df_lof_per_indiv <- left_join(df_age_end_follow_up, df_LoF_MAF01__Hauck_ITPnHAI, by = c("person_id" = "s")) %>% 
  mutate(AoU_LoF_MAF01__Hauck_ITPnHAI = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0)) # Replace NAs with zero

df_lof_1percent_per_indiv <- left_join(df_lof_per_indiv, df_LoF_MAF1__Hauck_ITPnHAI, by = c("person_id" = "s")) %>% 
  mutate(AoU_LoF_MAF1__Hauck_ITPnHAI = replace_na(AoU_LoF_MAF1__Hauck_ITPnHAI, 0)) # Replace NAs with zero


# SLE PGS 196
df_main_PGS_196 <- left_join(df_lof_1percent_per_indiv, df_PGS_196, by = c("person_id" = "s"))

# SLE PGS 4917
df_main_PGS_4917 <- left_join(df_main_PGS_196, df_PGS_4917, by = c("person_id" = "s"))

# 10 PCs
df_PCs <- left_join(df_main_PGS_4917, df_all_covariables, by = c("person_id"))





# Diagnostics
# Annotate ITP
df_main_itp <- left_join(df_PCs, df_ITP_693_indivs, by = c("person_id")) %>% 
  rename(TPI = TPI_693)

# Annotate AIHA
df_main_aiha <- left_join(df_main_itp, df_AIHA_indivs, by = c("person_id")) %>% 
  rename(HAI = AHAI)

# Annotate SLE
df_main_sle <- left_join(df_main_aiha, df_sle_indivs, by = c("person_id"))

# Annotate LLC
df_main_llc <- left_join(df_main_sle, df_llc_indivs, by = c("person_id"))

# Annotate lymphome
df_main_lymphome <- left_join(df_main_llc, df_lymphome_indivs, by = c("person_id"))

# Annotate other_leuk
df_main_other_leuk <- left_join(df_main_lymphome, df_other_leuk_indivs, by = c("person_id"))

# Annotate CHIP
df_ch_age_sample <- left_join(df_main_other_leuk, df_CHIP, by = c("person_id"))

# Span age of inclusion - age end of follow-up
df_span_age_inclusion_end_follow_up <- df_ch_age_sample %>% 
  mutate(span_age_end_follow_up_age_inclusion = age_end_follow_up - biosample_collection_age/365.25)

# Age at recruitment
df_age_recruitment <- df_span_age_inclusion_end_follow_up %>% 
  mutate(`Age.at.recruitment` = biosample_collection_age/365.25)





# Binary cols
# ITP or AIHA 
df_ITP_or_AIHA_binary <- df_age_recruitment %>%
  mutate("binary_ITP_or_AIHA_indivs" = case_when(
    !is.na(TPI) ~ "1",
    !is.na(HAI) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 


# AIHA 
df_AIHA_binary <- df_ITP_or_AIHA_binary %>%
  mutate("binary_AIHA" = case_when(
    !is.na(HAI) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 


# ITP
df_ITP_binary <- df_AIHA_binary %>%
  mutate("binary_ITP" = case_when(
    !is.na(TPI) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 


# SLE
df_sle_binary <- df_ITP_binary %>% 
  mutate("binary_SLE" = case_when(
    !is.na(SLE) ~ "1",
    TRUE ~ "0"
  ))


# LLC
df_llc_binary <- df_sle_binary %>% 
  mutate("binary_LLC" = case_when(
    !is.na(LLC) ~ "1",
    TRUE ~ "0"
  ))


# lymphoma
df_lymphoma_binary <- df_llc_binary %>% 
  mutate("binary_LYMPHOME" = case_when(
    !is.na(LYMPHOME) ~ "1",
    TRUE ~ "0"
  ))


# other leuk
df_other_leuk_binary <- df_lymphoma_binary %>% 
  mutate("binary_Other_LEUK" = case_when(
    !is.na(Other_LEUK) ~ "1",
    TRUE ~ "0"  # keep original value for other cases
  )) 




# Annotate CH indivs
df_binary_CH_Age_sample <- df_other_leuk_binary %>% 
  mutate(binary_CH_Age_sample = ifelse(!is.na(CH_Age_sample), 1, 0))





# Lower age

# Create col with lower age value for cells containing multiple ages
df_lower_age_AIHA <- df_binary_CH_Age_sample %>%
  mutate(AIHA_lower_age = sapply(strsplit(as.character(HAI), ";"), function(x) min(as.numeric(x))))

df_lower_age_ITP <- df_lower_age_AIHA %>%
  mutate(TPI_lower_age = sapply(strsplit(as.character(TPI), ";"), function(x) min(as.numeric(x))))

df_lower_age_CAI <- df_lower_age_ITP %>% 
  mutate(CAI_lower_age = pmin(TPI_lower_age,AIHA_lower_age, na.rm = T))

df_lower_age_SLE <- df_lower_age_CAI %>%
  mutate(SLE_lower_age = sapply(strsplit(as.character(SLE), ";"), function(x) min(as.numeric(x))))

df_lower_age_lymphoma <- df_lower_age_SLE %>%
  mutate(LYMPHOME_lower_age = sapply(strsplit(as.character(LYMPHOME), ";"), function(x) min(as.numeric(x))))

df_lower_age_lymphoma_Other_LEUK <- df_lower_age_lymphoma %>% 
  mutate(Other_LEUK_lower_age = sapply(strsplit(as.character(Other_LEUK), ";"), function(x) min(as.numeric(x))))

df_lower_age_lymphoma_Other_LEUK_LLC <- df_lower_age_lymphoma_Other_LEUK %>% 
  mutate(LLC_lower_age = sapply(strsplit(as.character(LLC), ";"), function(x) min(as.numeric(x))))

df_lower_age_lymphoma_Other_LEUK_LLC_HC <- df_lower_age_lymphoma_Other_LEUK_LLC %>% 
  mutate(CH_Age_sample_lower_age = sapply(strsplit(as.character(CH_Age_sample), ";"), function(x) min(as.numeric(x))))





# write file
output <- df_lower_age_lymphoma_Other_LEUK_LLC_HC 

names(output)

# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')













# ____AoU_LoF_MAF01__Hauck_ITPnHAI --------------------------------------------

######### Preprocessing raw data
### LOF annotation


### Correct the premature /n occurring in some lines in all files lof_carriers_chr*.tsv
for file in lof_carriers_chr*.tsv; do

sed -i ':a;N;$!ba;s/\n / /g' "$file"

done





### Merge parts
df_final <- data.table()

for (chr in 1:22){
  
  df_temp <- paste0("/home/rstudio/workspace/workspace-bucket/stennio/lof_genotypes/lof_carriers_chr", chr, ".tsv") %>% fread()
  df_temp$chr = chr
  df_final <- rbind(df_final, df_temp)
  
}

chrx <- fread("/home/rstudio/workspace/workspace-bucket/stennio/lof_genotypes/lof_carriers_chrX.tsv")
chrx$chr = "x"
df_final <- rbind(df_final, chrx)




# write file
output <- df_final 

# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/001_all_chr_lof_carriers.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')

######### Preprocessing raw data






### Filter genes in ITP or AIHA
df_main <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/001_all_chr_lof_carriers.tsv") %>% 
  mutate(Gene = gene_symbol)

hauck_list <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/IEI_Hauck_2024_FASLG.csv")




# Identifying genes associated with (ITP OR AIHA)  in the dataset

# Create new columns for further filtering
df_genes_associatied_col <- df_main %>%
  mutate(itp_associated = "no",
         aiha_associated = "no" ) 

hauck_ITP <- hauck_list %>%
  filter(ITP != "") %>%
  select(ITP) %>%
  distinct(.) %>%
  pull()
hauck_ITP <- unlist(strsplit(hauck_ITP, "[;,]")) # Splits FASL;FASLG cell


hauck_AIHA <- hauck_list %>%
  filter(AIHA != "") %>%
  select(AIHA) %>%
  distinct(.) %>%
  pull()
hauck_AIHA <- unlist(strsplit(hauck_AIHA, "[;,]")) # Splits FASL;FASLG cell



# Annotate genes associated with ITP
df_genes_in_ITP <- df_genes_associatied_col %>%
  mutate(
    itp_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_ITP), "yes", "no")
    })
  )


# Annotate genes associated with AIHA
df_genes_in_AIHA <- df_genes_in_ITP %>%
  mutate(
    aiha_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_AIHA), "yes", "no")
    })
  )


# Keep only Hauck ITP or AIHA genes
df_genes_in_CAI <- df_genes_in_AIHA %>%
  filter(itp_associated == "yes" | aiha_associated == "yes")

## Patient with 13 lof variants
# patient_13_variants <- df_genes_in_CAI %>% 
#   filter(s == "5685880")
# 
# 
# # write file
# output <- patient_13_variants 
# 
# # Write output
# print('writing final table')
# write.table(output,
#             file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/004_indiv_13_varfiants.tsv" ,
#             sep = "\t", row.names = F, col.names = T , quote = F)
# 
# print('final table written')



### Create keys to merge with AF table

# Create col pos
df_pos <- df_genes_in_CAI %>% 
  separate(locus, 
           into = c("chrom", "pos"), 
           sep = ":", remove = F)

# Create cols ref and alt
df_ref_alt <- df_pos %>% 
  mutate(alleles_treated = str_remove_all(alleles, '\\[|\\]|"')) %>%  # Removes [ ] and "
  separate(alleles_treated, 
           into = c("ref", "alt"), 
           sep = ",",
           remove = F)

# Correct chrx col
df_chrX <- df_ref_alt %>% 
  mutate(chr = str_replace(chr, "x", "X"))


# Create hg38ID keys
df_keys <- df_chrX %>% 
  unite(hg38ID, chr, pos, ref, alt, sep = "-") %>% 
  select(-chrom, -alleles_treated, -itp_associated, -aiha_associated, -Gene)




###  Prepare AF table for merge

# AF
df_lof_AF <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/results/2_variants/20260506_lof_variants_517genes_autoimmune.tsv")
# The df has duplicated variant lines because the multiple vals in the gene col was splitted to unique values in the col gene_symbol



# Filter MAF < 0.1%
df_maf_0p1 <- df_lof_AF %>% 
  filter(allele_frequency <= 0.001)


# Select cols
df_maf_cols <- df_maf_0p1 %>% 
  select(vid, allele_frequency, genes) %>% 
  distinct(vid, .keep_all = T)





### Fusion tables
df_lof_fusion <- left_join(df_keys, df_maf_cols, by = c("hg38ID" = "vid"))

# nas <- df_lof_fusion %>% 
#   filter(is.na(allele_frequency))
# 2277 variants have no allele_frequency


multiple_indivs <- df_lof_fusion %>% 
  group_by(s) %>% 
  filter(n() > 1) %>% 
  ungroup()






### Create table LoF variants CAI MAF < 0.1%
df_CAI_maf_0p1 <- df_lof_fusion %>% 
  filter(!is.na(allele_frequency))





### Annotate number of copies for each variant
df_number_copies_each_variant <- df_CAI_maf_0p1 %>%
  mutate(variant_quantity = case_when(
    GT == "0/1" ~ "1",
    GT == "0|1" ~ "1",
    GT == "1|0" ~ "1",
    GT == "1" ~ "2",
    GT == "1/1" ~ "2"
  ))


### Annotate the number of LoF variants per indiv
df_lof_per_indiv <- df_number_copies_each_variant %>% 
  group_by(s) %>% 
  mutate(AoU_LoF_MAF01__Hauck_ITPnHAI = sum(as.numeric(variant_quantity))) %>% 
  ungroup()




# write file
output <- df_lof_per_indiv 

# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/002_AoU_LoF_MAF01__Hauck_ITPnHAI.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# ____AoU_LoF_MAF1__Hauck_ITPnHAI --------------------------------------------



### Filter genes in ITP or AIHA
df_main <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/001_all_chr_lof_carriers.tsv") %>% 
  mutate(Gene = gene_symbol)

hauck_list <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/IEI_Hauck_2024_FASLG.csv")




# Identifying genes associated with (ITP OR AIHA)  in the dataset

# Create new columns for further filtering
df_genes_associatied_col <- df_main %>%
  mutate(itp_associated = "no",
         aiha_associated = "no" ) 

hauck_ITP <- hauck_list %>%
  filter(ITP != "") %>%
  select(ITP) %>%
  distinct(.) %>%
  pull()
hauck_ITP <- unlist(strsplit(hauck_ITP, "[;,]")) # Splits FASL;FASLG cell


hauck_AIHA <- hauck_list %>%
  filter(AIHA != "") %>%
  select(AIHA) %>%
  distinct(.) %>%
  pull()
hauck_AIHA <- unlist(strsplit(hauck_AIHA, "[;,]")) # Splits FASL;FASLG cell



# Annotate genes associated with ITP
df_genes_in_ITP <- df_genes_associatied_col %>%
  mutate(
    itp_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_ITP), "yes", "no")
    })
  )


# Annotate genes associated with AIHA
df_genes_in_AIHA <- df_genes_in_ITP %>%
  mutate(
    aiha_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_AIHA), "yes", "no")
    })
  )


# Keep only Hauck ITP or AIHA genes
df_genes_in_CAI <- df_genes_in_AIHA %>%
  filter(itp_associated == "yes" | aiha_associated == "yes")

## Patient with 13 lof variants
# patient_13_variants <- df_genes_in_CAI %>% 
#   filter(s == "5685880")
# 
# 
# # write file
# output <- patient_13_variants 
# 
# # Write output
# print('writing final table')
# write.table(output,
#             file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/004_indiv_13_varfiants.tsv" ,
#             sep = "\t", row.names = F, col.names = T , quote = F)
# 
# print('final table written')



### Create keys to merge with AF table

# Create col pos
df_pos <- df_genes_in_CAI %>% 
  separate(locus, 
           into = c("chrom", "pos"), 
           sep = ":", remove = F)

# Create cols ref and alt
df_ref_alt <- df_pos %>% 
  mutate(alleles_treated = str_remove_all(alleles, '\\[|\\]|"')) %>%  # Removes [ ] and "
  separate(alleles_treated, 
           into = c("ref", "alt"), 
           sep = ",",
           remove = F)

# Correct chrx col
df_chrX <- df_ref_alt %>% 
  mutate(chr = str_replace(chr, "x", "X"))


# Create hg38ID keys
df_keys <- df_chrX %>% 
  unite(hg38ID, chr, pos, ref, alt, sep = "-") %>% 
  select(-chrom, -alleles_treated, -itp_associated, -aiha_associated, -Gene)




###  Prepare AF table for merge

# AF
df_lof_AF <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/results/2_variants/20260506_lof_variants_517genes_autoimmune.tsv")
# The df has duplicated variant lines because the multiple vals in the gene col was splitted to unique values in the col gene_symbol




# Filter MAF < 1%
df_maf_1 <- df_lof_AF %>% 
  filter(allele_frequency <= 0.01)


# Select cols
df_maf_cols <- df_maf_1 %>% 
  select(vid, allele_frequency, genes) %>% 
  distinct(vid, .keep_all = T)





### Fusion tables
df_lof_fusion <- left_join(df_keys, df_maf_cols, by = c("hg38ID" = "vid"))

# nas <- df_lof_fusion %>% 
#   filter(is.na(allele_frequency))
# 2277 variants have no allele_frequency






### Create table LoF variants CAI MAF < 1%
df_CAI_maf_1 <- df_lof_fusion %>% 
  filter(!is.na(allele_frequency))





### Annotate number of copies for each variant
df_number_copies_each_variant <- df_CAI_maf_1 %>%
  mutate(variant_quantity = case_when(
    GT == "0/1" ~ "1",
    GT == "0|1" ~ "1",
    GT == "1|0" ~ "1",
    GT == "1" ~ "2",
    GT == "1/1" ~ "2"
  ))


### Annotate the number of LoF variants per indiv
df_lof_per_indiv <- df_number_copies_each_variant %>% 
  group_by(s) %>% 
  mutate(AoU_LoF_MAF1__Hauck_ITPnHAI = sum(as.numeric(variant_quantity))) %>% 
  ungroup()




# write file
output <- df_lof_per_indiv 

# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/005_AoU_LoF_MAF1__Hauck_ITPnHAI.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# ____HC variants -------------------------------------------------------------

#all_hc_variants <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/data/CH/all_CHIP_calls_02082023.txt")

df_chip_per_person <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/data/CH/CHIP_calls_per_person_02082023.txt")


df_all_covariables <- fread("/home/rstudio/workspace/workspace-bucket/data/cohort/20260618_all_covariables.tsv") %>% 
  select(person_id, biosample_collection_age) %>% 
  distinct(person_id, .keep_all = T)


# df_first_condition <- fread("/home/rstudio/workspace/workspace-bucket/Project_Rare_Variants/data/cohort/archive/lastcondition_allparticipants.tsv") %>% 
#   select(person_id, first_condition_date)
# We can't use this because we don't know if the indivs had the clone at the time of first condition


# Filter HC indivs
df_hc_indivs <- df_chip_per_person %>% 
  filter(CHIP == 1) %>% 
  select(person_id, CHIP)



# Annotate age at blood draw
df_age_blood_draw <- left_join(df_hc_indivs, df_all_covariables, by = c("person_id")) %>% 
  rename(CH_Age_sample = biosample_collection_age) %>% 
  select(-CHIP)
# 4630 indivs





nas <- df_age_blood_draw %>% 
  filter(is.na(CH_Age_sample))


# write file
output <- df_age_blood_draw 

# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/003_AoU_CHIP.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')







# ___001 VARIANTS RARE ---------------------------------------------------


# ____001 CAI - variants dans genes associees aux CAI (liste Hauck ITP + AIHA ensemble mais sans SLE) --------
# maf <0.1 et maf <1%
# LOF seul, Revel 0.9 et 0.5 seul, LOF + Revel 0.9
# Faire en régression logistique (ajusté sur sexe, âge à la dernière visite/décès et 10 PC (pas besoin de faire de chi2/Fisher)
# Faire également le test en ne considerant pas les variants TET2 (donc la meme liste mais sans variant TET2)


# _____001a MAF < 0.1% -------------------------------------------------------------------
# Figure 1A

### MAF < 0.1, LOF seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)






model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)






### MAF < 0.1, REVEL 0.5 seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF01_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF01_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)









### MAF < 0.1, REVEL 0.9 seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df
nrow(df_subset)




# Print to console
summary(glm(
  
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF01_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)



# Calculate OR, CI
model_logistic_regression <- glm(
  
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF01_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)












### MAF < 0.1, LOF + Revel 0.9, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df
nrow(df_subset)




# Print to console
summary(glm(
  
  as.factor(binary_ITP_or_AIHA_indivs) ~  MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)



# Calculate OR, CI
model_logistic_regression <- glm(
  
  as.factor(binary_ITP_or_AIHA_indivs) ~  MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)









# _________AoU Cox-------------------------------------------------------------------

# follow-up = age_inclusion until CAI for cases


### MAF < 0.1, LOF seul, TET2 included



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_stratified_score <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )




# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)










# _____001b MAF < 1% -------------------------------------------------------------------


### MAF 1 %, LOF seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF1__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF1__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)









### MAF < 1%, REVEL 0.5 seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF1_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF1_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)









### MAF < 1% , REVEL 0.9 seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df
nrow(df_subset)




# Print to console
summary(glm(
  
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF1_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)



# Calculate OR, CI
model_logistic_regression <- glm(
  
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF1_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)












### MAF < 1% , LOF + Revel 0.9, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df
nrow(df_subset)




# Print to console
summary(glm(
  
  as.factor(binary_ITP_or_AIHA_indivs) ~  MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)



# Calculate OR, CI
model_logistic_regression <- glm(
  
  as.factor(binary_ITP_or_AIHA_indivs) ~  MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)











# ____002 ITP variants dans genes associees a TPI ------------------------------



# _____001a MAF < 0.1% ----------------------------------------------------
# Figure 1A


### MAF < 0.1, LOF seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







### MAF < 0.1, REVEL 0.5 seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF01_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF01_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







### MAF < 0.1, REVEL 0.9 seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF01_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF01_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)













### MAF < 0.1, REVEL 0.9 or LoF variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










# _________AoU Cox-------------------------------------------------------------------

# follow-up = age_inclusion until CAI for cases


### MAF < 0.1, LOF seul, TET2 included



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP == 1, 
                                     TPI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_stratified_score <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    lower_age_years = TPI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP == 1,
                                      lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )




# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)















# _____001b MAF < 1% ------------------------------------------------------




### MAF < 1, LoF only variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF1__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF1__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)











### MAF < 1, REVEL > 0.5 variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF1_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF1_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)









### MAF < 1, REVEL > 0.9 variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF1_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF1_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








### MAF < 1, REVEL > 0.9 or LoF variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)





















# ____003 AIHA - variants dans genes associees a HAI -----------------------------
# Figure 1A

# _____001a MAF < 0.1% ----------------------------------------------------


### MAF < 0.1, LOF seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










### MAF < 0.1, REVEL 0.5 seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF01_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF01_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)











### MAF < 0.1, REVEL 0.5 seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF01_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF01_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








### MAF < 0.1, REVEL 0.9 seul, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF01_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF01_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)









### MAF < 0.1, REVEL 0.9 or LoF variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)












# _________AoU Cox-------------------------------------------------------------------

# follow-up = age_inclusion until CAI for cases


### MAF < 0.1, LOF seul, TET2 included



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_AIHA == 1, 
                                     AIHA_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_stratified_score <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    lower_age_years = AIHA_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_AIHA == 1,
                                      lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )




# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_AIHA) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)












# _____001b MAF < 1% ------------------------------------------------------







### MAF < 1%,  LoF only  variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF1__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF1__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










### MAF < 1%,  REVEL > 0.5 only  variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF1_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF1_revel_0.5_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)












### MAF < 1%,  REVEL > 0.9 only  variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF1_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF1_revel_0.9_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










### MAF < 1%,  REVEL > 0.9  or LoF variants, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



















# ____003 âge (moyenne, SD et t-test) au diagnostic de la première CAI pour les patients avec LOF et sans LOF  -----------------------------------------------------------------
# 
# A faire @Stennio Faria (liste Hauck ITP+AIHA maf 0.1%): 
# âge (moyenne, SD et t-test) au diagnostic de la première CAI pour les patients avec LOF et sans LOF 
# 
# Nombre de patients avec 1, 2 et 3 variants


# _____003a CAI and LoF mean diag age comparsion --------------------------------------------------


# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)




df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# Subset CAI
df_cai_indivs <- df %>% 
  filter(binary_ITP_or_AIHA_indivs == 1)



### Subset CAI AND LOF
df_cai_and_lof <- df_cai_indivs %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0)



# Mean age CAI AND LOF
mean_age_df_cai_and_lof <- (mean(df_cai_and_lof$CAI_lower_age)) /365.25
mean_age_df_cai_and_lof
# 66.50212 years


# SD CAI AND LOF
sd_age_df_cai_and_lof <- (sd(df_cai_and_lof$CAI_lower_age)) /365.25
sd_age_df_cai_and_lof
# 10.66283 years




### Subset CAI NOT LOF
df_cai_not_lof <- df_cai_indivs %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI == 0)

# Mean age CAI NOT LOF
mean_age_df_cai_not_lof <- (mean(df_cai_not_lof$CAI_lower_age)) /365.25
mean_age_df_cai_not_lof

# SD CAI NOT LOF
sd_age_df_cai_not_lof <- (sd(df_cai_not_lof$CAI_lower_age)) /365.25
sd_age_df_cai_not_lof




### t-test age comparison
t_test <- t.test(x = df_cai_and_lof$CAI_lower_age, y = df_cai_not_lof$CAI_lower_age)
t_test




### Number of variants per indiv in CAI patients
cai_variants_per_patient <- table(df_cai_indivs$WES_500k_LoF_MAF01__Hauck_ITPnHAI) %>% as.data.frame() %>% 
  rename("No_variants" = Var1,
         "No_patients" = Freq)
cai_variants_per_patient



# ______ Special task: 37 indivs x 33 indivs in my data -----------------------------------------------------
# 
# Salut!
#   Quand je regarde la liste des analyses par gène d’Estelle j’ai 37 variants mais toi tu n’as que 33 individus avec CAI qui ont un variant. Est-ce que c’est parce que certains ont eu ITP et AIHA?
# 
#   Non tout est bon:
#   33 patients avec variants dont 1 avec 2 variants et 2 avec PTI + AHAI
# 
# Ca fait 36 variants si tu analyse séparement PTI et AHAI
# 
# Tu peux juste me donner les variants de celui qui en a 2 et des deux qui ont PTI + AHAI stp?
#   





# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)




df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


df_bygenes <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/002_RareVariants_allscores_bygenes_20260414_no_HC.txt", select = c("Participant.ID", "Gene", "WES_500k_LoF_MAF01"))


# Subset CAI
df_cai_indivs <- df %>% 
  filter(binary_ITP_or_AIHA_indivs == 1)


indivs_with_tpi_and_aiha <- df %>% 
  filter(!is.na(TPI) & !is.na(HAI)) %>% 
  select(`Participant.ID`, TPI, HAI, WES_500k_LoF_MAF01__Hauck_ITPnHAI)
# 24 indivs

indivs_with_variants <- indivs_with_tpi_and_aiha %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0)

indivs_with_variants
# Participant.ID
# 1279173
# 4107363


indiv_with_2_vars <- df_cai_indivs %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI >1)
indiv_with_2_vars
# Participant.ID
# 1671995

# variants in the individuals:

indivs_genes <- df_bygenes %>% 
  filter(`Participant.ID` == 1279173 | 
           `Participant.ID` == 4107363 | 
           `Participant.ID` == 1671995) %>% 
  filter(WES_500k_LoF_MAF01 > 0)

indivs_genes






# ITP
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)




df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )



df_ITP_indivs <- df %>% 
  filter(binary_ITP == 1)


### Subset ITP AND LOF
df_ITP_and_lof <- df_ITP_indivs %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0)



# Mean age ITP AND LOF
mean_age_df_ITP_and_lof <- (mean(df_ITP_and_lof$TPI_lower_age)) /365.25
mean_age_df_ITP_and_lof



# SD ITP AND LOF
sd_age_df_ITP_and_lof <- (sd(df_ITP_and_lof$TPI_lower_age)) /365.25
sd_age_df_ITP_and_lof





### Subset ITP NOT LOF
df_ITP_not_lof <- df_ITP_indivs %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI == 0)

# Mean age ITP NOT LOF
mean_age_df_ITP_not_lof <- (mean(df_ITP_not_lof$TPI_lower_age)) /365.25
mean_age_df_ITP_not_lof

# SD ITP NOT LOF
sd_age_df_ITP_not_lof <- (sd(df_ITP_not_lof$TPI_lower_age)) /365.25
sd_age_df_ITP_not_lof




### t-test age comparison
t_test <- t.test(x = df_ITP_and_lof$TPI_lower_age, y = df_ITP_not_lof$TPI_lower_age)
t_test




### Number of variants per indiv in ITP patients
variants_per_patient <- table(df_ITP_indivs$WES_500k_LoF_MAF01__Hauck_ITPnHAI) %>% as.data.frame() %>% 
  rename("No_variants" = Var1,
         "No_patients" = Freq)
variants_per_patient






# AIHA
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)




df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )



df_AIHA_indivs <- df %>% 
  filter(binary_AIHA == 1)


### Subset AIHA AND LOF
df_AIHA_and_lof <- df_AIHA_indivs %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0)



# Mean age AIHA AND LOF
mean_age_df_AIHA_and_lof <- (mean(df_AIHA_and_lof$AIHA_lower_age)) /365.25
mean_age_df_AIHA_and_lof



# SD AIHA AND LOF
sd_age_df_AIHA_and_lof <- (sd(df_AIHA_and_lof$AIHA_lower_age)) /365.25
sd_age_df_AIHA_and_lof





### Subset AIHA NOT LOF
df_AIHA_not_lof <- df_AIHA_indivs %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI == 0)

# Mean age AIHA NOT LOF
mean_age_df_AIHA_not_lof <- (mean(df_AIHA_not_lof$AIHA_lower_age)) /365.25
mean_age_df_AIHA_not_lof

# SD AIHA NOT LOF
sd_age_df_AIHA_not_lof <- (sd(df_AIHA_not_lof$AIHA_lower_age)) /365.25
sd_age_df_AIHA_not_lof




### t-test age comparison
t_test <- t.test(x = df_AIHA_and_lof$AIHA_lower_age, y = df_AIHA_not_lof$AIHA_lower_age)
t_test




### Number of variants per indiv in AIHA patients
variants_per_patient <- table(df_AIHA_indivs$WES_500k_LoF_MAF01__Hauck_ITPnHAI) %>% as.data.frame() %>% 
  rename("No_variants" = Var1,
         "No_patients" = Freq)
variants_per_patient











# _____003b Missense variants --------------------------------------------------

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)


main_table <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


df_bygenes <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/002_RareVariants_allscores_bygenes_20260414_no_HC.txt", select = c("Participant.ID", "Gene", "WES_500k_Miss5_MAF01"))

hauck_list <- read_xlsx("/home/stn/Documents/GitHub/Doc/analysis/genes_lists/IEI_Hauck_2024_FASLG.xlsx")




# Identifying genes associated with (ITP OR AIHA)  in the dataset

# Create new columns for further filtering
df_genes_associatied_col <- df_bygenes %>%
  mutate(itp_associated = "no",
         aiha_associated = "no" ) 

hauck_ITP <- hauck_list %>%
  filter(ITP != "") %>%
  select(ITP) %>%
  distinct(.) %>%
  pull()
hauck_ITP <- unlist(strsplit(hauck_ITP, "[;,]")) # Splits FASL;FASLG cell


hauck_AIHA <- hauck_list %>%
  filter(AIHA != "") %>%
  select(AIHA) %>%
  distinct(.) %>%
  pull()
hauck_AIHA <- unlist(strsplit(hauck_AIHA, "[;,]")) # Splits FASL;FASLG cell



# Annotate genes associated with ITP
df_genes_in_ITP <- df_genes_associatied_col %>%
  mutate(
    itp_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_ITP), "yes", "no")
    })
  )


# Annotate genes associated with AIHA
df_genes_in_AIHA <- df_genes_in_ITP %>%
  mutate(
    aiha_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_AIHA), "yes", "no")
    })
  )


# Keep only Hauck ITP or AIHA genes
df_genes_in_CAI <- df_genes_in_AIHA %>%
  filter(itp_associated == "yes" | aiha_associated == "yes")


# Annotate the sum of ITP or AIHA variants for each indiv
df_total_variants_per_indiv <- df_genes_in_CAI %>% 
  group_by(`Participant.ID`) %>% 
  mutate("total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA" = sum(WES_500k_Miss5_MAF01, na.rm = T)) %>% 
  ungroup() %>% 
  distinct(`Participant.ID`, .keep_all = T) %>% 
  select(-itp_associated, -aiha_associated, -Gene, -WES_500k_Miss5_MAF01) # distinct() will keep only the gene name of the first occurrence for each indiv. Remove the col to avoid confusion




# Annotate the sum of variants for each indiv in the main table
df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA <- left_join(main_table, df_total_variants_per_indiv, by = c("Participant.ID"))



# Replace NAs with zero: indivs without variant annotation for a gene should be interpreted as having zero variants in that gene.
df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA <- df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA %>% 
  mutate(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA = ifelse(is.na(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA), 0, total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA))


df_cai_indivs <- df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA %>% 
  filter(binary_ITP_or_AIHA_indivs == 1)




### Subset CAI AND MS
df_cai_and_MS <- df_cai_indivs %>% 
  filter(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA > 0)


# Mean age CAI AND MS
mean_age_df_cai_and_MS <- (mean(df_cai_and_MS$CAI_lower_age)) /365.25
mean_age_df_cai_and_MS



# SD CAI AND MS
sd_age_df_cai_and_MS <- (sd(df_cai_and_MS$CAI_lower_age)) /365.25
sd_age_df_cai_and_MS








### Subset CAI NOT MS
df_cai_not_MS <- df_cai_indivs %>% 
  filter(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA == 0)

# Mean age CAI NOT MS
mean_age_df_cai_not_MS <- (mean(df_cai_not_MS$CAI_lower_age)) /365.25
mean_age_df_cai_not_MS
# 64.82602


# SD CAI NOT MS
sd_age_df_cai_not_MS <- (sd(df_cai_not_MS$CAI_lower_age)) /365.25
sd_age_df_cai_not_MS
# 10.37886



### t-test age comparison
t_test <- t.test(x = df_cai_and_MS$CAI_lower_age, y = df_cai_not_MS$CAI_lower_age)
t_test




### Number of variants per indiv in CAI patients
cai_variants_per_patient <- table(df_cai_indivs$total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA) %>% as.data.frame() %>% 
  rename("No_variants" = Var1,
         "No_patients" = Freq)
cai_variants_per_patient








# ITP
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)


main_table <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


df_bygenes <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/002_RareVariants_allscores_bygenes_20260414_no_HC.txt", select = c("Participant.ID", "Gene", "WES_500k_Miss5_MAF01"))

hauck_list <- read_xlsx("/home/stn/Documents/GitHub/Doc/analysis/genes_lists/IEI_Hauck_2024_FASLG.xlsx")




# Identifying genes associated with (ITP OR AIHA)  in the dataset

# Create new columns for further filtering
df_genes_associatied_col <- df_bygenes %>%
  mutate(itp_associated = "no",
         aiha_associated = "no" ) 

hauck_ITP <- hauck_list %>%
  filter(ITP != "") %>%
  select(ITP) %>%
  distinct(.) %>%
  pull()
hauck_ITP <- unlist(strsplit(hauck_ITP, "[;,]")) # Splits FASL;FASLG cell


hauck_AIHA <- hauck_list %>%
  filter(AIHA != "") %>%
  select(AIHA) %>%
  distinct(.) %>%
  pull()
hauck_AIHA <- unlist(strsplit(hauck_AIHA, "[;,]")) # Splits FASL;FASLG cell



# Annotate genes associated with ITP
df_genes_in_ITP <- df_genes_associatied_col %>%
  mutate(
    itp_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_ITP), "yes", "no")
    })
  )


# Annotate genes associated with AIHA
df_genes_in_AIHA <- df_genes_in_ITP %>%
  mutate(
    aiha_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_AIHA), "yes", "no")
    })
  )


# Keep only Hauck ITP or AIHA genes
df_genes_in_CAI <- df_genes_in_AIHA %>%
  filter(itp_associated == "yes" | aiha_associated == "yes")


# Annotate the sum of ITP or AIHA variants for each indiv
df_total_variants_per_indiv <- df_genes_in_CAI %>% 
  group_by(`Participant.ID`) %>% 
  mutate("total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA" = sum(WES_500k_Miss5_MAF01, na.rm = T)) %>% 
  ungroup() %>% 
  distinct(`Participant.ID`, .keep_all = T) %>% 
  select(-itp_associated, -aiha_associated, -Gene, -WES_500k_Miss5_MAF01) # distinct() will keep only the gene name of the first occurrence for each indiv. Remove the col to avoid confusion




# Annotate the sum of variants for each indiv in the main table
df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA <- left_join(main_table, df_total_variants_per_indiv, by = c("Participant.ID"))



# Replace NAs with zero: indivs without variant annotation for a gene should be interpreted as having zero variants in that gene.
df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA <- df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA %>% 
  mutate(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA = ifelse(is.na(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA), 0, total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA))







# Subset ITP
df_ITP_indivs <- df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA %>% 
  filter(binary_ITP == 1)




### Subset ITP AND MS
df_ITP_and_MS <- df_ITP_indivs %>% 
  filter(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA > 0)


# Mean age ITP AND MS
mean_age_df_ITP_and_MS <- (mean(df_ITP_and_MS$TPI_lower_age)) /365.25
mean_age_df_ITP_and_MS



# SD ITP AND MS
sd_age_df_ITP_and_MS <- (sd(df_ITP_and_MS$TPI_lower_age)) /365.25
sd_age_df_ITP_and_MS








### Subset ITP NOT MS
df_ITP_not_MS <- df_ITP_indivs %>% 
  filter(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA == 0)

# Mean age ITP NOT MS
mean_age_df_cai_not_MS <- (mean(df_ITP_not_MS$TPI_lower_age)) /365.25
mean_age_df_cai_not_MS



# SD CAI NOT MS
sd_age_df_cai_not_MS <- (sd(df_ITP_not_MS$TPI_lower_age)) /365.25
sd_age_df_cai_not_MS




### t-test age comparison
t_test <- t.test(x = df_ITP_and_MS$TPI_lower_age, y = df_ITP_not_MS$TPI_lower_age)
t_test




### Number of variants per indiv in CAI patients
cai_variants_per_patient <- table(df_ITP_indivs$total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA) %>% as.data.frame() %>% 
  rename("No_variants" = Var1,
         "No_patients" = Freq)
cai_variants_per_patient







# AIHA
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)


main_table <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


df_bygenes <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/002_RareVariants_allscores_bygenes_20260414_no_HC.txt", select = c("Participant.ID", "Gene", "WES_500k_Miss5_MAF01"))

hauck_list <- read_xlsx("/home/stn/Documents/GitHub/Doc/analysis/genes_lists/IEI_Hauck_2024_FASLG.xlsx")




# Identifying genes associated with (ITP OR AIHA)  in the dataset

# Create new columns for further filtering
df_genes_associatied_col <- df_bygenes %>%
  mutate(itp_associated = "no",
         aiha_associated = "no" ) 

hauck_ITP <- hauck_list %>%
  filter(ITP != "") %>%
  select(ITP) %>%
  distinct(.) %>%
  pull()
hauck_ITP <- unlist(strsplit(hauck_ITP, "[;,]")) # Splits FASL;FASLG cell


hauck_AIHA <- hauck_list %>%
  filter(AIHA != "") %>%
  select(AIHA) %>%
  distinct(.) %>%
  pull()
hauck_AIHA <- unlist(strsplit(hauck_AIHA, "[;,]")) # Splits FASL;FASLG cell



# Annotate genes associated with ITP
df_genes_in_ITP <- df_genes_associatied_col %>%
  mutate(
    itp_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_ITP), "yes", "no")
    })
  )


# Annotate genes associated with AIHA
df_genes_in_AIHA <- df_genes_in_ITP %>%
  mutate(
    aiha_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_AIHA), "yes", "no")
    })
  )


# Keep only Hauck ITP or AIHA genes
df_genes_in_CAI <- df_genes_in_AIHA %>%
  filter(itp_associated == "yes" | aiha_associated == "yes")


# Annotate the sum of ITP or AIHA variants for each indiv
df_total_variants_per_indiv <- df_genes_in_CAI %>% 
  group_by(`Participant.ID`) %>% 
  mutate("total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA" = sum(WES_500k_Miss5_MAF01, na.rm = T)) %>% 
  ungroup() %>% 
  distinct(`Participant.ID`, .keep_all = T) %>% 
  select(-itp_associated, -aiha_associated, -Gene, -WES_500k_Miss5_MAF01) # distinct() will keep only the gene name of the first occurrence for each indiv. Remove the col to avoid confusion




# Annotate the sum of variants for each indiv in the main table
df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA <- left_join(main_table, df_total_variants_per_indiv, by = c("Participant.ID"))



# Replace NAs with zero: indivs without variant annotation for a gene should be interpreted as having zero variants in that gene.
df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA <- df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA %>% 
  mutate(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA = ifelse(is.na(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA), 0, total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA))







# Subset AIHA
df_AIHA_indivs <- df_total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA %>% 
  filter(binary_AIHA == 1)




### Subset AIHA AND MS
df_AIHA_and_MS <- df_AIHA_indivs %>% 
  filter(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA > 0)


# Mean age AIHA AND MS
mean_age_df_AIHA_and_MS <- (mean(df_AIHA_and_MS$AIHA_lower_age)) /365.25
mean_age_df_AIHA_and_MS



# SD ITP AND MS
sd_age_df_AIHA_and_MS <- (sd(df_AIHA_and_MS$AIHA_lower_age)) /365.25
sd_age_df_AIHA_and_MS








### Subset ITP NOT MS
df_AIHA_not_MS <- df_AIHA_indivs %>% 
  filter(total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA == 0)

# Mean age ITP NOT MS
mean_age_df_AIHA_not_MS <- (mean(df_AIHA_not_MS$AIHA_lower_age)) /365.25
mean_age_df_AIHA_not_MS



# SD CAI NOT MS
sd_age_df_AIHA_not_MS <- (sd(df_AIHA_not_MS$AIHA_lower_age)) /365.25
sd_age_df_AIHA_not_MS




### t-test age comparison
t_test <- t.test(x = df_AIHA_and_MS$AIHA_lower_age, y = df_AIHA_not_MS$AIHA_lower_age)
t_test




### Number of variants per indiv in patients
variants_per_patient <- table(df_AIHA_indivs$total_no_variants_WES_500k_Miss5_MAF01_ITPnAIHA) %>% as.data.frame() %>% 
  rename("No_variants" = Var1,
         "No_patients" = Freq)
variants_per_patient








# ____004 glm AlphaMS CAI/ITP/AIHA  ----------------------------------------------

# Supplemental table 2

# CAI

### MAF < 0.1, AlphaMS, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df_main <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


df_alphaMS <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260220_Summary_age_delays_diag_CH_PRS_REVEL_LOF_alphaMS.txt", 
                    select = c("Participant.ID","list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count")
)


df <- left_join(df_main, df_alphaMS, by = c("Participant.ID"))


# nas <- df %>% 
#   filter(is.na(list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count))
# 3573 indivs




# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








# ITP

### MAF < 0.1, AlphaMS, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df_main <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


df_alphaMS <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260220_Summary_age_delays_diag_CH_PRS_REVEL_LOF_alphaMS.txt", 
                    select = c("Participant.ID","list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count")
)


df <- left_join(df_main, df_alphaMS, by = c("Participant.ID"))


# nas <- df %>% 
#   filter(is.na(list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count))
# 3573 indivs




# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)






# AIHA

### MAF < 0.1, AlphaMS, TET2 included

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df_main <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


df_alphaMS <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260220_Summary_age_delays_diag_CH_PRS_REVEL_LOF_alphaMS.txt", 
                    select = c("Participant.ID","list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count")
)


df <- left_join(df_main, df_alphaMS, by = c("Participant.ID"))


# nas <- df %>% 
#   filter(is.na(list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count))
# 3573 indivs




# Subset if needed
df_subset <- df 
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  list_ITPnHAI_maf_MAF0.1_alphaMS.likely.pathogenic_count +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



















# ____005 AoU replication liste genes AHAI/ITP et maf <0.1% ---------------
# -Faire deux modèle de Cox et courbe de survie dans AoU: incidence cumulative depuis 
# 1) l’inclusion et 2) la naissance pour CAI (ITP+AHAI) en fonction des LOF



## Survival analysis
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)





df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260611_AoU_master_tab_ITP_AIHA_FUP_age_date.tsv") 



# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age,
                                     age_end_follow_up))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years,
                                      age_end_follow_up_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. <= 3.03 ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. > 3.03 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)






### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score,
  data=df_subset)


summary(cox_model)



# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ stratified_score,
                  data = df_subset)

fit_km



### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "outcome CAI x binary SLE PRS stratified genetic method.",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

# Add legend
legend("bottomright",
       legend = c("Non-carriers", "Carriers"),
       col = c("blue", "red"),
       lwd = 2)
###










# ____1.4 UKB LoF effect Cox since inclusion  --------------------------------------------------
# Logistic regression: effect of LOF variants in general population without SLE as term




# _______CAI -----------------------------------------------------------------

# Fig 1D

### CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )



# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Subset if needed
df_subset <- df_subset_years %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)



# Non-adjusted
cox_model_unadjusted <- coxph(
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI,
  data = df_subset
)

summary(cox_model_unadjusted)




# 
# ### Univariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI,
#   data=df_subset)
# 
# 
# summary(cox_model)
# 
# 
# 
# 
# 
# 
# ### Puis tracer incidence cumulative AIC en fonction de LOF temps 0= naissance. evenement: diag AIC

# First, create binary variable (carrier vs non-carrier)
df_subset$LoF_binary <- ifelse(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0, 1, 0)
# 
# 
# 
# plot(survfit(Surv(time = CAI_lower_age_years, time2 = age_end_follow_up_years, event = binary_ITP_or_AIHA_indivs) ~ LoF_binary, data = df_subset))
# 
# 
# 
# 
# # Kaplan-Meier curve (equivalent to Cox model with binary variable)
# fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ LoF_binary,
#                   data = df_subset)
# 
# fit_km
# 
# 
# 
# 
# 
# ### Graph formatting
# # Plot without X axis initially
# plot(fit_km, 
#      fun = "event",                # Cumulative incidence (1 - survival)
#      col = c("blue", "red"),
#      lwd = 2,
#      xlab = "Age (years)",
#      ylab = "Cumulative incidence",
#      main = "Cumulative incidence curve (Kaplan-Meier). LoF variants",
#      xaxt = "n")                   # Suppresses automatic X axis
# 
# # Add manual X axis with 5-year intervals
# max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
# breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)
# 
# axis(side = 1, at = breaks_x, labels = breaks_x)
# 
# # Add vertical grid at the breaks
# abline(v = breaks_x, col = "lightgray", lty = 2)
# 
# # Add legend
# legend("bottomright",
#        legend = c("Non-carriers", "Carriers"),
#        col = c("blue", "red"),
#        lwd = 2)
# ###
# 
# 
# 
# 
# 
# # Log-rank test to compare curves
# logrank_test <- survdiff(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ LoF_binary,
#                          data = df_subset)
# print(logrank_test)
# 
# 
# 
# 
# 
# 
# ### Save the plot as PNG
# setwd("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/009_first_article/001_cox_curves")
# png("cumulative_incidence_curve_km.png", width = 8, height = 6, units = "in", res = 300)
# 
# # Re-run the plotting code to save
# plot(fit_km, 
#      fun = "event",
#      col = c("blue", "red"),
#      lwd = 2,
#      xlab = "Age (years)",
#      ylab = "Cumulative incidence",
#      main = "Cumulative incidence curve (Kaplan-Meier). LoF variants",
#      xaxt = "n")
# 
# max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
# breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)
# 
# axis(side = 1, at = breaks_x, labels = breaks_x)
# abline(v = breaks_x, col = "lightgray", lty = 2)
# 
# legend("bottomright",
#        legend = c("Non-carriers", "Carriers"),
#        col = c("blue", "red"),
#        lwd = 2)
# 
# dev.off()
# 
# # Confirm save location
# message("Plot saved as: cumulative_incidence_curve_km.png in ", getwd())
# 




# Export data to create graph on graphpad prism

output <- df_subset %>% 
  select(LoF_binary, binary_ITP_or_AIHA_indivs, follow_up_duration_years)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/001_Fig_1p4_incidence_cumulative_naissance_LOF_CAI.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')




# Export data to create graph on graphpad prism
# Data to use age of inclusion instead of birth

# age_inclusion = as.numeric(as.Date(blood_draw_date) - as.Date(formatted_DATE_BIRTH))) # 

output <- df_subset %>% 
  select(LoF_binary, binary_ITP_or_AIHA_indivs, follow_up_duration_years)






# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/001_incidence_cumulative_naissance_LOF_CAI_for_prism.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')














# _______ (legacy) ITP -----------------------------------------------------------------
# Legacy : results not used. Only CAI results was used

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP == 1, 
                                     TPI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    TPI_lower_age_years = TPI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP == 1,
                                      TPI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Subset if needed
df_subset <- df_subset_years %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)






### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI,
  data=df_subset)


summary(cox_model)





### Puis tracer incidence cumulative AIC en fonction de LOF temps 0= naissance. evenement: diag AIC
# First, create binary variable (carrier vs non-carrier)
df_subset$LoF_binary <- ifelse(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0, 1, 0)



plot(survfit(Surv(time = TPI_lower_age_years, time2 = age_end_follow_up_years, event = binary_ITP) ~ LoF_binary, data = df_subset))




# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP) ~ LoF_binary,
                  data = df_subset)

fit_km





### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Age (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). LoF variants",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

# Add legend
legend("bottomright",
       legend = c("Non-carriers", "Carriers"),
       col = c("blue", "red"),
       lwd = 2)
###









# _______ (legacy) AIHA -----------------------------------------------------------------
# Legacy : results not used. Only CAI results was used

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_AIHA == 1, 
                                     AIHA_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    AIHA_lower_age_years = AIHA_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_AIHA == 1,
                                      AIHA_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Subset if needed
df_subset <- df_subset_years %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_AIHA) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)






### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_AIHA) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI,
  data=df_subset)


summary(cox_model)



























#___2) PRS SLE to predict risk of CAI ----------------------------------------
# Supplemental table 4

# Faire avec AIHA et ITP séparé 



# ______2.1 CAI---------------------------------------------------------------------
# mkdir -p $UKBDATA/027_first_article_OptSurvCutR/001_section_2
# cd $UKBDATA/027_first_article_OptSurvCutR/001_section_2


# Fig 2B


### CAI


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# 
# # Subset if needed
# df_subset <- df
# nrow(df_subset)
# 
# 
# summary(glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
#     omop_gender_concept_id +
#     age_end_follow_up  + 
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset, family = binomial)
#   
# )
# 
# 
# 
# 
# model_logistic_regression <- glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
#     omop_gender_concept_id +
#     age_end_follow_up  + 
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset, family = binomial)
# 
# 
# results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
#   mutate(
#     OR = exp(estimate),
#     CI_inf = exp(conf.low),
#     CI_sup = exp(conf.high)
#   ) %>%
#   select(term, OR, CI_inf, CI_sup, p.value, statistic)
# 
# print(results, digits = 3)
# 
# 
# 
# summary_model <- summary(model_logistic_regression)
# 
# # Extract coefficients with p-values
# coef_table <- summary_model$coefficients
# 
# # View the full table with p-values (exact)
# print(coef_table, digits = 20)



# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. <= 3.03 ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. > 3.03 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)



# Non-adjusted
cox_model_unadjusted <- coxph(
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score,
  data = df_subset
)

summary(cox_model_unadjusted)





### Export pfor Prism
output <- df_subset %>% 
  mutate(binary_SLE_PRS_stratified_3p03 = stratified_score) %>% 
  select(binary_SLE_PRS_stratified_3p03, binary_ITP_or_AIHA_indivs, follow_up_duration_years)



# Fig 2B

# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/articles/first_year_phd_pubs/002_first_article/para_julie/Fig2B_cumulative_incidence_SLE_PGS_CAI.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')











# _________ Adjusted for LOF --------------------------------------------------


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)


model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)





# _________ Exclusion of ITP/AIHA + SLE --------------------------------------------------

# Supplemental table 4

df_binary_CAI_and_SLE <- df |> 
  mutate(binary_CAI_and_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & binary_SLE == 1 ~ "1",
    TRUE ~ "0"
  ))




# Subset if needed
df_subset <- df_binary_CAI_and_SLE |> 
  filter(binary_CAI_and_SLE == 0)
nrow(df_subset)



summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)

















### Defining cutpoint


# SLURM: Genetic method

cat << \EOF > $UKBSCRIPTS/prs_cutpoints_gen.sh
#!/bin/bash
#SBATCH --job-name=opt_gen_prscai
#SBATCH --output=opt_gen_prscai_%j.out
#SBATCH --error=opt_gen_prscai_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=8G




# Carregar módulos
module load r/4.4.0

# Configurações de performance
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_MAX_NUM_DLLS=1000
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Executar com logging detalhado
Rscript --vanilla --verbose -e '
library(dplyr)
library(OptSurvCutR)
library(survival)

# Configurar opções do R
options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(42)

# Carregar dados
cat("Carregando dados...\n")


# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/lustre09/project/6097258/stnfaria/ukb_rare_variants/000_whole_cohort/data/027_first_article_OptSurvCutR/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )


# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age,
                                     age_end_follow_up))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years,
                                      age_end_follow_up_years)
  )



# Subset if needed
df_subset <- df_subset_years %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)


df_subset_cutpoints <- df_subset_years %>% 
  select(Standard.PRS.for.systemic.lupus.erythematosus..SLE., follow_up_duration, binary_ITP_or_AIHA_indivs) %>% 
  filter(follow_up_duration > 0) 




# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "Standard.PRS.for.systemic.lupus.erythematosus..SLE.",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 1,
  method = "genetic",   # MORE STABLE
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result

'


EOF



sbatch $UKBSCRIPTS/prs_cutpoints_gen.sh


# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: Standard.PRS.for.systemic.lupus.erythematosus..SLE.
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 42.6966
# ✔ Recommended Cut-point(s): 3.03



















# _________AoU PGS000196 ---------------------------------------------------------------------


# ______________ Cutpoint -------------------------------------------------



# Core dependencies
install.packages(c("remotes", "survival"))

# Optional but highly recommended for multi-cut genetic optimisation
install.packages("rgenoud")

# Install the package from GitHub
remotes::install_github("paytonyau/OptSurvCutR")


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

library(OptSurvCutR)



### Defining cutpoint


df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 


# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age,
                                     age_end_follow_up))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years,
                                      age_end_follow_up_years)
  )



# Subset if needed
df_subset_cutpoints <- df_subset_years %>% 
  select(PGS000196, follow_up_duration, binary_ITP_or_AIHA_indivs) %>% 
  filter(follow_up_duration > 0) 




# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "PGS000196",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 1,
  method = "genetic",   # MORE STABLE
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result



# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: PGS000196
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 42.1626
# ✔ Recommended Cut-point(s): 2.573
# Hint: Use `summary()` for clinical details and Cox regression.





# ______________ Cox with cutpoint -------------------------------------------------



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    PGS000196 <= 2.573 ~ "0",
    PGS000196 > 2.573 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)







### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)




# _________________ Export for prism -------------------------------------------------





output <- df_subset %>% 
  mutate(binary_PGS000196_2p573 = stratified_score) %>% 
  select(binary_PGS000196_2p573, binary_ITP_or_AIHA_indivs, follow_up_duration_years)




# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/002_data_for_julie_fig2_code_2p1p2p3/001_Fig2A_AoU_CAI_PGS000196.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')



# ______________ Cox non-binary -------------------------------------------------



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    PGS000196 <= 2.573 ~ "0",
    PGS000196 > 2.573 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)







### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ PGS000196 +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)










# ______________ Logistic regression------------------------------------------------------

### CAI

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










#### Export for Prism

# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. <= 3.03 ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. > 3.03 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)




output <- df_subset %>% 
  mutate(binary_SLE_PRS_stratified_3p03 = stratified_score) %>% 
  select(binary_SLE_PRS_stratified_3p03, binary_ITP_or_AIHA_indivs, follow_up_duration_years, age_inclusion)




# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/002_incidence_cumulative_naissance_SLE_CAI_for_prism.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')











# _________AoU PGS004917  ---------------------------------------------------------------------


# ______________ Cutpoint -------------------------------------------------


### Install package
# Core dependencies
install.packages(c("remotes", "survival"))

# Optional but highly recommended for multi-cut genetic optimisation
install.packages("rgenoud")

# Install the package from GitHub
remotes::install_github("paytonyau/OptSurvCutR")




### Run the program
required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

library(OptSurvCutR)



### Defining cutpoint


df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 


# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age,
                                     age_end_follow_up))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years,
                                      age_end_follow_up_years)
  )



# Subset if needed
df_subset_cutpoints <- df_subset_years %>% 
  select(PGS004917, follow_up_duration, binary_ITP_or_AIHA_indivs) %>% 
  filter(follow_up_duration > 0) 




# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "PGS004917",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 1,
  method = "genetic",   # MORE STABLE
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: PGS004917
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 55.3413
# ✔ Recommended Cut-point(s): 17.593
# Hint: Use `summary()` for clinical details and Cox regression.






# ______________ Cox with cutpoint -------------------------------------------------



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    PGS004917 <= 17.593 ~ "0",
    PGS004917 > 17.593 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)







### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)





# _________________ Export for prism -------------------------------------------------





output <- df_subset %>% 
  mutate(binary_PGS004917_17p593 = stratified_score) %>% 
  select(binary_PGS004917_17p593, binary_ITP_or_AIHA_indivs, follow_up_duration_years)




# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/002_data_for_julie_fig2_code_2p1p2p3/002_Fig2A_AoU_CAI_PGS004917.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')





# ______________ Cox non-binary -------------------------------------------------




required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    PGS004917 <= 17.593 ~ "0",
    PGS004917 > 17.593 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)







### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ PGS004917 +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)





# ______________ Logistic regression------------------------------------------------------


### CAI

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)













# _________AoU PGS000196 adjusted for LoF ---------------------------------------------------------------------



### CAI

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS000196 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS000196 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








# _________AoU PGS004917 adjusted for LoF ---------------------------------------------------------------------



### CAI

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS004917 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS004917 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







# _________AoU PGS000196 after exclusion of patiens with ITP/AIHA + SLE---------------------------------------------------------------------



### CAI

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





df_binary_CAI_and_SLE <- df |> 
  mutate(binary_CAI_and_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & binary_SLE == 1 ~ "1",
    TRUE ~ "0"
  ))




# Subset if needed
df_subset <- df_binary_CAI_and_SLE |> 
  filter(binary_CAI_and_SLE == 0)
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







# _________AoU PGS004917 after exclusion of patiens with ITP/AIHA + SLE---------------------------------------------------------------------



### CAI

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




df_binary_CAI_and_SLE <- df |> 
  mutate(binary_CAI_and_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & binary_SLE == 1 ~ "1",
    TRUE ~ "0"
  ))



# Subset if needed
df_subset <- df_binary_CAI_and_SLE |> 
  filter(binary_CAI_and_SLE == 0)
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)


















# ______2.2 AIHA---------------------------------------------------------------------

### AIHA

# Logistic regression


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)




summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)





# _________Adjusted for LOF --------------------------------------------------


summary(glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)


model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)






# _________ Exclusion of ITP/AIHA + SLE --------------------------------------------------


df_binary_CAI_and_SLE <- df |> 
  mutate(binary_CAI_and_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & binary_SLE == 1 ~ "1",
    TRUE ~ "0"
  ))




# Subset if needed
df_subset <- df_binary_CAI_and_SLE |> 
  filter(binary_CAI_and_SLE == 0)
nrow(df_subset)



summary(glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)






# _________AoU PGS000196 ---------------------------------------------------------------------



### AIHA

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)






# ______________ Cox with cutpoint -------------------------------------------------



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_AIHA == 1, 
                                     AIHA_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    lower_age_years = AIHA_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_AIHA == 1,
                                      lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    PGS000196 <= 2.573 ~ "0",
    PGS000196 > 2.573 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)







### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_AIHA) ~ stratified_score +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)




# _________________ Export for prism -------------------------------------------------



output <- df_subset %>% 
  mutate(binary_PGS000196_2p573 = stratified_score) %>% 
  select(binary_PGS000196_2p573, binary_AIHA, follow_up_duration_years)




# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/002_data_for_julie_fig2_code_2p1p2p3/005_Fig2A_AoU_AIHA_PGS000196.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')








# _________AoU PGS004917  ---------------------------------------------------------------------


### AIHA

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)









# ______________ Cox with cutpoint -------------------------------------------------



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_AIHA == 1, 
                                     AIHA_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    lower_age_years = AIHA_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_AIHA == 1,
                                      lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    PGS004917 <= 17.593 ~ "0",
    PGS004917 > 17.593 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)







### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_AIHA) ~ stratified_score +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)





# _________________ Export for prism -------------------------------------------------





output <- df_subset %>% 
  mutate(binary_PGS004917_17p593 = stratified_score) %>% 
  select(binary_PGS004917_17p593, binary_AIHA, follow_up_duration_years)




# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/002_data_for_julie_fig2_code_2p1p2p3/006_Fig2A_AoU_AIHA_PGS004917.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')








# _________AoU PGS000196 adjusted for LoF ---------------------------------------------------------------------



### AIHA

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  PGS000196 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  PGS000196 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








# _________AoU PGS004917 adjusted for LoF ---------------------------------------------------------------------



### AIHA

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  PGS004917 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  PGS004917 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







# _________AoU PGS000196 after exclusion of patiens with ITP/AIHA + SLE---------------------------------------------------------------------



### AIHA

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 


df_binary_CAI_and_SLE <- df |> 
  mutate(binary_CAI_and_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & binary_SLE == 1 ~ "1",
    TRUE ~ "0"
  ))




# Subset if needed
df_subset <- df_binary_CAI_and_SLE |> 
  filter(binary_CAI_and_SLE == 0)
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







# _________AoU PGS004917 after exclusion of patiens with SLE---------------------------------------------------------------------



### AIHA

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





df_binary_CAI_and_SLE <- df |> 
  mutate(binary_CAI_and_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & binary_SLE == 1 ~ "1",
    TRUE ~ "0"
  ))




# Subset if needed
df_subset <- df_binary_CAI_and_SLE |> 
  filter(binary_CAI_and_SLE == 0)
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










# ______2.3 ITP---------------------------------------------------------------------



### ITP


# Logistic regression

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)





summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)





# _________Adjusted for LOF --------------------------------------------------


summary(glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)


model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)





# _________ Exclusion of ITP/AIHA + SLE --------------------------------------------------


df_binary_CAI_and_SLE <- df |> 
  mutate(binary_CAI_and_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & binary_SLE == 1 ~ "1",
    TRUE ~ "0"
  ))




# Subset if needed
df_subset <- df_binary_CAI_and_SLE |> 
  filter(binary_CAI_and_SLE == 0)
nrow(df_subset)



summary(glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)






# _________AoU PGS000196 ---------------------------------------------------------------------



### ITP

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)





# ______________ Cox with cutpoint -------------------------------------------------



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP == 1, 
                                     TPI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    lower_age_years = TPI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP == 1,
                                      lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    PGS000196 <= 2.573 ~ "0",
    PGS000196 > 2.573 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)







### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP) ~ stratified_score +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)




# _________________ Export for prism -------------------------------------------------



output <- df_subset %>% 
  mutate(binary_PGS000196_2p573 = stratified_score) %>% 
  select(binary_PGS000196_2p573, binary_ITP, follow_up_duration_years)




# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/002_data_for_julie_fig2_code_2p1p2p3/003_Fig2A_AoU_ITP_PGS000196.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')









# _________AoU PGS004917  ---------------------------------------------------------------------


### ITP

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)






# ______________ Cox with cutpoint -------------------------------------------------



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP == 1, 
                                     TPI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))


# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    lower_age_years = TPI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP == 1,
                                      lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years)
  )



# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    PGS004917 <= 17.593 ~ "0",
    PGS004917 > 17.593 ~"1"
  ))



# Subset if needed
df_subset <- df_stratified_score %>% 
  filter(follow_up_duration > 0)
nrow(df_subset)







### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP) ~ stratified_score +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)





# _________________ Export for prism -------------------------------------------------





output <- df_subset %>% 
  mutate(binary_PGS004917_17p593 = stratified_score) %>% 
  select(binary_PGS004917_17p593, binary_ITP, follow_up_duration_years)




# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/002_data_for_julie_fig2_code_2p1p2p3/004_Fig2A_AoU_ITP_PGS004917.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')










# _________AoU PGS000196 adjusted for LoF ---------------------------------------------------------------------



### ITP

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  PGS000196 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  PGS000196 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








# _________AoU PGS004917 adjusted for LoF ---------------------------------------------------------------------



### ITP

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  PGS004917 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  PGS004917 +
    AoU_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







# _________AoU PGS000196 after exclusion of patiens with ITP/AIHA + SLE---------------------------------------------------------------------



### ITP

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 



df_binary_CAI_and_SLE <- df |> 
  mutate(binary_CAI_and_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & binary_SLE == 1 ~ "1",
    TRUE ~ "0"
  ))


# Subset if needed
df_subset <- df_binary_CAI_and_SLE |> 
  filter(binary_CAI_and_SLE == 0)
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







# _________AoU PGS004917 after exclusion of patiens with SLE---------------------------------------------------------------------



### ITP

# Logistic regression


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





df_binary_CAI_and_SLE <- df |> 
  mutate(binary_CAI_and_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & binary_SLE == 1 ~ "1",
    TRUE ~ "0"
  ))


# Subset if needed
df_subset <- df_binary_CAI_and_SLE |> 
  filter(binary_CAI_and_SLE == 0)
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)












# ______2.4 Assoc PRS age au diag -------------------------------------------
# 3.  @Stennio Faria: Parmi les patients avec AIC, association entre PRS SLE et âge au diag

### CAI

# Logistic regression
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )









# Subset if needed
df_subset <- df %>% 
  filter(binary_ITP_or_AIHA_indivs == 1)
nrow(df_subset)


summary(glm(
  as.factor(CAI_lower_age) ~  scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(CAI_lower_age) ~  scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







# ______2.5 PRS x age au diag -------------------------------------------
# Analysis: subset CAI binary = 1,
# exlcure indivs with lymphoiod malignancy avant CAI AND 
# exclure indivs qui ont Lof AND
# exclure ceux avec CH, refaire age x PRS



### CAI

# Logistic regression
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )



df_lymphoid_malignancies_before_CAI <- df %>% 
  mutate(binary_lymphoid_malignancies_before_CAI = case_when(
    !is.na(CAI_lower_age) & LYMPHOME_lower_age < CAI_lower_age ~ "1",
    !is.na(CAI_lower_age) & LLC_lower_age < CAI_lower_age  ~ "1",
    !is.na(CAI_lower_age) & Other_LEUK_lower_age < CAI_lower_age  ~ "1",
    TRUE ~ "0"
  )) 




# Subset if needed
df_subset <- df_lymphoid_malignancies_before_CAI %>% 
  
  filter(binary_ITP_or_AIHA_indivs == 1 & # subset CAI binary = 1
           binary_lymphoid_malignancies_before_CAI == 0 & # exlcure indivs with lymphoiod malignancy avant CAI AND
           WES_500k_LoF_MAF01__Hauck_ITPnHAI == 0 & # exclure indivs qui ont Lof AND
           is.na(CH_Age_sample) # exclure ceux avec CH
         
  )



nrow(df_subset)


summary(glm(
  as.factor(CAI_lower_age) ~  scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(CAI_lower_age) ~  scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)


















# ______ (other paper) 2.6 risk of CAI in SLE patients -------------------------------------


#### To put in a separate paper



# subset SLE, binary_CAI


# Logistic regression

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df %>% 
  filter(!is.na(SLE))
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










### subset SLE, binary_ITP

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df %>% 
  filter(!is.na(SLE))
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)









### subset SLE, binary_AIHA

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df %>% 
  filter(!is.na(SLE))
nrow(df_subset)


summary(glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










# LOF + REVEL 0.9 variants
# Variants LoF not HC
# No PRS score
# LoF_or_REVEL0.9 column created by the sum of LOF not HC column + REVEL0.9
# No risk of overlapping variants, since REVEL are missense


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Create LoF_not_CH_revel_0.9 column
annotated_table <- df  %>% 
  mutate("LoF_not_CH_revel_0.9" = WES_500k_LoF_MAF01__Hauck_ITPnHAI + `list_ITPnHAI_maf_MAF01_revel_0.9_count`)


# Subset if needed
df_subset <- annotated_table %>% 
  filter(!is.na(SLE))
nrow(df_subset)



# -Le risque de CAI chez les patients avec SLE. Ton modele doit être
# ITP_or_AIHA_indivs ~ LOF + PRS + sexe + age + PC1-10 dans le sous groupe avec SLE
summary(glm(as.factor(binary_ITP_or_AIHA_indivs) ~ LoF_not_CH_revel_0.9 + age_end_follow_up + omop_gender_concept_id + 
              Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
              Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
              Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
              Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
              Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
            data=df_subset, family = binomial))


sle_cai_logistic_regression <- glm(as.factor(binary_ITP_or_AIHA_indivs) ~ LoF_not_CH_revel_0.9 + age_end_follow_up + omop_gender_concept_id + 
                                     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
                                     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
                                     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
                                     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
                                     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
                                   data=df_subset, family = binomial)




results <- tidy(sle_cai_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










# ______2.7 Adjusted for LOF ---------------------------------------------------------------------
# mkdir -p $UKBDATA/027_first_article_OptSurvCutR/001_section_2
# cd $UKBDATA/027_first_article_OptSurvCutR/001_section_2

### CAI

# Logistic regression
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# Subset if needed
df_subset <- df
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)

model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)




### ITP
summary(glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)

model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)




### AIHA
summary(glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)
  
)

model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)












# ___3) PREDICT SLE AMONG AIC ------------------------------------


# ____3.1 subset to !SLE (any diag age) --------------------------------------------


### This analysis probably shjould be in 2), since it's done in the whole cohort and
# section 3) concerns CAI patients subset




# Logistic regression

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )





# Subset if needed
df_subset <- df %>% 
  filter(is.na(SLE))
nrow(df_subset)






### CAI
summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)





### ITP
summary(glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)






### AIHA
summary(glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)













# ____3.3 subset CAI glm SLE ~ PGS SLE --------------------------------------------
#Canvas 3.3
# @Stennio Faria: patients avec AIC en ayant enlevé ceux avec SLE avant AIC.


# Faire modele de Cox multivarié
# 
# SLE ~ PGS_SLE + sex_label + age inclusion + 
#     PC1 + PC2 + PC3 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + 
#     PC10, data = subset(UKB, AIC =="1"),
# 
# Le temps 0 est le diagnostic de l’AIC
# Puis séparer en PGS_SLE en 2 groupes par breakpoint defini plus faut. Comparaison des deux courbes par Cox univarié 




# Supplemental table 5

# Fig 2D


# Incidence cumulative

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)



df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )



# Outcome: SLE instead of CAI 

# Annotate indivs SLE < CAI
df_binary_SLE_before_CAI <- df %>% 
  mutate(binary_SLE_before_CAI = case_when(
    !is.na(SLE) & SLE_lower_age < CAI_lower_age ~ 1,
    TRUE ~ 0
  )) 



# Annotate CAI before SLE
df_binary_SLE_before_CAI <- df_binary_SLE_before_CAI |> 
  mutate(binary_CAI_before_SLE = case_when(
    binary_ITP_or_AIHA_indivs == 1 & CAI_lower_age < SLE_lower_age ~ 1,
    TRUE ~ 0 
  ))









### General population



# Subset if needed
df_subset <- df_binary_SLE_before_CAI 
nrow(df_subset)


# test
sle_indivs_after_filter <- df_subset |> 
  filter(!is.na(SLE)) |> 
  nrow()




summary(glm(
  binary_SLE ~  scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    age_end_follow_up +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset, family = binomial)
  
)



model_logistic_regression <- glm(
  binary_SLE ~  scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) +
    age_end_follow_up +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset, family = binomial)





results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)









### CAI indivs

# Subset if needed
df_subset <- df_binary_SLE_before_CAI %>% 
  #  filter(binary_CAI_before_SLE != 1) %>% 
  filter(binary_ITP_or_AIHA_indivs == 1)

nrow(df_subset)


# test
sle_indivs_after_filter <- df_subset |> 
  filter(!is.na(SLE)) |> 
  nrow()


# Fig 2D

summary(glm(
  as.numeric(binary_SLE) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    age_end_follow_up +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset, family = binomial)
  
)



model_logistic_regression <- glm(
  as.numeric(binary_SLE) ~  Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    age_end_follow_up +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)











#### Export for Prism

output <- df_subset %>% 
  select(Standard.PRS.for.systemic.lupus.erythematosus..SLE., binary_SLE, follow_up_duration_years)




# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/003_inc_cumul_section3_SLE_CAI_1SD_for_prism.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')










# _______AoU replication --------------------------------------------
# Figure 2E: Il faut rajouter les résultats de AoU. 
# Tu peux refaire les analyses toi avec la date de fin de suivi 
# (régression logistique SLE ~ PGS SLE + age derniere nouvelle/sexe/PC chez les patients avec ITP/AHAI). 




# ___________ PGS000196 ---------------------------------------------------------------------



# CAI
required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom", 
                       "ggplot2", 
                       "tidyplots",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}



df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 



# Annotate 
df_binary_SLE_before_CAI <- df |> 
  mutate(binary_SLE_before_CAI = case_when(
    binary_ITP_or_AIHA_indivs == 1 & SLE_lower_age < CAI_lower_age ~ 1,
    TRUE ~ 0 
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_SLE_before_CAI %>%
  mutate(follow_up_duration = ifelse(binary_SLE == 1,
                                     SLE_lower_age - CAI_lower_age,
                                     age_end_follow_up - CAI_lower_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    SLE_lower_age_years = SLE_lower_age / 365.25,
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_SLE == 1,
                                      SLE_lower_age_years - CAI_lower_age_years,
                                      age_end_follow_up_years - CAI_lower_age_years
    ))


### General population


# Subset if needed
df_subset <- df_binary_SLE_before_CAI

nrow(df_subset)




### Logistic regression


summary(glm(
  as.factor(binary_SLE) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)


model_logistic_regression <- glm(
  as.factor(binary_SLE) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)





### CAI indivs

# Subset if needed
df_subset <- df_binary_SLE_before_CAI %>% 
  #  filter(binary_CAI_before_SLE != 1) %>% # Eliminate SLE before CAI
  filter(binary_ITP_or_AIHA_indivs == 1) # Only CAI patients
nrow(df_subset)




### Logistic regression

summary(glm(
  as.factor(binary_SLE) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)


model_logistic_regression <- glm(
  as.factor(binary_SLE) ~  PGS000196 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)






# ______________Cox subset CAI patients -------------------------------------------------------------------

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Outcome: SLE instead of CAI 


# Annotate 
df_binary_SLE_before_CAI <- df |> 
  mutate(binary_SLE_before_CAI = case_when(
    binary_ITP_or_AIHA_indivs == 1 & SLE_lower_age < CAI_lower_age ~ 1,
    TRUE ~ 0 
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_SLE_before_CAI %>%
  mutate(follow_up_duration = ifelse(binary_SLE == 1,
                                     SLE_lower_age - CAI_lower_age,
                                     age_end_follow_up - CAI_lower_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    SLE_lower_age_years = SLE_lower_age / 365.25,
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_SLE == 1,
                                      SLE_lower_age_years - CAI_lower_age_years,
                                      age_end_follow_up_years - CAI_lower_age_years
    ))





# Cox not used: AoU doesn't have enough follow-up duration for Cox. Using logistic regression instead



# Subset if needed
df_subset <- df_subset_years %>%
  filter(follow_up_duration > 0) %>%
  filter(binary_SLE_before_CAI == 0) %>% # Eliminate CAI before SLE
  filter(binary_ITP_or_AIHA_indivs == 1)
nrow(df_subset)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_SLE) ~ PGS000196 +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)










# ___________PGS004917 ---------------------------------------------------------------------



# CAI
required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom", 
                       "ggplot2", 
                       "tidyplots",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}



df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate 
df_binary_SLE_before_CAI <- df |> 
  mutate(binary_SLE_before_CAI = case_when(
    binary_ITP_or_AIHA_indivs == 1 & SLE_lower_age < CAI_lower_age ~ 1,
    TRUE ~ 0 
  ))





### General population


# Subset if needed
df_subset <- df_binary_SLE_before_CAI

nrow(df_subset)




### Logistic regression

summary(glm(
  as.factor(binary_SLE) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)


model_logistic_regression <- glm(
  as.factor(binary_SLE) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)








### CAI indivs

# Subset if needed
df_subset <- df_binary_SLE_before_CAI %>% 
  #  filter(binary_CAI_before_SLE != 1) %>% # Eliminate SLE before CAI
  filter(binary_ITP_or_AIHA_indivs == 1) # Only CAI patients
nrow(df_subset)




### Logistic regression

summary(glm(
  as.factor(binary_SLE) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)


model_logistic_regression <- glm(
  as.factor(binary_SLE) ~  PGS004917 +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)







#### Export for Prism


# PGS000196
output <- df_subset %>% 
  select(PGS000196, binary_SLE, follow_up_duration_years)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/001_3p3_figure2E_AoU_PGS000196.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# PGS004917
output <- df_subset %>% 
  select(PGS004917, binary_SLE, follow_up_duration_years)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/001_3p3_figure2E_AoU_PGS004917.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')







# ______________Cox subset CAI patients -------------------------------------------------------------------

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    AoU_LoF_MAF1__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Outcome: SLE instead of CAI 


# Annotate CAI before SLE
df_binary_SLE_before_CAI <- df |> 
  mutate(binary_SLE_before_CAI = case_when(
    binary_ITP_or_AIHA_indivs == 1 & SLE_lower_age < CAI_lower_age ~ 1,
    TRUE ~ 0 
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_SLE_before_CAI %>%
  mutate(follow_up_duration = ifelse(binary_SLE == 1,
                                     SLE_lower_age - CAI_lower_age,
                                     age_end_follow_up - CAI_lower_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    SLE_lower_age_years = SLE_lower_age / 365.25,
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_SLE == 1,
                                      SLE_lower_age_years - CAI_lower_age_years,
                                      age_end_follow_up_years - CAI_lower_age_years
    ))





# Cox not used: AoU doesn't have enough follow-up duration for Cox. Using logistic regression instead



# Subset if needed
df_subset <- df_subset_years %>%
  filter(follow_up_duration > 0) %>%
  filter(binary_SLE_before_CAI == 0) %>% # Eliminate CAI before SLE
  filter(binary_ITP_or_AIHA_indivs == 1)
nrow(df_subset)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_SLE) ~ PGS004917 +
    first_condition_age +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10,
  data=df_subset)


summary(cox_model)












# _______3.3a def 1 cutpoint UKB ------------------------------------------------
# Cutpoints for the SLE PRS in section 3)
# Outcome: SLE
# Follow-up: CAI until SLE
# Subset: CAI indivs
# 2 cutpoints don't perform well

# Manuscript: section PGS for SLE



# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: Standard.PRS.for.systemic.lupus.erythematosus..SLE.
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 15.7532
# ✔ Recommended Cut-point(s): 1.662




# Genetic method


# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/lustre09/project/6097258/stnfaria/ukb_rare_variants/000_whole_cohort/data/027_first_article_OptSurvCutR/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )





# Outcome: SLE instead of CAI 

# Annotate indivs SLE < CAI
df_binary_SLE_before_CAI <- df %>% 
  mutate(binary_SLE_before_CAI = case_when(
    
    is.na(SLE) | is.na(CAI_lower_age) ~ 0,
    SLE_lower_age < CAI_lower_age ~ 1,
    TRUE ~ 0
    
  )) 
# 13 indivs with SLE before CAI


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_SLE_before_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_SLE == 1, 
                                     SLE_lower_age - CAI_lower_age,
                                     age_end_follow_up - CAI_lower_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    SLE_lower_age_years = SLE_lower_age / 365.25,
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_SLE == 1,
                                      SLE_lower_age_years - CAI_lower_age_years,
                                      age_end_follow_up_years - CAI_lower_age_years
    ))




# Subset if needed
df_subset_cutpoints <- df_subset_years %>% 
  filter(follow_up_duration > 0) %>% 
  filter(binary_SLE_before_CAI != 1) %>% 
  filter(binary_ITP_or_AIHA_indivs == 1) %>% 
  select(Standard.PRS.for.systemic.lupus.erythematosus..SLE., follow_up_duration, binary_SLE)
nrow(df_subset_cutpoints)



# genetic method
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "Standard.PRS.for.systemic.lupus.erythematosus..SLE.",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_SLE",
  num_cuts = 1,
  method = "genetic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result



# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: Standard.PRS.for.systemic.lupus.erythematosus..SLE.
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 15.7532
# ✔ Recommended Cut-point(s): 1.662








# _______3.3b cutpoint PGS000196 ------------------------------------------------
# Manuscript: section PGS for SLE


#── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: PGS000196
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 24.0283
# ✔ Recommended Cut-point(s): 3.683
# Hint: Use `summary()` for clinical details and Cox regression.



# Core dependencies
install.packages(c("remotes", "survival"))

# Optional but highly recommended for multi-cut genetic optimisation
install.packages("rgenoud")

# Install the package from GitHub
remotes::install_github("paytonyau/OptSurvCutR")


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

library(OptSurvCutR)




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Outcome: SLE instead of CAI 

# Annotate indivs SLE < CAI
df_binary_SLE_before_CAI <- df %>% 
  mutate(binary_SLE_before_CAI = case_when(
    is.na(SLE) & SLE_lower_age < CAI_lower_age ~ 1,
    TRUE ~ 0
  )) 



# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_SLE_before_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_SLE == 1, 
                                     SLE_lower_age - CAI_lower_age,
                                     age_end_follow_up - CAI_lower_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    SLE_lower_age_years = SLE_lower_age / 365.25,
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_SLE == 1,
                                      SLE_lower_age_years - CAI_lower_age_years,
                                      age_end_follow_up_years - CAI_lower_age_years
    ))




# Subset if needed
df_subset_cutpoints <- df_subset_years %>% 
  filter(follow_up_duration > 0) %>% 
  filter(binary_SLE_before_CAI != 1) %>% 
  filter(binary_ITP_or_AIHA_indivs == 1) %>% 
  select(PGS000196, follow_up_duration, binary_SLE)
nrow(df_subset_cutpoints)



# genetic method
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "PGS000196",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_SLE",
  num_cuts = 1,
  method = "genetic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result








# _______3.3b cutpoint PGS004917 ------------------------------------------------
# Manuscript: section PGS for SLE


# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: PGS004917
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 18.6301
# ✔ Recommended Cut-point(s): 18.139
# Hint: Use `summary()` for clinical details and Cox regression.




# Core dependencies
install.packages(c("remotes", "survival"))

# Optional but highly recommended for multi-cut genetic optimisation
install.packages("rgenoud")

# Install the package from GitHub
remotes::install_github("paytonyau/OptSurvCutR")


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

library(OptSurvCutR)





required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Outcome: SLE instead of CAI 

# Annotate indivs SLE < CAI
df_binary_SLE_before_CAI <- df %>% 
  mutate(binary_SLE_before_CAI = case_when(
    is.na(SLE) & SLE_lower_age < CAI_lower_age ~ 1,
    TRUE ~ 0
  )) 



# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_SLE_before_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_SLE == 1, 
                                     SLE_lower_age - CAI_lower_age,
                                     age_end_follow_up - CAI_lower_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    SLE_lower_age_years = SLE_lower_age / 365.25,
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_SLE == 1,
                                      SLE_lower_age_years - CAI_lower_age_years,
                                      age_end_follow_up_years - CAI_lower_age_years
    ))




# Subset if needed
df_subset_cutpoints <- df_subset_years %>% 
  filter(follow_up_duration > 0) %>% 
  filter(binary_SLE_before_CAI != 1) %>% 
  filter(binary_ITP_or_AIHA_indivs == 1) %>% 
  select(PGS004917, follow_up_duration, binary_SLE)
nrow(df_subset_cutpoints)



# genetic method
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "PGS004917",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_SLE",
  num_cuts = 1,
  method = "genetic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result




# ____3.4 CAI PRS stratified ---------------------------------------------------------------------
# -UKB: Calculer OR et 95 CI association Score PRS stratifié avec CAI en multivarié

### CAI

# Logistic regression
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )






# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify score
df_stratified_score <- df %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. <= mean_SLE_PGS + 1 * sd_SLE_PGS ~ "low_risk",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. > mean_SLE_PGS + 1 * sd_SLE_PGS ~ "high_risk"
  ))




# Subset if needed
df_subset <- df_stratified_score
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  stratified_score +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  stratified_score +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.codmponents...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)






# ____3.5 subset CAI Cox SLE ~ PGS 1.662 ---------------------------------------------------------------------
# -UKB: Calculer OR et 95 CI association Score PRS stratifié avec CAI en multivarié

### CAI


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )



# Annotate 
df_binary_SLE_before_CAI <- df |> 
  mutate(binary_SLE_before_CAI = case_when(
    binary_ITP_or_AIHA_indivs == 1 & SLE_lower_age < CAI_lower_age ~ 1,
    TRUE ~ 0 
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_SLE_before_CAI %>%
  mutate(follow_up_duration = ifelse(binary_SLE == 1,
                                     SLE_lower_age - CAI_lower_age,
                                     age_end_follow_up - CAI_lower_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    SLE_lower_age_years = SLE_lower_age / 365.25,
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_SLE == 1,
                                      SLE_lower_age_years - CAI_lower_age_years,
                                      age_end_follow_up_years - CAI_lower_age_years
    ))






# Stratify score
df_stratified_score <- df_subset_years %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. <= 1.662 ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. > 1.662 ~"1"
  ))




# Subset if needed
df_subset <- df_stratified_score %>%
  filter(follow_up_duration > 0) %>%
  filter(binary_SLE_before_CAI == 0) %>% # Eliminate CAI before SLE
  filter(binary_ITP_or_AIHA_indivs == 1)
nrow(df_subset)





### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_SLE) ~ stratified_score +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)



# Non-adjusted
cox_model_unadjusted <- coxph(
  Surv(follow_up_duration, binary_SLE) ~ stratified_score,
  data = df_subset
)

summary(cox_model_unadjusted)





### Export for Prism
output <- df_subset %>% 
  mutate(binary_SLE_PRS_stratified_1p662 = stratified_score) %>% 
  select(binary_SLE_PRS_stratified_1p662, binary_SLE, follow_up_duration_years)



# Fig 2F

# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/articles/first_year_phd_pubs/002_first_article/para_julie/Fig2F_survival_subset_CAI_SLE_PGS_1p662.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')



















# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify score
df_stratified_score <- df %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. <= mean_SLE_PGS + 1 * sd_SLE_PGS ~ "low_risk",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. > mean_SLE_PGS + 1 * sd_SLE_PGS ~ "high_risk"
  ))




# Subset if needed
df_subset <- df_stratified_score
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  stratified_score +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  stratified_score +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.codmponents...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








# ___4) Effect LOF PGS CH and LM ------------------------------------


# ____4.1 Interaction of terms --------------------------------------------



# Interaction of covariables 

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))

table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))




# Subset if needed
df_subset <- df_follow_up_duration %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)



### Cox model interactions

# LOF : PRS SLE
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    WES_500k_LoF_MAF01__Hauck_ITPnHAI:Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)





# LOF : LM
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    WES_500k_LoF_MAF01__Hauck_ITPnHAI:binary_lymphoid_malignancies_after_CAI_coded_as_zero +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)




# LOF : CH
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    WES_500k_LoF_MAF01__Hauck_ITPnHAI:binary_CH_Age_sample +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)





# SLE : LM
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    Standard.PRS.for.systemic.lupus.erythematosus..SLE.:binary_lymphoid_malignancies_after_CAI_coded_as_zero +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)







# SLE : CH
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    Standard.PRS.for.systemic.lupus.erythematosus..SLE.:binary_CH_Age_sample +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)







# LM : CH
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    binary_lymphoid_malignancies_after_CAI_coded_as_zero:binary_CH_Age_sample +
    WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)



# _____4.1b Global model: 4 variant types, no interaction --------------------------------------------
# Supplemental table 5

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)





# Exact pval
p_values <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
print(p_values, digits = 20)





# Export for Prism

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, 
         follow_up_duration, 
         WES_500k_LoF_MAF01__Hauck_ITPnHAI,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         binary_lymphoid_malignancies_after_CAI_coded_as_zero,
         binary_CH_Age_sample)


print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/articles/first_year_phd_pubs/002_first_article/para_rafael/003_fig3_UKB_4groups.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')








# ______AoU replication --------------------------------------------
# Supplemental table 5

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))

nrow(df_follow_up_duration)

# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)

# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)





cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    AoU_LoF_MAF01__Hauck_ITPnHAI +
    PGS004917 +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    first_condition_age + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)





# Exact pval
p_values <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
print(p_values, digits = 20)









# ______AoU rep: first condition age --------------------------------------------

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     (age_end_follow_up*365.25) - first_condition_age))

nrow(df_follow_up_duration)

# CAI avant first condition: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)

# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)





cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    AoU_LoF_MAF01__Hauck_ITPnHAI +
    PGS004917 +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)





# Exact pval
p_values <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
print(p_values, digits = 20)














# _____4.1c Global model: CAI_lower_age 4 variant types, no interaction --------------------------------------------
# Supplemental table 10


# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))







# Subset if needed
df_subset <- df_lymphoid_malignancies_after_CAI %>% 
  filter(binary_ITP_or_AIHA_indivs == 1) 


nrow(df_subset)




summary(glm(
  as.numeric(CAI_lower_age / 365.25) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)
  
)




model_logistic_regression <- glm(
  as.numeric(CAI_lower_age / 365.25) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








# _____4.1d Global model: No HC, no interaction --------------------------------------------
# Supplemental table 6


# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))

table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))




# CAI avant first condition: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(age_inclusion) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)

# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)





# Exact pval
p_values <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
print(p_values, digits = 20)




# Export for Prism

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, 
         follow_up_duration, 
         WES_500k_LoF_MAF01__Hauck_ITPnHAI,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         binary_lymphoid_malignancies_after_CAI_coded_as_zero)


print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/articles/first_year_phd_pubs/002_first_article/para_rafael/001_fig3_UKB_3groups.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')




# ________ glm ------------------------------------------------------------------


df_subset <- df_CAI_before_blood_draw 

nrow(df_subset)



summary(glm(
  binary_ITP_or_AIHA_indivs ~      WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    age_end_follow_up +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset, family = binomial)
  
)



model_logistic_regression <- glm(
  binary_ITP_or_AIHA_indivs ~      WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    age_end_follow_up +
    omop_gender_concept_id +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset, family = binomial)





results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)












# ________AoU replication --------------------------------------------
# Supplemental table 6



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    #    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    #    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 






# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))

nrow(df_follow_up_duration)


# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)

# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)





cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
    AoU_LoF_MAF01__Hauck_ITPnHAI +
    PGS004917 +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    first_condition_age + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)





# Exact pval
p_values <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
print(p_values, digits = 20)





# Export for Prism

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, 
         follow_up_duration, 
         AoU_LoF_MAF01__Hauck_ITPnHAI,
         PGS004917,
         binary_lymphoid_malignancies_after_CAI_coded_as_zero)


print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/data_for_rafael_figure/002_fig3_AoU_3groups.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')





# __________ glm ------------------------------------------------------------------


df_subset <- df_CAI_before_blood_draw 

nrow(df_subset)



summary(glm(
  binary_ITP_or_AIHA_indivs ~      AoU_LoF_MAF01__Hauck_ITPnHAI +
    PGS004917 +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    age_end_follow_up +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)



model_logistic_regression <- glm(
  binary_ITP_or_AIHA_indivs ~      AoU_LoF_MAF01__Hauck_ITPnHAI +
    PGS004917 +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    age_end_follow_up +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)





results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)












# ____4.2 Score LOF_SLE_LM_CH combined ----------------------------------------
# Runs the global score: 4 variants types script
# Creates score by summing the HR * colmun value for each of the 4 terms

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt") 

# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))

table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))




# Subset if needed
df_subset <- df_follow_up_duration %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)




# Creates score by summing the HR * colmun value for each of the 4 terms
df_score_LoF_SLE_LM_CH <- df_lymphoid_malignancies_after_CAI %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    CH_clean = replace_na(binary_CH_Age_sample, 0),
    
    # Create score
    score_LoF_SLE_LM_CH = (1.6721 * LOF_clean) + 
      (1.1573 * PRS_clean) +
      (13.5987 * LM_clean) +
      (2.1425 * CH_clean)
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean, -CH_clean)  # remove temporary columns



output <- df_score_LoF_SLE_LM_CH


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# _____Aou Replication -----------------------------------------------------------


# CAI
required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom", 
                       "ggplot2", 
                       "tidyplots",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 



# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    
    TRUE ~ "1"
    
  ))


table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - biosample_collection_age,
                                     age_end_follow_up - biosample_collection_age))



# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(biosample_collection_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
    PGS004917 +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    binary_CH_Age_sample +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)




# Creates score by summing the HR * colmun value for each of the 4 terms
df_score_LoF_SLE_LM_CH <- df_lymphoid_malignancies_after_CAI %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    CH_clean = replace_na(binary_CH_Age_sample, 0),
    
    # Create score
    score_LoF_SLE_LM_CH = (1.6721 * LOF_clean) + 
      (1.1573 * PRS_clean) +
      (13.5987 * LM_clean) +
      (2.1425 * CH_clean)
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean, -CH_clean)  # remove temporary columns



output <- df_score_LoF_SLE_LM_CH


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260622_AoU_master_table_LoF_SLE_LM_CH.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')











# _____4.2a UKB score - no HC -----------------------------------------------------------



# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt") 

# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))

table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))




# Subset if needed
df_subset <- df_follow_up_duration %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    Age.at.recruitment + omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)




# Creates score by summing the HR * colmun value for each of the 4 terms
df_score_LoF_SLE_LM <- df_lymphoid_malignancies_after_CAI %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    score_LoF_SLE_LM = (1.8448 * LOF_clean) + 
      (1.1549 * PRS_clean) +
      (14.2431 * LM_clean) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns



output <- df_score_LoF_SLE_LM


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260623_summary_age_delays_score_LOF_SLE_LM.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')










# _____4.2b AoU score - no HC -----------------------------------------------------------



# CAI
required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom", 
                       "ggplot2", 
                       "tidyplots",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 



# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    
    TRUE ~ "1"
    
  ))


table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - biosample_collection_age,
                                     age_end_follow_up - biosample_collection_age))






# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(biosample_collection_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
    PGS004917 +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)




# Creates score by summing the HR * colmun value for each of the 4 terms
df_score_LoF_SLE_LM <- df_lymphoid_malignancies_after_CAI %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    
    # Create score
    score_LoF_SLE_LM = (1.8448 * LOF_clean) + 
      (1.1549 * PRS_clean) +
      (14.2431 * LM_clean) 
    
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns



output <- df_score_LoF_SLE_LM


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260624_AoU_master_table_LoF_SLE_LM.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')








# _____4.2c AoU score - no HC - first_condition_age -----------------------------------------------------------



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    
    TRUE ~ "1"
    
  ))


table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))






# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
    PGS004917 +
    binary_lymphoid_malignancies_after_CAI_coded_as_zero + 
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)




# Creates score by summing the HR * colmun value for each of the 4 terms
df_score_LoF_SLE_LM <- df_lymphoid_malignancies_after_CAI %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    
    # Create score
    score_LoF_SLE_LM = (1.8448 * LOF_clean) + 
      (1.1549 * PRS_clean) +
      (14.2431 * LM_clean) 
    
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns



output <- df_score_LoF_SLE_LM


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260624_AoU_master_table_LoF_SLE_LM.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')










# ____4.3 Cumulative incidence with score LOF_SLE_LM_CH -----------------------------------------------------------


# _____4.3a Cutpoints OptSurvCutR -------------------------------------------------


cat << \EOF > $UKBSCRIPTS/cutpoints_genetic_method.sh
#!/bin/bash
#SBATCH --job-name=2_3_cp_gen
#SBATCH --output=2_3_cp_gen%j.out
#SBATCH --error=2_3_cp_gen%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G

# Carregar módulos
module load r/4.4.0

# Configurações de performance
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_MAX_NUM_DLLS=1000
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Executar com logging detalhado
Rscript --vanilla --verbose -e '
library(dplyr)
library(OptSurvCutR)
library(survival)

# Configurar opções do R
options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(42)

# Carregar dados
cat("Carregando dados...\n")


# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/lustre09/project/6097258/stnfaria/ukb_rare_variants/000_whole_cohort/data/027_first_article_OptSurvCutR/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )


# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))



# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))



df_subset_cutpoints <- df_subset_years %>% 
  select(score_LoF_SLE_LM_CH, follow_up_duration, binary_ITP_or_AIHA_indivs) %>% 
  filter(follow_up_duration > 0) 


# Opção 2: Limpar memória antes de rodar
gc()


# 2 cps

# Opção 3: Rodar com método sistemático em vez de genético
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "genetic",
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


# 3 cps 

# Opção 3: Rodar com método sistemático em vez de genético
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 3,
  method = "genetic",
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result

'


EOF



sbatch $UKBSCRIPTS/cutpoints_genetic_method.sh



# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: score_LoF_SLE_LM_CH
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 2471.7215
# ✔ Recommended Cut-point(s): 11.826 and 13.112




# _____4.3b Analysis ------------------------------------------------------------



# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )


# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))



# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))



# 3 categories
df_stratified_score_genetic_method <- df_subset_years %>% 
  mutate(stratified_score_genetic_method = case_when(
    score_LoF_SLE_LM_CH <= 11.826 ~ "low_risk",
    score_LoF_SLE_LM_CH > 11.826 & score_LoF_SLE_LM_CH <= 13.112 ~ "medium_risk",
    score_LoF_SLE_LM_CH > 13.112 ~ "high_risk"
  ))




# Subset if needed
df_subset <- df_stratified_score_genetic_method %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)



freqs <- table(df_subset$stratified_score_genetic_method) %>% as.data.frame()


### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)





### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
  data=df_subset)


summary(cox_model)




# Very positive
### Puis tracer incidence cumulative AIC en fonction de LOF temps 0= naissance. evenement: diag AIC
# First, create binary variable very positive or not
df_subset$very_positive_binary <- ifelse(df_subset$stratified_score_genetic_method == "very_positive" , 1, 0)



# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
                  data = df_subset)

fit_km



### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). Combined LOF SLE LM CH score - genetic method cutpoints",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

# Add legend
legend("bottomright",
       legend = c("Non-carriers", "Carriers"),
       col = c("blue", "red"),
       lwd = 2)
###





### Save the plot as PNG
setwd("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/009_first_article/001_cox_curves")
png("005_combined_genetic_method_LOF_SLE_LM_CH_score_survival.png", width = 8, height = 6, units = "in", res = 300)

# Re-run the plotting code to save



### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). Combined LOF SLE LM CH score - genetic method cutpoints",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

# Add legend
legend("bottomright",
       legend = c("Non-carriers", "Carriers"),
       col = c("blue", "red"),
       lwd = 2)
###


dev.off()

# Confirm save location
message("Plot saved as: cumulative_incidence_curve_km.png in ", getwd())








# Export data to creawte graph on graphpad prism

output <- df_subset %>% 
  select(follow_up_duration_years, binary_ITP_or_AIHA_indivs,  stratified_score_genetic_method)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/003_stratified_score_genetic_method_for_prism.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')









# ______AoU replication------------------------------------------------------------




# CAI
required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom", 
                       "ggplot2", 
                       "tidyplots",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260622_AoU_master_table_LoF_SLE_LM_CH.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    score_LoF_SLE_LM_CH,
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))



# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - biosample_collection_age,
                                     age_end_follow_up - biosample_collection_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = biosample_collection_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))





# 3 categories
df_stratified_score_genetic_method <- df_subset_years %>% 
  mutate(stratified_score_genetic_method = case_when(
    score_LoF_SLE_LM_CH <= 11.826 ~ "low_risk",
    score_LoF_SLE_LM_CH > 11.826 & score_LoF_SLE_LM_CH <= 13.112 ~ "medium_risk",
    score_LoF_SLE_LM_CH > 13.112 ~ "high_risk"
  ))


nrow(df_stratified_score_genetic_method)





# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_stratified_score_genetic_method %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(biosample_collection_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




table(df_subset$stratified_score_genetic_method)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method +
    Age.at.recruitment + 
    sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)





### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
  data=df_subset)


summary(cox_model)




# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
                  data = df_subset)

fit_km








# Definir margens (um pouco mais espaçoso para a legenda)
par(mar = c(4, 4, 3, 2) + 0.1)

# Definir paleta de cores mais atraente
cores <- c("high_risk" = "#D73027",    # Vermelho escuro
           "medium_risk" = "#FDBB84",  # Laranja claro
           "low_risk" = "#4575B4")     # Azul escuro

# Plotar o gráfico
plot(fit_km, 
     fun = "event",
     col = c("#D73027", "#4575B4", "#FDBB84"),  # Ordem: high_risk, low_risk, medium_risk
     lwd = 2.5,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). Combined LOF SLE LM CH score - genetic method cutpoints",
     xaxt = "n")

# Adicionar eixo X
axis(1, at = seq(0, max(df_subset$follow_up_duration_years, na.rm = TRUE), by = 5))

# Adicionar legenda com as cores corretas
legend("topleft",
       legend = c("High risk", "Medium risk", "Low risk"),
       col = c("#D73027", "#FDBB84", "#4575B4"),
       lwd = 2.5,
       lty = 1,
       bty = "o",
       cex = 0.9,
       title = "Risk category")






### Save the plot as PNG
setwd("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/009_first_article/001_cox_curves")
png("005_combined_genetic_method_LOF_SLE_LM_CH_score_survival.png", width = 8, height = 6, units = "in", res = 300)

# Re-run the plotting code to save



### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). Combined LOF SLE LM CH score - genetic method cutpoints",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

# Add legend
legend("bottomright",
       legend = c("Non-carriers", "Carriers"),
       col = c("blue", "red"),
       lwd = 2)
###


dev.off()

# Confirm save location
message("Plot saved as: cumulative_incidence_curve_km.png in ", getwd())








# Export data to creawte graph on graphpad prism

output <- df_subset %>% 
  select(follow_up_duration_years, binary_ITP_or_AIHA_indivs,  stratified_score_genetic_method)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/003_stratified_score_genetic_method_for_prism.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')










# _____ (legacy)4.3c UKB stratified score - no HC ------------------------------------------------------------


# ______4.3d Cutpoints OptSurvCutR -------------------------------------------------



# Configure R options
options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(42)


# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)




df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260623_summary_age_delays_score_LOF_SLE_LM.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM
  )


# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))



# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))



df_subset_cutpoints <- df_subset_years %>% 
  select(score_LoF_SLE_LM, follow_up_duration, binary_ITP_or_AIHA_indivs) %>% 
  filter(follow_up_duration > 0) 


# Opção 2: Limpar memória antes de rodar
gc()


# 2 cps

# Opção 3: Rodar com método sistemático em vez de genético
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "genetic",
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result



# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: score_LoF_SLE_LM
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 2424.2298
# ✔ Recommended Cut-point(s): 12.465 and 13.758


# ______4.3e Analysis ------------------------------------------------------------

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260623_summary_age_delays_score_LOF_SLE_LM.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM
  )


# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))



# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))



# 3 categories
df_stratified_score_genetic_method <- df_subset_years %>% 
  mutate(stratified_score_genetic_method = case_when(
    score_LoF_SLE_LM <= 12.465 ~ "low_risk",
    score_LoF_SLE_LM > 12.465 & score_LoF_SLE_LM <= 13.758 ~ "medium_risk",
    score_LoF_SLE_LM > 13.758 ~ "high_risk"
  ))




# Subset if needed
df_subset <- df_stratified_score_genetic_method %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)



table(df_subset$stratified_score_genetic_method)


### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)





### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
  data=df_subset)


summary(cox_model)








# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
                  data = df_subset)

fit_km








# Definir margens (um pouco mais espaçoso para a legenda)
par(mar = c(4, 4, 3, 2) + 0.1)

# Definir paleta de cores mais atraente
cores <- c("high_risk" = "#D73027",    # Vermelho escuro
           "medium_risk" = "#FDBB84",  # Laranja claro
           "low_risk" = "#4575B4")     # Azul escuro

# Plotar o gráfico
plot(fit_km, 
     fun = "event",
     col = c("#D73027", "#4575B4", "#FDBB84"),  # Ordem: high_risk, low_risk, medium_risk
     lwd = 2.5,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). Combined LOF SLE LM score - genetic method cutpoints",
     xaxt = "n")

# Adicionar eixo X
axis(1, at = seq(0, max(df_subset$follow_up_duration_years, na.rm = TRUE), by = 5))

# Adicionar legenda com as cores corretas
legend("topleft",
       legend = c("High risk", "Medium risk", "Low risk"),
       col = c("#D73027", "#FDBB84", "#4575B4"),
       lwd = 2.5,
       lty = 1,
       bty = "o",
       cex = 0.9,
       title = "Risk category")











# Export data to creawte graph on graphpad prism

output <- df_subset %>% 
  select(follow_up_duration_years, binary_ITP_or_AIHA_indivs,  stratified_score_genetic_method)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/003_stratified_score_genetic_method_for_prism.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# ______AoU replication ------------------------------------------------------------




# CAI
required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom", 
                       "ggplot2", 
                       "tidyplots",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260624_AoU_master_table_LoF_SLE_LM.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    score_LoF_SLE_LM,
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))



# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - biosample_collection_age,
                                     age_end_follow_up - biosample_collection_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = biosample_collection_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))





# 3 categories
df_stratified_score_genetic_method <- df_subset_years %>% 
  mutate(stratified_score_genetic_method = case_when(
    score_LoF_SLE_LM <= 12.465 ~ "low_risk",
    score_LoF_SLE_LM > 12.465 & score_LoF_SLE_LM <= 13.758 ~ "medium_risk",
    score_LoF_SLE_LM > 13.758 ~ "high_risk"
  ))


nrow(df_stratified_score_genetic_method)





# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_stratified_score_genetic_method %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(biosample_collection_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)

# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)



table(df_subset$stratified_score_genetic_method)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method +
    Age.at.recruitment + 
    sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)





### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
  data=df_subset)


summary(cox_model)




# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
                  data = df_subset)

fit_km







# Use smaller margins
par(mar = c(4, 4, 2, 1) + 0.1)  # Reduced from c(4,4,3,2)

# Then plot
plot(fit_km, 
     fun = "event",
     col = c("#D73027", "#4575B4", "#FDBB84"),
     lwd = 2.5,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). Combined LOF SLE LM CH score - genetic method cutpoints",
     xaxt = "n")










### Save the plot as PNG
setwd("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/009_first_article/001_cox_curves")
png("005_combined_genetic_method_LOF_SLE_LM_CH_score_survival.png", width = 8, height = 6, units = "in", res = 300)

# Re-run the plotting code to save



### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). Combined LOF SLE LM CH score - genetic method cutpoints",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

# Add legend
legend("bottomright",
       legend = c("Non-carriers", "Carriers"),
       col = c("blue", "red"),
       lwd = 2)
###


dev.off()

# Confirm save location
message("Plot saved as: cumulative_incidence_curve_km.png in ", getwd())








# Export data to creawte graph on graphpad prism

output <- df_subset %>% 
  select(follow_up_duration_years, binary_ITP_or_AIHA_indivs,  stratified_score_genetic_method)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/003_stratified_score_genetic_method_for_prism.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')





# ______4.3f AoU rep: first condition age ------------------------------------------------------------




# CAI

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260624_AoU_master_table_LoF_SLE_LM.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    score_LoF_SLE_LM,
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 






# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))



# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))




# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = first_condition_age / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))





# 3 categories
df_stratified_score_genetic_method <- df_subset_years %>% 
  mutate(stratified_score_genetic_method = case_when(
    score_LoF_SLE_LM <= 12.465 ~ "low_risk",
    score_LoF_SLE_LM > 12.465 & score_LoF_SLE_LM <= 13.758 ~ "medium_risk",
    score_LoF_SLE_LM > 13.758 ~ "high_risk"
  ))


nrow(df_stratified_score_genetic_method)





# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_stratified_score_genetic_method %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)

# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)



table(df_subset$stratified_score_genetic_method)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method +
    Age.at.recruitment + 
    sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)





### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
  data=df_subset)


summary(cox_model)




# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ stratified_score_genetic_method,
                  data = df_subset)

fit_km







# Use smaller margins
par(mar = c(4, 4, 2, 1) + 0.1)  # Reduced from c(4,4,3,2)

# Then plot
plot(fit_km, 
     fun = "event",
     col = c("#D73027", "#4575B4", "#FDBB84"),
     lwd = 2.5,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). Combined LOF SLE LM CH score - genetic method cutpoints",
     xaxt = "n")










### Save the plot as PNG
setwd("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/009_first_article/001_cox_curves")
png("005_combined_genetic_method_LOF_SLE_LM_CH_score_survival.png", width = 8, height = 6, units = "in", res = 300)

# Re-run the plotting code to save



### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "Cumulative incidence curve (Kaplan-Meier). Combined LOF SLE LM CH score - genetic method cutpoints",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

# Add legend
legend("bottomright",
       legend = c("Non-carriers", "Carriers"),
       col = c("blue", "red"),
       lwd = 2)
###


dev.off()

# Confirm save location
message("Plot saved as: cumulative_incidence_curve_km.png in ", getwd())








# Export data to creawte graph on graphpad prism

output <- df_subset %>% 
  select(follow_up_duration_years, binary_ITP_or_AIHA_indivs,  stratified_score_genetic_method)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/003_stratified_score_genetic_method_for_prism.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')














# ____ 4.4 Second score - no HC -----------------------------------------------------------
# Supplemental table 7


#Fig 3B : 1.5 SD

# _______4.4a UKB -----------------------------------------------------------



# _________LM + LoF + PGS -----------------------------------------------------------

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt") 

# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))

table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))



# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))


# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_subset_years %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(Age.at.recruitment*365.25) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




# ___________1 SD score ----------------------------------------------------

### Score with 1 SD

# Creates score :
# -High risk if LM (regardless of LOF or PGS)
# -Médium if LoF or PGS >= 1 SD
# -Low for the other


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0 | Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= mean_SLE_PGS + 1 * sd_SLE_PGS ~ "medium_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns





### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model


# Multivariate

# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "high_risk", "medium_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]







# 
# ### Univariate
# 
# # Reorder so that high_risk is the reference (first level)
# df_subset$second_score <- factor(df_subset$second_score, 
#                                  levels = c("high_risk", "low_risk", "medium_risk"))
# 
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score,
#   data=df_subset)
# 
# 
# summary(cox_model)
# summary(cox_model)$coefficients[, "Pr(>|z|)"]
# 
# 
# 
# 
# # Reorder so that low_risk is the reference (first level)
# df_subset$second_score <- factor(df_subset$second_score, 
#                                  levels = c("low_risk", "high_risk", "medium_risk"))
# 
# model <- coxph(Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
#                  second_score, data = df_subset)
# 
# summary(model)
# summary(model)$coefficients[, "Pr(>|z|)"]





# ___________1.5 SD ----------------------------------------------------



# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0 | Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= mean_SLE_PGS + 1.5 * sd_SLE_PGS ~ "medium_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns





### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model




### Univariate - Unadjusted HR
# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "medium_risk", "high_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]



# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "high_risk", "low_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]










### Multivariate

# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "medium_risk", "high_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]







# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "high_risk", "low_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]







### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, second_score)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/articles/first_year_phd_pubs/002_first_article/para_julie/Fig3/001_fig3B_UKB_3riks_1p5SD.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# ___________2 SD score----------------------------------------------------


# Creates score :
# -High risk if LM (regardless of LOF or PGS)
# -Médium if LoF or PGS >= 2 SD
# -Low for the other


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0 | Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= mean_SLE_PGS + 2 * sd_SLE_PGS ~ "medium_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns







### Cumulative incidence 

# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)






### Cox model


# Multivariate

# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "high_risk", "medium_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]









### Univariate
# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk"))
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score,
  data=df_subset)


summary(cox_model)


summary(cox_model)$coefficients[, "Pr(>|z|)"]

### Check significance between low and medium
# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "high_risk",  "medium_risk"))

# Now run the model
model <- coxph(Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
                 second_score, data = df_subset)

summary(model)



summary(model)$coefficients[, "Pr(>|z|)"]







### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, second_score)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/to_prism/44a_UKB_3_risks_2SD.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')




# _________LM -----------------------------------------------------------

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt") 

# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))



# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))


# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_subset_years %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(Age.at.recruitment*365.25) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




### Cox model


# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_lymphoid_malignancies_after_CAI_coded_as_zero +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]






### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_lymphoid_malignancies_after_CAI_coded_as_zero,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]





### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, binary_lymphoid_malignancies_after_CAI_coded_as_zero)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/to_prism/44a_UKB_LM_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# _________LoF -----------------------------------------------------------

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt") 

# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))



# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))


# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_subset_years %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(Age.at.recruitment*365.25) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




### Cox model


# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]






### Univariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]






### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, WES_500k_LoF_MAF01__Hauck_ITPnHAI)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/to_prism/44a_UKB_LoF_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')










# _________SLE PRS binarised -----------------------------------------------------------

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt") 

# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))




# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))



# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))


# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_subset_years %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(Age.at.recruitment*365.25) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




# ___________1 SD score ----------------------------------------------------

### Score with 1 SD

# Creates score :
# -High risk if LM (regardless of LOF or PGS)
# -Médium if LoF or PGS >= 1 SD
# -Low for the other


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    binary_SLE = case_when(
      Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= mean_SLE_PGS + 1 * sd_SLE_PGS ~ "1",
      TRUE ~ "0"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns





### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model


# Multivariate

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_SLE +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]





### Univariate


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_SLE,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]





### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, binary_SLE)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/to_prism/44a_UKB_SLE_1SD_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')




# ___________2 SD score----------------------------------------------------


# Creates score :
# -High risk if LM (regardless of LOF or PGS)
# -Médium if LoF or PGS >= 2 SD
# -Low for the other


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    binary_SLE = case_when(
      Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= mean_SLE_PGS + 2 * sd_SLE_PGS ~ "1",
      TRUE ~ "0"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns







### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model


# Multivariate

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_SLE +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]





### Univariate


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_SLE,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]






### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, binary_SLE)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/to_prism/44a_UKB_SLE_2SD_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')











# _______4.4b AoU -----------------------------------------------------------

# _________LM + LoF + PGS -----------------------------------------------------------
required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    
    TRUE ~ "1"
    
  ))


table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))



# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




# ___________1 SD score ----------------------------------------------------

### Score with 1 SD

# Creates score :
# -High risk if LM (regardless of LOF or PGS)
# -Médium if LoF or PGS >= 1 SD
# -Low for the other


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      AoU_LoF_MAF01__Hauck_ITPnHAI > 0 | PGS004917 >= mean_SLE_PGS + 1 * sd_SLE_PGS ~ "medium_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns





### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




#### score frequency distribution
# 
# library(ggplot2)
# 
# ggplot(df_second_score, aes(x = PGS000196)) +
#   geom_bar(fill = "steelblue", color = "purple", alpha = 0.7) +
#   labs(title = "Frequência de cada valor de PGS000196",
#        x = "PGS000196",
#        y = "Frequência") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid = element_blank())  # Remove todas as grades










### Cox model


# Multivariate

# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "high_risk", "medium_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




### Univariate

# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "high_risk", "medium_risk"))

model <- coxph(Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
                 second_score, data = df_subset)

summary(model)
summary(model)$coefficients[, "Pr(>|z|)"]






### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, second_score)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/001_44_second_score_no_HC/44a_AoU_3_risks_1SD.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')







# ___________1.5 SD  ----------------------------------------------------



# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      AoU_LoF_MAF01__Hauck_ITPnHAI > 0 | PGS004917 >= mean_SLE_PGS + 1.5 * sd_SLE_PGS ~ "medium_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns





### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model




### Univariate - Unadjusted HR
# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "medium_risk", "high_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]



# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "high_risk", "low_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]







# Multivariate

# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "medium_risk", "high_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]





# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "high_risk", "low_risk"))




cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]









### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, second_score)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/004_fig_3/002_fig3B_AoU_3riks_1p5SD.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')



r






# ___________2 SD score----------------------------------------------------


# Creates score :
# -High risk if LM (regardless of LOF or PGS)
# -Médium if LoF or PGS >= 2 SD
# -Low for the other



# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      AoU_LoF_MAF01__Hauck_ITPnHAI > 0 | PGS004917 >= mean_SLE_PGS + 2 * sd_SLE_PGS ~ "medium_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns





### Cumulative incidence 

# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model


# Multivariate

# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "high_risk", "medium_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]









### Univariate

# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "high_risk", "medium_risk"))

model <- coxph(Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ 
                 second_score, data = df_subset)

summary(model)
summary(model)$coefficients[, "Pr(>|z|)"]








### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, second_score)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/001_44_second_score_no_HC/44a_AoU_3_risks_2SD.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')







# _________LM -----------------------------------------------------------

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    
    TRUE ~ "1"
    
  ))


table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))






# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




### Cox model


# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_lymphoid_malignancies_after_CAI_coded_as_zero +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]





### Univariate

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_lymphoid_malignancies_after_CAI_coded_as_zero,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]







### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, binary_lymphoid_malignancies_after_CAI_coded_as_zero)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/001_44_second_score_no_HC/44a_AoU_LM_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# _________LoF -----------------------------------------------------------


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    
    TRUE ~ "1"
    
  ))


table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))






# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




### Cox model


# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]





### Univariate

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]







### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, AoU_LoF_MAF01__Hauck_ITPnHAI)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/001_44_second_score_no_HC/44a_AoU_LoF_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# _________SLE PRS binarised -----------------------------------------------------------


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    
    TRUE ~ "1"
    
  ))


table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))





# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




# ___________1 SD ----------------------------------------------------


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)


# Create score
df_SLE_binarised <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    binary_SLE = case_when(
      PGS004917 >= mean_SLE_PGS + 1 * sd_SLE_PGS ~ "1",
      TRUE ~ "0"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns





### Cox model

df_subset <- df_SLE_binarised


# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_SLE +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]





### Univariate

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_SLE,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, binary_SLE)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/001_44_second_score_no_HC/44a_AoU_SLE_1SD_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')




# ___________2 SD ----------------------------------------------------


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)


# Create score
df_SLE_binarised <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    
    # Create score
    binary_SLE = case_when(
      PGS004917 >= mean_SLE_PGS + 2 * sd_SLE_PGS ~ "1",
      TRUE ~ "0"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean)  # remove temporary columns






### Cox model

df_subset <- df_SLE_binarised



# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_SLE +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]





### Univariate

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ binary_SLE,
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]






### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, binary_SLE)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/001_44_second_score_no_HC/44a_AoU_SLE_2SD_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')






# ____4.5 Second score with HC -----------------------------------------------------------
# Supplemental table 9


# _______4.5a UKB -----------------------------------------------------------

# _________LM + LoF + PGS + CH -----------------------------------------------------------



# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt") 

# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    TRUE ~ "1"
  ))

table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - age_inclusion,
                                     age_end_follow_up - age_inclusion))



# Convert ages from days to years BEFORE creating variables
df_subset_years <- df_follow_up_duration %>%
  mutate(
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    age_inclusion_years = age_inclusion / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - age_inclusion_years,
                                      age_end_follow_up_years - age_inclusion_years
    ))


# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_subset_years %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(Age.at.recruitment*365.25) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




# ___________1 SD score ----------------------------------------------------

### Score with 1 SD

# Creates score :
# -High risk if LM (regardless of LOF or PGS)
# -Médium if LoF or PGS >= 1 SD
# -Low for the other


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    CH_clean = replace_na(binary_CH_Age_sample, 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0 | Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= mean_SLE_PGS + 1 * sd_SLE_PGS ~ "medium_risk",
      CH_clean > 0 ~ "CH_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean, -CH_clean)  # remove temporary columns






### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model


# Multivariate

# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk", "CH_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "low_risk", "high_risk", "CH_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that CH_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("CH_risk", "medium_risk", "low_risk", "high_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]






# ___________1.5 SD  ----------------------------------------------------



# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    CH_clean = replace_na(binary_CH_Age_sample, 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0 | Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= mean_SLE_PGS + 1.5 * sd_SLE_PGS ~ "medium_risk",
      CH_clean > 0 ~ "CH_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean, -CH_clean)  # remove temporary columns






### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model



### Univariate - Unadjusted HR
# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "medium_risk", "high_risk", "CH_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "high_risk", "low_risk", "CH_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]



# Reorder so that CH_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("CH_risk", "high_risk", "low_risk", "medium_risk" ))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]













# Multivariate

# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "medium_risk", "high_risk", "CH_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "high_risk", "low_risk", "CH_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that CH_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("CH_risk", "high_risk", "low_risk", "medium_risk" ))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]







### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, second_score)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/articles/first_year_phd_pubs/002_first_article/para_julie/Fig3/003_fig3D_UKB_4riks_1p5SD.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')










# ___________1.75 SD  ----------------------------------------------------



# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    CH_clean = replace_na(binary_CH_Age_sample, 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0 | Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= mean_SLE_PGS + 1.75 * sd_SLE_PGS ~ "medium_risk",
      CH_clean > 0 ~ "CH_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean, -CH_clean)  # remove temporary columns






### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model


# Multivariate

# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk", "CH_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "low_risk", "high_risk", "CH_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that CH_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("CH_risk", "medium_risk", "low_risk", "high_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]













# ___________2 SD score----------------------------------------------------


# Creates score :
# -High risk if LM (regardless of LOF or PGS)
# -Médium if LoF or PGS >= 2 SD
# -Low for the other


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)


# Create score
df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(Standard.PRS.for.systemic.lupus.erythematosus..SLE., 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    CH_clean = replace_na(binary_CH_Age_sample, 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0 | Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= mean_SLE_PGS + 2 * sd_SLE_PGS ~ "medium_risk",
      CH_clean > 0 ~ "CH_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean, -CH_clean)  # remove temporary columns







### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model


# Multivariate

# Reorder so that high_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("high_risk", "low_risk", "medium_risk", "CH_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "low_risk", "high_risk", "CH_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that CH_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("CH_risk", "medium_risk", "low_risk", "high_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + 
    omop_gender_concept_id + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]









### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, second_score)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/articles/first_year_phd_pubs/002_first_article/para_julie/45a_UKB_4_risks_2SD.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')





# _______4.5b (legacy) AoU -----------------------------------------------------------
# Legacy: not used in the results because CH risk pval is not significant

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate lymphoid malignancies lower age
df_lymphoid_malignancies_lower_age <- df %>% 
  mutate(lymphoid_malignancies_lower_age = pmin(LYMPHOME_lower_age, Other_LEUK_lower_age, LLC_lower_age, na.rm = T))


# Annotate binary lymphoid malignancies after CAI
df_lymphoid_malignancies_after_CAI <- df_lymphoid_malignancies_lower_age %>% 
  
  mutate(binary_lymphoid_malignancies_after_CAI_coded_as_zero = case_when(
    
    !is.na(CAI_lower_age) & lymphoid_malignancies_lower_age > CAI_lower_age ~ "0",
    
    is.na(lymphoid_malignancies_lower_age) ~ "0",
    
    TRUE ~ "1"
    
  ))


table(df_lymphoid_malignancies_after_CAI$binary_lymphoid_malignancies_after_CAI, 
      useNA = "ifany")


# Annotate follow-up durations into a single column
df_follow_up_duration <- df_lymphoid_malignancies_after_CAI %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - first_condition_age,
                                     age_end_follow_up - first_condition_age))






# CAI avant prise sang: enlever
df_CAI_before_blood_draw <- df_follow_up_duration %>% 
  mutate(CAI_before_blood_draw = case_when(
    !is.na(CAI_lower_age) & as.numeric(CAI_lower_age) < as.numeric(first_condition_age) ~ "1",
    is.na(CAI_lower_age) ~ "0",
    TRUE ~ "0"
  ))

table(df_CAI_before_blood_draw$CAI_before_blood_draw)


nrow(df_CAI_before_blood_draw)



# Subset if needed
df_subset <- df_CAI_before_blood_draw %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_blood_draw != 1)

nrow(df_subset)




# ___________1 SD score ----------------------------------------------------


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)



df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    CH_clean = replace_na(binary_CH_Age_sample, 0),
    
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      AoU_LoF_MAF01__Hauck_ITPnHAI > 0 | PGS004917 >= mean_SLE_PGS + 1 * sd_SLE_PGS ~ "medium_risk",
      CH_clean > 0 ~ "CH_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean, -CH_clean)  # remove temporary columns






### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model


# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "medium_risk", "high_risk", "CH_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "high_risk", "low_risk", "CH_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]










# ___________1.5 SD ----------------------------------------------------


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)



df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    CH_clean = replace_na(binary_CH_Age_sample, 0),
    
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      AoU_LoF_MAF01__Hauck_ITPnHAI > 0 | PGS004917 >= mean_SLE_PGS + 1.5 * sd_SLE_PGS ~ "medium_risk",
      CH_clean > 0 ~ "CH_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean, -CH_clean)  # remove temporary columns





### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)





### Cox model


# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "medium_risk", "high_risk", "CH_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "high_risk", "low_risk", "CH_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]










# ___________2 SD score ----------------------------------------------------


# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)



df_second_score <- df_subset %>% 
  
  mutate(
    
    # Replace NAs with 0 to avoid NAs in the results column
    LOF_clean = replace_na(AoU_LoF_MAF01__Hauck_ITPnHAI, 0),
    PRS_clean = replace_na(PGS004917, 0),
    LM_clean = replace_na(as.numeric(binary_lymphoid_malignancies_after_CAI_coded_as_zero), 0),
    CH_clean = replace_na(binary_CH_Age_sample, 0),
    
    # Create score
    second_score = case_when(
      LM_clean > 0 ~ "high_risk",
      AoU_LoF_MAF01__Hauck_ITPnHAI > 0 | PGS004917 >= mean_SLE_PGS + 2 * sd_SLE_PGS ~ "medium_risk",
      CH_clean > 0 ~ "CH_risk",
      TRUE ~ "low_risk"
    ) 
  ) %>%
  select(-LOF_clean, -PRS_clean, -LM_clean, -CH_clean)  # remove temporary columns





### Cumulative incidence 


# Subset if needed
df_subset <- df_second_score %>% 
  filter(follow_up_duration > 0) 
nrow(df_subset)




### Cox model


# Reorder so that low_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("low_risk", "medium_risk", "high_risk", "CH_risk"))


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]




# Reorder so that medium_risk is the reference (first level)
df_subset$second_score <- factor(df_subset$second_score, 
                                 levels = c("medium_risk", "high_risk", "low_risk", "CH_risk"))

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ second_score +
    Age.at.recruitment + sex_genetic + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset)

summary(cox_model)
summary(cox_model)$coefficients[, "Pr(>|z|)"]






### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration, second_score)


# Write output
print('writing final table')
write.table(output,
            file = "/home/rstudio/workspace/workspace-bucket/stennio/003_first_article_replication/002_data_for_julie_fig2_code_2p1p2p3/45a_AoU_4_risks_2SD.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')








# ___5) CAI in high risk indivs ------------------------------------
# Supplemental table 12

# ____5.1 CAI after Lymphoma --------------------



library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LYMPHOME_lower_age,
                                     age_end_follow_up - LYMPHOME_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LYMPHOME_lower_age_years = LYMPHOME_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LYMPHOME_lower_age_years,
                                      age_end_follow_up_years - LYMPHOME_lower_age_years)
  )



# Annotate CAI before LYMPHOME
df_CAI_before_LYMPHOME <- df_follow_up_duration_years %>%
  mutate(CAI_before_LYMPHOME_indivs = case_when(
    CAI_lower_age < LYMPHOME_lower_age ~ "1",
    TRUE ~ "0"
  )) 





# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_LYMPHOME %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LYMPHOME_indivs != 1) # Remove indivs with CAI before LYMPHOME
# 
# nrow(df_subset)
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + omop_gender_concept_id + 
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)




### Logistic regression


# Subset if needed
df_subset <- df_CAI_before_LYMPHOME %>% 
  filter(CAI_before_LYMPHOME_indivs != 1) |> # Remove indivs with CAI before LYMPHOME
  filter(binary_LYMPHOME == 1)
nrow(df_subset)




## Covars : LoF and SLE PRS
summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)












# Cutpoints: Genetic and systematic

cat << \EOF > $UKBSCRIPTS/2cp_lymph.sh
#!/bin/bash
#SBATCH --job-name=2cp_lymph
#SBATCH --output=2cp_lymph%j.out
#SBATCH --error=2cp_lymph%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G



# Load modules
module load r/4.4.0

# Performance settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_MAX_NUM_DLLS=1000
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run with detailed logging
Rscript --vanilla --verbose -e '
library(dplyr)
library(OptSurvCutR)
library(survival)


# Configure R options
options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(42)

# Load data
cat("Loading data...\n")

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/lustre09/project/6097258/stnfaria/ukb_rare_variants/000_whole_cohort/data/027_first_article_OptSurvCutR/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )





# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LYMPHOME_lower_age,
                                     age_end_follow_up - LYMPHOME_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LYMPHOME_lower_age_years = LYMPHOME_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LYMPHOME_lower_age_years,
                                      age_end_follow_up_years - LYMPHOME_lower_age_years)
  )



# Annotate CAI before LYMPHOME
df_CAI_before_LYMPHOME <- df_follow_up_duration_years %>%
  mutate(CAI_before_LYMPHOME_indivs = case_when(
    CAI_lower_age < LYMPHOME_lower_age ~ "1",
    TRUE ~ "0"
  )) 



# Subset if needed
df_subset_cutpoints <- df_CAI_before_LYMPHOME %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_LYMPHOME_indivs != 1) %>% 
  select(score_LoF_SLE_LM_CH, follow_up_duration, binary_ITP_or_AIHA_indivs)
nrow(df_subset_cutpoints)



# Genetic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "genetic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


# Systematic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "systematic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


'


EOF



sbatch $UKBSCRIPTS/2cp_lymph.sh


# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: score_LoF_SLE_LM_CH
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 37.3844
# ✔ Recommended Cut-point(s): 16.248 and 16.748





# _______AoU replication --------------------
# Not significant LoF and CH in AoU, probably because of lack of statistical power

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LYMPHOME_lower_age,
                                     age_end_follow_up - LYMPHOME_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LYMPHOME_lower_age_years = LYMPHOME_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up /365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LYMPHOME_lower_age_years,
                                      age_end_follow_up_years - LYMPHOME_lower_age_years)
  )



# Annotate CAI before LYMPHOME
df_CAI_before_LYMPHOME <- df_follow_up_duration_years %>%
  mutate(CAI_before_LYMPHOME_indivs = case_when(
    CAI_lower_age < LYMPHOME_lower_age ~ "1",
    TRUE ~ "0"
  )) 






# 
# # Subset if needed
# df_subset <- df_CAI_before_LYMPHOME %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LYMPHOME_indivs != 1) # Remove indivs with CAI before LYMPHOME
# 
# nrow(df_subset)
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + sex_genetic + 
#     PGS004917 +
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)




### Logistic regression
# Subset if needed
df_subset <- df_CAI_before_LYMPHOME %>% 
  filter(CAI_before_LYMPHOME_indivs != 1) |> # Remove indivs with CAI before LYMPHOME
  filter(binary_LYMPHOME == 1)

nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)









# ____5.2 CAI after Other_Leuk --------------------


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - Other_LEUK_lower_age,
                                     age_end_follow_up - Other_LEUK_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    Other_LEUK_lower_age_years = Other_LEUK_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - Other_LEUK_lower_age_years,
                                      age_end_follow_up_years - Other_LEUK_lower_age_years)
  )



# Annotate CAI before Other_LEUK
df_CAI_before_Other_LEUK <- df_follow_up_duration_years %>%
  mutate(CAI_before_other_leuk_indivs = case_when(
    CAI_lower_age < Other_LEUK_lower_age ~ "1",
    TRUE ~ "0"
  )) 

# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_Other_LEUK %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_other_leuk_indivs != 1) # Remove indivs with CAI before Other_LEUK
# 
# nrow(df_subset)
# 
# 
#
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + omop_gender_concept_id + 
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 



### Logistic regression
# Subset if needed
df_subset <- df_CAI_before_Other_LEUK %>% 
  filter(CAI_before_other_leuk_indivs != 1) |>  # Remove indivs with CAI before Other_LEUK
  filter(binary_LYMPHOME == 1)

nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)





# Cutpoints: Genetic and systematic

cat << \EOF > $UKBSCRIPTS/2cp_otleuk.sh
#!/bin/bash
#SBATCH --job-name=2cp_otleuk
#SBATCH --output=2cp_otleuk%j.out
#SBATCH --error=2cp_otleuk%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G



# Load modules
module load r/4.4.0

# Performance settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_MAX_NUM_DLLS=1000
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run with detailed logging
Rscript --vanilla --verbose -e '
library(dplyr)
library(OptSurvCutR)
library(survival)


# Configure R options
options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(42)

# Load data
cat("Loading data...\n")

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/lustre09/project/6097258/stnfaria/ukb_rare_variants/000_whole_cohort/data/027_first_article_OptSurvCutR/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )







# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - Other_LEUK_lower_age,
                                     age_end_follow_up - Other_LEUK_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    Other_LEUK_lower_age_years = Other_LEUK_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - Other_LEUK_lower_age_years,
                                      age_end_follow_up_years - Other_LEUK_lower_age_years)
  )



# Annotate CAI before Other_LEUK
df_CAI_before_Other_LEUK <- df_follow_up_duration_years %>%
  mutate(CAI_before_other_leuk_indivs = case_when(
    CAI_lower_age < Other_LEUK_lower_age ~ "1",
    TRUE ~ "0"
  )) 



# Subset if needed
df_subset_cutpoints <- df_CAI_before_Other_LEUK %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_other_leuk_indivs != 1) %>% 
  select(score_LoF_SLE_LM_CH, follow_up_duration, binary_ITP_or_AIHA_indivs)
nrow(df_subset_cutpoints)



# Genetic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "genetic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


# Systematic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "systematic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


'


EOF



sbatch $UKBSCRIPTS/2cp_otleuk.sh



# 
# ── Optimal Cut-point Analysis for Survival Data (Genetic) 

# • Predictor: score_LoF_SLE_LM_CH
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 30.2628
# ✔ Recommended Cut-point(s): 16.19 and 16.358










# _______AoU replication --------------------
# Not significant LoF and CH in AoU, probably because of lack of statistical power

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 



# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - Other_LEUK_lower_age,
                                     age_end_follow_up - Other_LEUK_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    Other_LEUK_lower_age_years = Other_LEUK_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - Other_LEUK_lower_age_years,
                                      age_end_follow_up_years - Other_LEUK_lower_age_years)
  )



# Annotate CAI before Other_LEUK
df_CAI_before_Other_LEUK <- df_follow_up_duration_years %>%
  mutate(CAI_before_other_leuk_indivs = case_when(
    CAI_lower_age < Other_LEUK_lower_age ~ "1",
    TRUE ~ "0"
  )) 





# 
# # Subset if needed
# df_subset <- df_CAI_before_Other_LEUK %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_other_leuk_indivs != 1) # Remove indivs with CAI before Other_LEUK
# 
# nrow(df_subset)
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + sex_genetic + 
#     PGS004917 +
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 








### Logistic regression
# Subset if needed
df_subset <- df_CAI_before_Other_LEUK %>% 
  filter(CAI_before_other_leuk_indivs != 1) |> # Remove indivs with CAI before Other_LEUK
  filter(binary_Other_LEUK == 1)

nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)











# ____5.3 CAI after LLC --------------------


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LLC_lower_age,
                                     age_end_follow_up - LLC_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LLC_lower_age_lower_age_years = LLC_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LLC_lower_age_lower_age_years,
                                      age_end_follow_up_years - LLC_lower_age_lower_age_years)
  )



# Annotate CAI before LLC
df_CAI_before_LLC <- df_follow_up_duration_years %>%
  mutate(CAI_before_LLC_indivs = case_when(
    CAI_lower_age < LLC_lower_age ~ "1",
    TRUE ~ "0"
  )) 


# 
# # Subset if needed
# df_subset <- df_CAI_before_LLC %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LLC_indivs != 1) # Remove indivs with CAI before LLC
# 
# nrow(df_subset)
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + omop_gender_concept_id + 
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)




### Logistic regression
# Subset if needed
df_subset <- df_CAI_before_LLC %>% 
  filter(CAI_before_LLC_indivs != 1) |> # Remove indivs with CAI before LLC
  filter(binary_LLC == 1)
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)






# Cutpoints: Genetic and systematic

cat << \EOF > $UKBSCRIPTS/2cp_llc.sh
#!/bin/bash
#SBATCH --job-name=2cp_llc
#SBATCH --output=2cp_llc%j.out
#SBATCH --error=2cp_llc%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G



# Load modules
module load r/4.4.0

# Performance settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_MAX_NUM_DLLS=1000
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run with detailed logging
Rscript --vanilla --verbose -e '
library(dplyr)
library(OptSurvCutR)
library(survival)


# Configure R options
options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(42)

# Load data
cat("Loading data...\n")

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/lustre09/project/6097258/stnfaria/ukb_rare_variants/000_whole_cohort/data/027_first_article_OptSurvCutR/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )




# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LLC_lower_age,
                                     age_end_follow_up - LLC_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LLC_lower_age_lower_age_years = LLC_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LLC_lower_age_lower_age_years,
                                      age_end_follow_up_years - LLC_lower_age_lower_age_years)
  )



# Annotate CAI before LLC
df_CAI_before_LLC <- df_follow_up_duration_years %>%
  mutate(CAI_before_LLC_indivs = case_when(
    CAI_lower_age < LLC_lower_age ~ "1",
    TRUE ~ "0"
  )) 



# Subset if needed
df_subset_cutpoints <- df_CAI_before_LLC %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_LLC_indivs != 1) %>% 
  select(score_LoF_SLE_LM_CH, follow_up_duration, binary_ITP_or_AIHA_indivs)
nrow(df_subset_cutpoints)



# Genetic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "genetic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


# Systematic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "systematic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


'


EOF



sbatch $UKBSCRIPTS/2cp_llc.sh




# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: score_LoF_SLE_LM_CH
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 14.3328
# ✔ Recommended Cut-point(s): 13.385 and 14.069








# _______AoU replication --------------------
# Not significant LoF and CH in AoU, probably because of lack of statistical power

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 



# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LLC_lower_age,
                                     age_end_follow_up - LLC_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LLC_lower_age_lower_age_years = LLC_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LLC_lower_age_lower_age_years,
                                      age_end_follow_up_years - LLC_lower_age_lower_age_years)
  )



# Annotate CAI before LLC
df_CAI_before_LLC <- df_follow_up_duration_years %>%
  mutate(CAI_before_LLC_indivs = case_when(
    CAI_lower_age < LLC_lower_age ~ "1",
    TRUE ~ "0"
  )) 

# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_LLC %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LLC_indivs != 1) # Remove indivs with CAI before LLC
# 
# nrow(df_subset)
# 
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + sex_genetic + 
#     PGS004917 +
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)





### Logistic regression
# Subset if needed
df_subset <- df_CAI_before_LLC %>% 
  filter(CAI_before_LLC_indivs != 1) |> # Remove indivs with CAI before LLC
  filter(binary_LLC == 1)

nrow(df_subset)



summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)













# ____5.4 CAI after CH  --------------------


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - CH_Age_sample_lower_age,
                                     age_end_follow_up - CH_Age_sample_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    CH_lower_age_lower_age_years = CH_Age_sample_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - CH_lower_age_lower_age_years)
  )



# Annotate CAI before CH
df_CAI_before_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_CH_indivs = case_when(
    CAI_lower_age < CH_Age_sample_lower_age ~ "1",
    TRUE ~ "0"
  )) 



# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + omop_gender_concept_id + 
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)





### Logistic regression
# Subset if needed
df_subset <- df_CAI_before_CH %>% 
  filter(CAI_before_CH_indivs != 1) |> # Remove indivs with CAI before CH
  filter(binary_CH_Age_sample == 1)

nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)






# Cutpoints: Genetic and systematic

cat << \EOF > $UKBSCRIPTS/2cp_ch.sh
#!/bin/bash
#SBATCH --job-name=2cp_ch
#SBATCH --output=2cp_ch%j.out
#SBATCH --error=2cp_ch%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G



# Load modules
module load r/4.4.0

# Performance settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_MAX_NUM_DLLS=1000
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run with detailed logging
Rscript --vanilla --verbose -e '
library(dplyr)
library(OptSurvCutR)
library(survival)


# Configure R options
options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(42)

# Load data
cat("Loading data...\n")

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/lustre09/project/6097258/stnfaria/ukb_rare_variants/000_whole_cohort/data/027_first_article_OptSurvCutR/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )





# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - CH_Age_sample_lower_age,
                                     age_end_follow_up - CH_Age_sample_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    CH_lower_age_lower_age_years = CH_Age_sample_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - CH_lower_age_lower_age_years)
  )



# Annotate CAI before CH
df_CAI_before_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_CH_indivs = case_when(
    CAI_lower_age < CH_Age_sample_lower_age ~ "1",
    TRUE ~ "0"
  )) 



# Subset if needed
df_subset_cutpoints <- df_CAI_before_CH %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_CH_indivs != 1) %>% 
  select(score_LoF_SLE_LM_CH, follow_up_duration, binary_ITP_or_AIHA_indivs)
nrow(df_subset_cutpoints)



# Genetic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "genetic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


# Systematic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "systematic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


'


EOF


sbatch $UKBSCRIPTS/2cp_ch.sh


# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: score_LoF_SLE_LM_CH
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 157.7528
# ✔ Recommended Cut-point(s): 14.173 and 15.055







# _______AoU replication --------------------
# Not significant LoF and CH in AoU, probably because of lack of statistical power

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - CH_Age_sample_lower_age,
                                     age_end_follow_up - CH_Age_sample_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    CH_lower_age_lower_age_years = CH_Age_sample_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - CH_lower_age_lower_age_years)
  )



# Annotate CAI before CH
df_CAI_before_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_CH_indivs = case_when(
    CAI_lower_age < CH_Age_sample_lower_age ~ "1",
    TRUE ~ "0"
  )) 


# 
# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + sex_genetic + 
#     PGS004917 +
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)







### Logistic regression
# Subset if needed
df_subset <- df_CAI_before_CH %>% 
  filter(CAI_before_CH_indivs != 1) |> # Remove indivs with CAI before CH
  filter(binary_CH_Age_sample == 1)
nrow(df_subset)



summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)













# ____5.5 LYMHPOME OR OTHER_LEUK OR LLC OR CH after CAI --------------------


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# Annotate LM or CH lower age
df_LM_CH_lower_age <- df %>% 
  mutate(LM_CH_lower_age = pmin(LYMPHOME_lower_age, 
                                Other_LEUK_lower_age, 
                                LLC_lower_age,
                                CH_Age_sample_lower_age,
                                na.rm = T))


# Annotate binary LM or CH
df_binary_LM_CH <- df_LM_CH_lower_age |> 
  mutate(binary_LM_CH = case_when(
    !is.na(LM_CH_lower_age) ~ "1",
    TRUE ~ "0"
  ))



# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_LM_CH %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LM_CH_lower_age,
                                     age_end_follow_up - LM_CH_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LM_CH_lower_age_lower_age_years = LM_CH_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LM_CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - LM_CH_lower_age_lower_age_years)
  )



# Annotate CAI before LM_CH
df_CAI_before_LM_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_LM_CH_indivs = case_when(
    CAI_lower_age < LM_CH_lower_age ~ "1",
    TRUE ~ "0"
  )) 


# 
# # Subset if needed
# df_subset <- df_CAI_before_LM_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ WES_500k_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + omop_gender_concept_id + 
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)






### Logistic regression
# Subset if needed
df_subset <- df_CAI_before_LM_CH %>% 
  filter(CAI_before_LM_CH_indivs != 1) |> # Remove indivs with CAI before LM_CH
  filter(binary_LM_CH == 1)
nrow(df_subset)


summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
    omop_gender_concept_id +
    age_end_follow_up  + 
    scale(Standard.PRS.for.systemic.lupus.erythematosus..SLE.) + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)






# Cutpoints: Genetic and systematic

cat << \EOF > $UKBSCRIPTS/2cp_4ml.sh
#!/bin/bash
#SBATCH --job-name=2cp_4ml
#SBATCH --output=2cp_4ml%j.out
#SBATCH --error=2cp_4ml%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G



# Load modules
module load r/4.4.0

# Performance settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export R_MAX_NUM_DLLS=1000
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run with detailed logging
Rscript --vanilla --verbose -e '
library(dplyr)
library(OptSurvCutR)
library(survival)


# Configure R options
options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(42)

# Load data
cat("Loading data...\n")

# CAI
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)
library(OptSurvCutR)




df <- fread("/lustre09/project/6097258/stnfaria/ukb_rare_variants/000_whole_cohort/data/027_first_article_OptSurvCutR/20260514_summary_age_delays_score_LOF_SLE_LM_CH.tsv")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         
         # Score
         score_LoF_SLE_LM_CH
  )




# Annotate LM or CH lower age
df_LM_CH_lower_age <- df %>% 
  mutate(LM_CH_lower_age = pmin(LYMPHOME_lower_age, 
                                Other_LEUK_lower_age, 
                                LLC_lower_age,
                                CH_Age_sample_lower_age,
                                na.rm = T))


# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df_LM_CH_lower_age %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LM_CH_lower_age,
                                     age_end_follow_up - LM_CH_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LM_CH_lower_age_lower_age_years = LM_CH_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LM_CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - LM_CH_lower_age_lower_age_years)
  )



# Annotate CAI before LM_CH
df_CAI_before_LM_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_LM_CH_indivs = case_when(
    CAI_lower_age < LM_CH_lower_age ~ "1",
    TRUE ~ "0"
  )) 



# Subset if needed
df_subset_cutpoints <- df_CAI_before_LM_CH %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_LM_CH_indivs != 1) %>% 
  select(score_LoF_SLE_LM_CH, follow_up_duration, binary_ITP_or_AIHA_indivs)
nrow(df_subset_cutpoints)



# Genetic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "genetic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


# Systematic

# Option 2: Clear memory before running
gc()

# Option 3: Run with systematic method instead of genetic
cuts_result <- find_cutpoint(
  data = df_subset_cutpoints,
  predictor = "score_LoF_SLE_LM_CH",
  outcome_time = "follow_up_duration",
  outcome_event = "binary_ITP_or_AIHA_indivs",
  num_cuts = 2,
  method = "systematic",   
  criterion = "logrank",
  nmin = 50,
  quiet = FALSE,
  seed = 42
)

cuts_result


'


EOF


sbatch $UKBSCRIPTS/2cp_4ml.sh



# ── Optimal Cut-point Analysis for Survival Data (Genetic) 
# • Predictor: score_LoF_SLE_LM_CH
# • Criterion: logrank
# • Optimal Log-Rank Statistic: 311.2146
# ✔ Recommended Cut-point(s): 11.827 and 16.194










# _____________Special task: columns LOF, PRS_SLE, LM, CH, ITP binary, AIHA binary et age at first AIC ------------------------------------------------------------
# 2025 05 25
# D’ailleurs c’est pour la semaine prochaine mais je te l’écris là car je serai en clinique:
#   
#   Peux-tu m’envoyer les colonnes LOF, PRS_SLE, LM, CH, ITP binary, AIHA binary et age at first AIC pour les patients avec CAI stp? Je vais essayer un truc sur prism


# Annotate LM or CH lower age
df_LM_lower_age <- df_LM_CH_lower_age %>% 
  mutate(LM_lower_age = pmin(LYMPHOME_lower_age, 
                             Other_LEUK_lower_age, 
                             LLC_lower_age,
                             na.rm = T))


df_subset <- df_LM_lower_age %>% 
  filter(!is.na(CAI_lower_age)) 
nrow(df_subset)




output <- df_subset %>% 
  select(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
         Standard.PRS.for.systemic.lupus.erythematosus..SLE., 
         LM_lower_age, 
         CH_Age_sample, 
         binary_ITP, 
         binary_AIHA, 
         CAI_lower_age,
  )




# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260525_UKB_LOF_PRS_SLE_LM_CH_ITP_binary_AIHA_binary_CAI_lower_age_CAI_subset.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')





# Write output: table for ROC in prism
# Annotate LM or CH lower age
df_LM_lower_age <- df_LM_CH_lower_age %>% 
  mutate(LM_lower_age = pmin(LYMPHOME_lower_age, 
                             Other_LEUK_lower_age, 
                             LLC_lower_age,
                             na.rm = T))


df_subset <- df_LM_lower_age
nrow(df_subset)



output <- df_subset %>% 
  select(WES_500k_LoF_MAF01__Hauck_ITPnHAI, 
         Standard.PRS.for.systemic.lupus.erythematosus..SLE., 
         LM_lower_age, 
         CH_Age_sample, 
         binary_ITP, 
         binary_AIHA, 
         CAI_lower_age,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.,
         age_end_follow_up,
         CAI_lower_age,
         binary_ITP_or_AIHA_indivs
  ) %>% 
  mutate(omop_gender_concept_id = ifelse(is.na(omop_gender_concept_id), "", omop_gender_concept_id))





print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/001_UKB_data_for_ROC.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')







# _______AoU replication --------------------
# Not significant LoF and CH in AoU, probably because of lack of statistical power

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 



# Annotate LM or CH lower age
df_LM_CH_lower_age <- df %>% 
  mutate(LM_CH_lower_age = pmin(LYMPHOME_lower_age, 
                                Other_LEUK_lower_age, 
                                LLC_lower_age,
                                CH_Age_sample_lower_age,
                                na.rm = T))


# Annotate binary LM or CH
df_binary_LM_CH <- df_LM_CH_lower_age |> 
  mutate(binary_LM_CH = case_when(
    !is.na(LM_CH_lower_age) ~ "1",
    TRUE ~ "0"
  ))





# Annotate LM or CH lower age
df_LM_CH_lower_age <- df_binary_LM_CH %>% 
  mutate(LM_CH_lower_age = pmin(LYMPHOME_lower_age, 
                                Other_LEUK_lower_age, 
                                LLC_lower_age,
                                CH_Age_sample_lower_age,
                                na.rm = T))


# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df_LM_CH_lower_age %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LM_CH_lower_age,
                                     age_end_follow_up - LM_CH_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LM_CH_lower_age_lower_age_years = LM_CH_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LM_CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - LM_CH_lower_age_lower_age_years)
  )



# Annotate CAI before LM_CH
df_CAI_before_LM_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_LM_CH_indivs = case_when(
    CAI_lower_age < LM_CH_lower_age ~ "1",
    TRUE ~ "0"
  )) 


# 
# # Subset if needed
# df_subset <- df_CAI_before_LM_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ AoU_LoF_MAF01__Hauck_ITPnHAI +
#     Age.at.recruitment + sex_genetic + 
#     PGS004917 +
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)








### Logistic regression
# Subset if needed
df_subset <- df_CAI_before_LM_CH %>% 
  filter(CAI_before_LM_CH_indivs != 1) |> # Remove indivs with CAI before LM_CH
  filter(binary_LM_CH == 1)

nrow(df_subset)



summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  AoU_LoF_MAF01__Hauck_ITPnHAI +
    sex_genetic +
    age_end_follow_up  + 
    PGS004917 + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)












# ________________________(legacy) 5.6 PRS SLE and LoF binary summed- LM and CH --------------------
# PRS SLE stratified using the cutpoint of 2): 3.03 and binarised
# LoF binarised. Cutpoint: 0.1
# Annotate sum PRS SLE binarised + LoF binarised
# Subset to remove CAI before LM or CH
# Run cumulative incidence




library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# Annotate LM or CH lower age
df_LM_CH_lower_age <- df %>% 
  mutate(LM_CH_lower_age = pmin(LYMPHOME_lower_age, 
                                Other_LEUK_lower_age, 
                                LLC_lower_age,
                                CH_Age_sample_lower_age,
                                na.rm = T))


# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df_LM_CH_lower_age %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LM_CH_lower_age,
                                     age_end_follow_up - LM_CH_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LM_CH_lower_age_lower_age_years = LM_CH_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LM_CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - LM_CH_lower_age_lower_age_years)
  )



# Annotate CAI before LM_CH
df_CAI_before_LM_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_LM_CH_indivs = case_when(
    CAI_lower_age < LM_CH_lower_age ~ "1",
    TRUE ~ "0"
  )) 




# Stratify PRS score
df_stratified_score <- df_CAI_before_LM_CH %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. <= 3.03 ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. > 3.03 ~"1"
  ))




# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))



# Subset if needed
df_subset <- df_nouvelle_colonne %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH

nrow(df_subset)



table(df_subset$prs_n_lof)

table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)

table(subset(df_subset, stratified_score == "1")$lof_binaire)





### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
    Age.at.recruitment + omop_gender_concept_id  +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)













# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ prs_n_lof,
                  data = df_subset)

fit_km



### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     lwd = 2,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "outcome: SLE. T0: CAI lower age. SLE PRS 2 cutpoints. Subset: CAI indivs without primary SLE ",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

###










# ________________________(legacy) 5.7 PRS SLE and LoF binary summed- LM only --------------------
# PRS SLE stratified using the cutpoint of 2): 3.03 and binarised
# LoF binarised. Cutpoints: 0.1
# Annotate sum PRS SLE binarised + LoF binarised
# Subset to remove CAI before LM or CH
# Run cumulative incidence




library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# Annotate LM or CH lower age
df_LM_lower_age <- df %>% 
  mutate(LM_lower_age = pmin(LYMPHOME_lower_age, 
                             Other_LEUK_lower_age, 
                             LLC_lower_age,
                             na.rm = T))


# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df_LM_lower_age %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LM_lower_age,
                                     age_end_follow_up - LM_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LM_lower_age_lower_age_years = LM_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LM_lower_age_lower_age_years,
                                      age_end_follow_up_years - LM_lower_age_lower_age_years)
  )



# Annotate CAI before LM_CH
df_CAI_before_LM <- df_follow_up_duration_years %>%
  mutate(CAI_before_LM_indivs = case_when(
    CAI_lower_age < LM_lower_age ~ "1",
    TRUE ~ "0"
  )) 




# Stratify PRS score
df_stratified_score <- df_CAI_before_LM %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. <= 3.03 ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. > 3.03 ~"1"
  ))




# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))



# Subset if needed
df_subset <- df_nouvelle_colonne %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_LM_indivs != 1) # Remove indivs with CAI before LM

nrow(df_subset)



### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
    Age.at.recruitment + omop_gender_concept_id  +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)













# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ prs_n_lof,
                  data = df_subset)

fit_km



### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     lwd = 2,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "outcome: SLE. T0: CAI lower age. SLE PRS 2 cutpoints. Subset: CAI indivs without primary SLE ",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

###





# Export for Prism

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, prs_n_lof, follow_up_duration)


print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/5p7_PRS_SLE_LOF_binary_summed_LM_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')







# ________________________(legacy) 5.8 PRS SLE and LoF binary summed- CH only --------------------
# PRS SLE stratified using the cutpoint of 2): 3.03 and binarised
# LoF binarised. Cutpoint: 0.1
# Annotate sum PRS SLE binarised + LoF binarised
# Subset to remove CAI before CH only
# Run cumulative incidence


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - CH_Age_sample_lower_age,
                                     age_end_follow_up - CH_Age_sample_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    CH_lower_age_lower_age_years = CH_Age_sample_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - CH_lower_age_lower_age_years)
  )



# Annotate CAI before LM_CH
df_CAI_before_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_CH_indivs = case_when(
    CAI_lower_age < CH_Age_sample_lower_age ~ "1",
    TRUE ~ "0"
  )) 




# Stratify PRS score
df_stratified_score <- df_CAI_before_CH %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. <= 3.03 ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. > 3.03 ~"1"
  ))




# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))




# Subset if needed
df_subset <- df_nouvelle_colonne %>% 
  filter(follow_up_duration > 0) %>% 
  filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH

nrow(df_subset)




### Cox model

# Multivariate
cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
    Age.at.recruitment + omop_gender_concept_id  +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset)


summary(cox_model)










# Kaplan-Meier curve (equivalent to Cox model with binary variable)
fit_km <- survfit(Surv(follow_up_duration_years, binary_ITP_or_AIHA_indivs) ~ prs_n_lof,
                  data = df_subset)

fit_km



### Graph formatting
# Plot without X axis initially
plot(fit_km, 
     fun = "event",                # Cumulative incidence (1 - survival)
     lwd = 2,
     xlab = "Time of follow-up (years)",
     ylab = "Cumulative incidence",
     main = "outcome: SLE. T0: CAI lower age. SLE PRS 2 cutpoints. Subset: CAI indivs without primary SLE ",
     xaxt = "n")                   # Suppresses automatic X axis

# Add manual X axis with 5-year intervals
max_age <- max(df_subset$follow_up_duration_years, na.rm = TRUE)
breaks_x <- seq(0, ceiling(max_age/5)*5, by = 5)

axis(side = 1, at = breaks_x, labels = breaks_x)

# Add vertical grid at the breaks
abline(v = breaks_x, col = "lightgray", lty = 2)

###





# Export for Prism

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, prs_n_lof, follow_up_duration)


print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/5p8_PRS_SLE_LOF_binary_summed_CH_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')







# ____ (legacy) 5.6 PRS_n_LoF binary summed- LM and CH --------------------
# Legacy : not used in final results 


# PRS SLE stratified using the cutpoint of 2): 3.03 and binarised
# LoF binarised. Cutpoint: 0.1
# Annotate sum PRS SLE binarised + LoF binarised
# Subset to remove CAI before LM or CH
# Run cumulative incidence



# Supplemental table 13

# Fig 5B: Section 5.7 1|2 SD SLE  PGS and LoF binary summed- LM only, 1.5 SD


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )


# Annotate LM or CH lower age
df_LM_CH_lower_age <- df %>% 
  mutate(LM_CH_lower_age = pmin(LYMPHOME_lower_age, 
                                Other_LEUK_lower_age, 
                                LLC_lower_age,
                                CH_Age_sample_lower_age,
                                na.rm = T))



# Annotate binary LM or CH
df_binary_LM_CH <- df_LM_CH_lower_age |> 
  mutate(binary_LM_CH = case_when(
    !is.na(LM_CH_lower_age) ~ "1",
    TRUE ~ "0"
  ))




# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_LM_CH %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LM_CH_lower_age,
                                     age_end_follow_up - LM_CH_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LM_CH_lower_age_lower_age_years = LM_CH_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LM_CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - LM_CH_lower_age_lower_age_years)
  )



# Annotate CAI before LM_CH
df_CAI_before_LM_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_LM_CH_indivs = case_when(
    CAI_lower_age < LM_CH_lower_age ~ "1",
    TRUE ~ "0"
  )) 



# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_LM_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
#
# 
# ______1 SD --------------------------------------------------------------------
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# table(df_subset$prs_n_lof)
# 
# table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)
# 
# table(subset(df_subset, stratified_score == "1")$lof_binaire)
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + omop_gender_concept_id  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 




### Logistic regression

# Subset if needed
df_subset <- df_CAI_before_LM_CH %>% 
  filter(binary_LM_CH == 1) |> 
  filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH

nrow(df_subset)



# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))



# Subset if needed
df_subset <- df_nouvelle_colonne



## Logistic regression

summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)






# ______1.5 SD --------------------------------------------------------------------
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# table(df_subset$prs_n_lof)
# 
# table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)
# 
# table(subset(df_subset, stratified_score == "1")$lof_binaire)
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + omop_gender_concept_id  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 




### Logistic regression

# Subset if needed
df_subset <- df_CAI_before_LM_CH %>% 
  filter(binary_LM_CH == 1) |> 
  filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH

nrow(df_subset)



# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1.5 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1.5 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))



# Subset if needed
df_subset <- df_nouvelle_colonne



## Logistic regression

summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)
























# ______2 SD --------------------------------------------------------------------

# 
# ### Cox model
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_LM_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# table(df_subset$prs_n_lof)
# 
# table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)
# 
# table(subset(df_subset, stratified_score == "1")$lof_binaire)
# 
# 
# 
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + omop_gender_concept_id  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 








### Logistic regression

# Subset if needed
df_subset <- df_CAI_before_LM_CH %>% 
  filter(binary_LM_CH == 1) |> 
  filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH

nrow(df_subset)



# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))



# Subset if needed
df_subset <- df_nouvelle_colonne



## Logistic regression

summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)










# __________ (legacy) AoU replication --------------------
# Legacy: not used because results are not significant in AoU v8

required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 





# Annotate LM or CH lower age
df_LM_CH_lower_age <- df %>% 
  mutate(LM_CH_lower_age = pmin(LYMPHOME_lower_age, 
                                Other_LEUK_lower_age, 
                                LLC_lower_age,
                                CH_Age_sample_lower_age,
                                na.rm = T))



# Annotate binary LM or CH
df_binary_LM_CH <- df_LM_CH_lower_age |> 
  mutate(binary_LM_CH = case_when(
    !is.na(LM_CH_lower_age) ~ "1",
    TRUE ~ "0"
  ))




# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_LM_CH %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LM_CH_lower_age,
                                     age_end_follow_up - LM_CH_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LM_CH_lower_age_lower_age_years = LM_CH_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LM_CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - LM_CH_lower_age_lower_age_years)
  )



# Annotate CAI before LM_CH
df_CAI_before_LM_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_LM_CH_indivs = case_when(
    CAI_lower_age < LM_CH_lower_age ~ "1",
    TRUE ~ "0"
  )) 







# _____________1 SD --------------------------------------------------------------------

# 
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_LM_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     PGS004917 < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
#     PGS004917 >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     AoU_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     AoU_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + sex_genetic + 
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 
# 


# _____________2 SD --------------------------------------------------------------------
# 
# 
# 
# ### Cox model
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     PGS004917 < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
#     PGS004917 >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     AoU_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     AoU_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + sex_genetic + 
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 
# 






### Logistic regression

# Subset if needed
df_subset <- df_CAI_before_LM_CH %>%
  filter(binary_LM_CH == 1) %>%
  filter(CAI_before_LM_CH_indivs != 1) # Remove indivs with CAI before LM_CH

nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)




# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    PGS004917 < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
    PGS004917 >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    AoU_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    AoU_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))




# Subset if needed
df_subset <- df_nouvelle_colonne




summary(glm(
  binary_ITP_or_AIHA_indivs ~  prs_n_lof +
    age_end_follow_up +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)



model_logistic_regression <- glm(
  binary_ITP_or_AIHA_indivs ~  prs_n_lof +
    age_end_follow_up +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)





results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)












# ____5.7 1|2 SD SLE  PGS and LoF binary summed- LM only --------------------
# PRS SLE stratified using mean + 2 SD and binarised
# LoF binarised. Cutpoint: 0.1
# Annotate sum PRS SLE binarised + LoF binarised
# Subset to remove CAI before LM or CH
# Run cumulative incidence


# Fig 5B: 1.5 SD

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )

df <- df |> 
  mutate(`Age.at.recruitment` = `Age.at.recruitment` * 365.25)


# Annotate LM lower age
df_LM_lower_age <- df %>% 
  mutate(LM_lower_age = pmin(LYMPHOME_lower_age, 
                             Other_LEUK_lower_age, 
                             LLC_lower_age,
                             na.rm = T))


# Annotate binary LM or CH
df_binary_LM <- df_LM_lower_age |> 
  mutate(binary_LM = case_when(
    !is.na(LM_lower_age) ~ "1",
    TRUE ~ "0"
  ))






# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_LM %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LM_lower_age,
                                     age_end_follow_up - LM_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LM_lower_age_lower_age_years = LM_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LM_lower_age_lower_age_years,
                                      age_end_follow_up_years - LM_lower_age_lower_age_years)
  )



# Annotate CAI before LM
df_CAI_before_LM <- df_follow_up_duration_years %>%
  mutate(CAI_before_LM_indivs = case_when(
    CAI_lower_age < LM_lower_age ~ "1",
    TRUE ~ "0"
  )) 






# ______1 SD --------------------------------------------------------------------




# Subset if needed
df_subset <- df_CAI_before_LM %>%
  filter(follow_up_duration > 0) %>%
  filter(binary_LM == 1) |> 
  filter(CAI_before_LM_indivs != 1) |> # Remove indivs with CAI before LM
  filter(!is.na(Standard.PRS.for.systemic.lupus.erythematosus..SLE.)) # Remove indivs without SLE PGS

  
nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)




# Stratify PRS score
df_stratified_score <- df_subset %>%
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>%
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>%
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score)) |> 
  mutate(prs_n_lof = case_when(
    prs_n_lof > 0 ~ "1",
    TRUE ~ "0"
  ))




# Subset if needed
df_subset <- df_nouvelle_colonne






table(df_subset$prs_n_lof)

table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)

table(subset(df_subset, stratified_score == "1")$lof_binaire)






### Cox model

# Multivariate
cox_model <- coxph(

  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
    Age.at.recruitment + omop_gender_concept_id  +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset)


summary(cox_model)




# 
# 
# ### Logistic regression
# 
# # Subset if needed
# df_subset <- df_CAI_before_LM %>%
#   filter(binary_LM == 1) |>
#   filter(CAI_before_LM_indivs != 1) # Remove indivs with CAI before LM_CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>%
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>%
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>%
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# table(df_subset$prs_n_lof)
# 
# table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)
# 
# table(subset(df_subset, stratified_score == "1")$lof_binaire)
# 
# 
# 
# 
# ## Logistic regression
# 
# summary(glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
#     omop_gender_concept_id +
#     age_end_follow_up  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
#   data=df_subset, family = binomial)
# 
# )
# 
# 
# 
# 
# model_logistic_regression <- glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
#     omop_gender_concept_id +
#     age_end_follow_up  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
#   data=df_subset, family = binomial)
# 
# 
# results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
#   mutate(
#     OR = exp(estimate),
#     CI_inf = exp(conf.low),
#     CI_sup = exp(conf.high)
#   ) %>%
#   select(term, OR, CI_inf, CI_sup, p.value, statistic)
# 
# print(results, digits = 3)
# 
# 
# 
# summary_model <- summary(model_logistic_regression)
# 
# # Extract coefficients with p-values
# coef_table <- summary_model$coefficients
# 
# # View the full table with p-values (exact)
# print(coef_table, digits = 20)
# 






# ______1.5 SD --------------------------------------------------------------------





# Subset if needed
df_subset <- df_CAI_before_LM %>%
  filter(follow_up_duration > 0) %>%
  filter(binary_LM == 1) |> 
  filter(CAI_before_LM_indivs != 1) |> # Remove indivs with CAI before LM
  filter(!is.na(Standard.PRS.for.systemic.lupus.erythematosus..SLE.)) # Remove indivs without SLE PGS


nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)




# Stratify PRS score
df_stratified_score <- df_subset %>%
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1.5 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1.5 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>%
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))




# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>%
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score)) |> 
  mutate(prs_n_lof = case_when(
    prs_n_lof > 0 ~ "1",
    TRUE ~ "0"
  ))



# Subset if needed
df_subset <- df_nouvelle_colonne



table(df_subset$prs_n_lof)

table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)

table(subset(df_subset, stratified_score == "1")$lof_binaire)







### Cox model



### Univariate - unadjusted HR

cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof,
  data=df_subset)


summary(cox_model)






### Multivariate - adjusted HR
cox_model <- coxph(

  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
    Age.at.recruitment + omop_gender_concept_id  +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset)


summary(cox_model)





# 
# ### Logistic regression
# 
# # Subset if needed
# df_subset <- df_CAI_before_LM %>%
#   filter(binary_LM == 1) |> 
#   filter(CAI_before_LM_indivs != 1) # Remove indivs with CAI before LM_CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1.5 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1.5 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# table(df_subset$prs_n_lof)
# 
# table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)
# 
# table(subset(df_subset, stratified_score == "1")$lof_binaire)
# 
# 
# 
# 
# ## Logistic regression
# 
# 
# 
# ### Univariate - Unadjusted HR
# summary(glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof, 
#   data=df_subset, family = binomial)
#   
# )
# 
# 
# 
# 
# model_logistic_regression <- glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof, 
#   data=df_subset, family = binomial)
# 
# 
# results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
#   mutate(
#     OR = exp(estimate),
#     CI_inf = exp(conf.low),
#     CI_sup = exp(conf.high)
#   ) %>%
#   select(term, OR, CI_inf, CI_sup, p.value, statistic)
# 
# print(results, digits = 3)
# 
# 
# 
# summary_model <- summary(model_logistic_regression)
# 
# # Extract coefficients with p-values
# coef_table <- summary_model$coefficients
# 
# # View the full table with p-values (exact)
# print(coef_table, digits = 20)
# 
# 
# 
# 
# 
# 
# 
# 
# ### Multivariate - Adjusted HR
# 
# summary(glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
#     omop_gender_concept_id +
#     age_end_follow_up  + 
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset, family = binomial)
#   
# )
# 
# 
# 
# 
# model_logistic_regression <- glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
#     omop_gender_concept_id +
#     age_end_follow_up  + 
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset, family = binomial)
# 
# 
# results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
#   mutate(
#     OR = exp(estimate),
#     CI_inf = exp(conf.low),
#     CI_sup = exp(conf.high)
#   ) %>%
#   select(term, OR, CI_inf, CI_sup, p.value, statistic)
# 
# print(results, digits = 3)
# 
# 
# 
# summary_model <- summary(model_logistic_regression)
# 
# # Extract coefficients with p-values
# coef_table <- summary_model$coefficients
# 
# # View the full table with p-values (exact)
# print(coef_table, digits = 20)
# 




### Exporting data for Prism:

output <- df_subset %>% 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration_years, prs_n_lof)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/articles/first_year_phd_pubs/002_first_article/para_julie/Fig3/004_fig5B_UKB_PRSandLoF_1p5SD.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')







#______2 SD --------------------------------------------------------------------



### Cox


# Subset if needed
df_subset <- df_CAI_before_LM %>%
  filter(follow_up_duration > 0) %>%
  filter(binary_LM == 1) |> 
  filter(CAI_before_LM_indivs != 1) |> # Remove indivs with CAI before LM
  filter(!is.na(Standard.PRS.for.systemic.lupus.erythematosus..SLE.)) # Remove indivs without SLE PGS


nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify PRS score
df_stratified_score <- df_subset %>%
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>%
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))




# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>%
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score)) |> 
  mutate(prs_n_lof = case_when(
    prs_n_lof > 0 ~ "1",
    TRUE ~ "0"
  ))



# Subset if needed
df_subset <- df_nouvelle_colonne



table(df_subset$prs_n_lof)

table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)

table(subset(df_subset, stratified_score == "1")$lof_binaire)





### Cox model

# Multivariate
cox_model <- coxph(

  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
    Age.at.recruitment + omop_gender_concept_id  +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset)


summary(cox_model)








# 
# 
# ### Logistic regression
# 
# # Subset if needed
# df_subset <- df_CAI_before_LM %>%
#   filter(binary_LM == 1) |> 
#   filter(CAI_before_LM_indivs != 1) # Remove indivs with CAI before LM_CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# table(df_subset$prs_n_lof)
# 
# table(df_subset$WES_500k_LoF_MAF01__Hauck_ITPnHAI)
# 
# table(subset(df_subset, stratified_score == "1")$lof_binaire)
# 
# 
# 
# 
# ## Logistic regression
# 
# summary(glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
#     omop_gender_concept_id +
#     age_end_follow_up  + 
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset, family = binomial)
#   
# )
# 
# 
# 
# 
# model_logistic_regression <- glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
#     omop_gender_concept_id +
#     age_end_follow_up  + 
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset, family = binomial)
# 
# 
# results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
#   mutate(
#     OR = exp(estimate),
#     CI_inf = exp(conf.low),
#     CI_sup = exp(conf.high)
#   ) %>%
#   select(term, OR, CI_inf, CI_sup, p.value, statistic)
# 
# print(results, digits = 3)
# 
# 
# 
# summary_model <- summary(model_logistic_regression)
# 
# # Extract coefficients with p-values
# coef_table <- summary_model$coefficients
# 
# # View the full table with p-values (exact)
# print(coef_table, digits = 20)
# 









# __________ (legacy) AoU replication --------------------


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




# Annotate LM lower age
df_LM_lower_age <- df %>% 
  mutate(LM_lower_age = pmin(LYMPHOME_lower_age, 
                             Other_LEUK_lower_age, 
                             LLC_lower_age,
                             na.rm = T))


# Annotate binary LM or CH
df_binary_LM <- df_LM_lower_age |> 
  mutate(binary_LM = case_when(
    !is.na(LM_lower_age) ~ "1",
    TRUE ~ "0"
  ))



# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df_binary_LM %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - LM_lower_age,
                                     age_end_follow_up - LM_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    LM_lower_age_lower_age_years = LM_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - LM_lower_age_lower_age_years,
                                      age_end_follow_up_years - LM_lower_age_lower_age_years)
  )



# Annotate CAI before LM
df_CAI_before_LM <- df_follow_up_duration_years %>%
  mutate(CAI_before_LM_indivs = case_when(
    CAI_lower_age < LM_lower_age ~ "1",
    TRUE ~ "0"
  )) 






# _____________1 SD --------------------------------------------------------------------

# 
# # Subset if needed
# df_subset <- df_CAI_before_LM %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_LM_indivs != 1) # Remove indivs with CAI before LM
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)
# 
# 
# 

# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     PGS004917 < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
#     PGS004917 >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     AoU_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     AoU_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + sex_genetic + 
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 





# _____________2 SD --------------------------------------------------------------------

# 
# 
# ### Cox model
# 
# # Subset if needed
# df_subset <- df_CAI_before_LM %>% 
#   filter(binary_LM == 1) %>% 
#   filter(follow_up_duration > 0) %>%
#   filter(CAI_before_LM_indivs != 1) # Remove indivs with CAI before LM
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)
# 
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     PGS004917 < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
#     PGS004917 >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     AoU_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     AoU_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# 
# 
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + sex_genetic + 
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)





### Logistic regression


# Subset if needed
df_subset <- df_CAI_before_LM %>% 
  filter(binary_LM == 1) %>% 
  filter(CAI_before_LM_indivs != 1) # Remove indivs with CAI before LM

nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)




# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    PGS004917 < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
    PGS004917 >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    AoU_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    AoU_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))




# Subset if needed
df_subset <- df_nouvelle_colonne




### Logistic regression

summary(glm(
  binary_ITP_or_AIHA_indivs ~  prs_n_lof +
    age_end_follow_up +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)



model_logistic_regression <- glm(
  binary_ITP_or_AIHA_indivs ~  prs_n_lof +
    age_end_follow_up +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)





results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)











# ____5.8 1|2 SD SLE  PGS and LoF binary summed- CH only --------------------
# PRS SLE stratified using mean + 2 SD and binarised
# LoF binarised. Cutpoint: 0.1
# Annotate sum PRS SLE binarised + LoF binarised
# Subset to remove CAI before LM or CH
# Run cumulative incidence




library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(survminer)
library(ggsurvfit)


df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         binary_LYMPHOME,
         binary_Other_LEUK,
         binary_LLC,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )




# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - CH_Age_sample_lower_age,
                                     age_end_follow_up - CH_Age_sample_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    CH_lower_age_lower_age_years = CH_Age_sample_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - CH_lower_age_lower_age_years)
  )



# Annotate CAI before CH
df_CAI_before_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_CH_indivs = case_when(
    CAI_lower_age < CH_Age_sample_lower_age ~ "1",
    TRUE ~ "0"
  )) 






# ______ LoF only --------------------------------------------------------------------





# Annotate LOF binary
df_lof_binaire <- df_CAI_before_CH %>%
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI == 0 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0 ~"1"
  ))


# Subset if needed
df_subset <- df_lof_binaire %>%
  filter(follow_up_duration > 0) %>%
  filter(binary_CH_Age_sample == 1) |> 
  filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH

nrow(df_subset)





### Cox model

# Univariate


cox_model <- coxph(
  
  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ lof_binaire,
  data=df_subset)


summary(cox_model)



# Multivariate
cox_model <- coxph(

  Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ lof_binaire +
    Age.at.recruitment + omop_gender_concept_id  +
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10,
  data=df_subset)


summary(cox_model)






### Exporting data Prism:
# Fig 5B

output <- df_subset %>% 
  mutate(follow_up_duration_years = follow_up_duration / 365.25) |> 
  select(binary_ITP_or_AIHA_indivs, follow_up_duration_years, lof_binaire)


# Write output
print('writing final table')
write.table(output,
            file = "/home/stn/Documents/GitHub/Doc/articles/first_year_phd_pubs/002_first_article/para_julie/Fig5B_high_risk_lof_only.tsv" ,
            sep = "\t", row.names = F, col.names = T , quote = F)

print('final table written')








# 
# 
# 
# ## Logistic regression
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(binary_CH_Age_sample == 1) |> 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# 
# 
# 
# summary(glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
#     omop_gender_concept_id +
#     age_end_follow_up  + 
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset, family = binomial)
#   
# )
# 
# 
# 
# 
# model_logistic_regression <- glm(
#   as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01__Hauck_ITPnHAI +
#     omop_gender_concept_id +
#     age_end_follow_up  + 
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset, family = binomial)
# 
# 
# results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
#   mutate(
#     OR = exp(estimate),
#     CI_inf = exp(conf.low),
#     CI_sup = exp(conf.high)
#   ) %>%
#   select(term, OR, CI_inf, CI_sup, p.value, statistic)
# 
# print(results, digits = 3)
# 
# 
# 
# summary_model <- summary(model_logistic_regression)
# 
# # Extract coefficients with p-values
# coef_table <- summary_model$coefficients
# 
# # View the full table with p-values (exact)
# print(coef_table, digits = 20)
# 












# 
# # ______1 SD --------------------------------------------------------------------
# 
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + omop_gender_concept_id  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 
# 





## Logistic regression


# Subset if needed
df_subset <- df_CAI_before_CH %>% 
  filter(binary_CH_Age_sample == 1) |> 
  filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH

nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))




# Subset if needed
df_subset <- df_nouvelle_colonne






summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)






# # ______1.5 SD --------------------------------------------------------------------
# 
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + omop_gender_concept_id  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 
# 





## Logistic regression


# Subset if needed
df_subset <- df_CAI_before_CH %>% 
  filter(binary_CH_Age_sample == 1) |> 
  filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH

nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1.5 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1.5 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))




# Subset if needed
df_subset <- df_nouvelle_colonne






summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)








# # ______1.75 SD --------------------------------------------------------------------
# 
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + omop_gender_concept_id  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 
# 





## Logistic regression


# Subset if needed
df_subset <- df_CAI_before_CH %>% 
  filter(binary_CH_Age_sample == 1) |> 
  filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH

nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 1.75 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 1.75 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))




# Subset if needed
df_subset <- df_nouvelle_colonne






summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)












# ______2 SD --------------------------------------------------------------------
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + omop_gender_concept_id  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 













## Logistic regression


# Subset if needed
df_subset <- df_CAI_before_CH %>% 
  filter(binary_CH_Age_sample == 1) |> 
  filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH

nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))




# Subset if needed
df_subset <- df_nouvelle_colonne






summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)













# ______3 SD --------------------------------------------------------------------
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
# 
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
#     Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + omop_gender_concept_id  +
#     Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
#     Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
#     Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
#     Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
#     Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 













## Logistic regression


# Subset if needed
df_subset <- df_CAI_before_CH %>% 
  filter(binary_CH_Age_sample == 1) |> 
  filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH

nrow(df_subset)




# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$Standard.PRS.for.systemic.lupus.erythematosus..SLE., na.rm = TRUE)



# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. < (mean_SLE_PGS + 3 * sd_SLE_PGS)  ~ "0",
    Standard.PRS.for.systemic.lupus.erythematosus..SLE. >= (mean_SLE_PGS + 3 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    WES_500k_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    WES_500k_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))




# Subset if needed
df_subset <- df_nouvelle_colonne






summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  prs_n_lof +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)


results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)








# __________ (legacy) AoU replication --------------------


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 






# Indivs with CAI diagnosis
# Annotate follow-up durations into a single column
df_follow_up_duration <- df %>% 
  mutate(follow_up_duration = ifelse(binary_ITP_or_AIHA_indivs == 1, 
                                     CAI_lower_age - CH_Age_sample_lower_age,
                                     age_end_follow_up - CH_Age_sample_lower_age))


# Convert ages from days to years BEFORE creating variables
df_follow_up_duration_years <- df_follow_up_duration %>%
  mutate(
    
    # Convert ages from days to years
    CAI_lower_age_years = CAI_lower_age / 365.25,
    CH_lower_age_lower_age_years = CH_Age_sample_lower_age / 365.25,
    age_end_follow_up_years = age_end_follow_up / 365.25,
    
    follow_up_duration_years = ifelse(binary_ITP_or_AIHA_indivs == 1,
                                      CAI_lower_age_years - CH_lower_age_lower_age_years,
                                      age_end_follow_up_years - CH_lower_age_lower_age_years)
  )



# Annotate CAI before CH
df_CAI_before_CH <- df_follow_up_duration_years %>%
  mutate(CAI_before_CH_indivs = case_when(
    CAI_lower_age < CH_Age_sample_lower_age ~ "1",
    TRUE ~ "0"
  )) 







# _____________1 SD --------------------------------------------------------------------

# 
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   binary_CH_Age_sample |> 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     PGS004917 < (mean_SLE_PGS + 1 * sd_SLE_PGS)  ~ "0",
#     PGS004917 >= (mean_SLE_PGS + 1 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     AoU_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     AoU_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# ### Cox model
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + sex_genetic + 
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 


# _____________2 SD --------------------------------------------------------------------

# 
# ### Cox model
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_CAI_before_CH %>% 
#   filter(follow_up_duration > 0) %>% 
#   binary_CH_Age_sample |> 
#   filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH
# 
# nrow(df_subset)
# 
# 
# 
# # Calculate SLE PGS mean and SD
# mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
# sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)
# 
# 
# # Stratify PRS score
# df_stratified_score <- df_subset %>% 
#   mutate(stratified_score = case_when(
#     PGS004917 < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
#     PGS004917 >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
#   ))
# 
# 
# 
# # Annotate LOF binary
# df_lof_binaire <- df_stratified_score %>% 
#   mutate(lof_binaire = case_when(
#     AoU_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
#     AoU_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
#   ))
# 
# 
# 
# # Annotate PRS binary + LoF
# df_nouvelle_colonne <- df_lof_binaire %>% 
#   mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))
# 
# 
# 
# 
# # Subset if needed
# df_subset <- df_nouvelle_colonne
# 
# 
# 
# 
# 
# # Multivariate
# cox_model <- coxph(
#   
#   Surv(follow_up_duration, binary_ITP_or_AIHA_indivs) ~ prs_n_lof +
#     Age.at.recruitment + sex_genetic + 
#     PC1 +
#     PC2 +
#     PC3 +
#     PC4 +
#     PC5 +
#     PC6 +
#     PC7 +
#     PC8 +
#     PC9 +
#     PC10, 
#   data=df_subset)
# 
# 
# summary(cox_model)
# 
# 
# 





### Logistic regression

# Subset if needed
df_subset <- df_CAI_before_CH %>% 
  filter(binary_CH_Age_sample == 1) |> 
  filter(CAI_before_CH_indivs != 1) # Remove indivs with CAI before CH

nrow(df_subset)



# Calculate SLE PGS mean and SD
mean_SLE_PGS <- mean(df_subset$PGS004917, na.rm = TRUE)
sd_SLE_PGS <- sd(df_subset$PGS004917, na.rm = TRUE)


# Stratify PRS score
df_stratified_score <- df_subset %>% 
  mutate(stratified_score = case_when(
    PGS004917 < (mean_SLE_PGS + 2 * sd_SLE_PGS)  ~ "0",
    PGS004917 >= (mean_SLE_PGS + 2 * sd_SLE_PGS ) ~"1"
  ))



# Annotate LOF binary
df_lof_binaire <- df_stratified_score %>% 
  mutate(lof_binaire = case_when(
    AoU_LoF_MAF01__Hauck_ITPnHAI <= 0.1 ~ "0",
    AoU_LoF_MAF01__Hauck_ITPnHAI > 0.1 ~"1"
  ))



# Annotate PRS binary + LoF
df_nouvelle_colonne <- df_lof_binaire %>% 
  mutate(prs_n_lof = as.numeric(lof_binaire) + as.numeric(stratified_score))




# Subset if needed
df_subset <- df_nouvelle_colonne





### Logistic regression

summary(glm(
  binary_ITP_or_AIHA_indivs ~  prs_n_lof +
    age_end_follow_up +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)



model_logistic_regression <- glm(
  binary_ITP_or_AIHA_indivs ~  prs_n_lof +
    age_end_follow_up +
    sex_genetic +
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)





results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)



summary_model <- summary(model_logistic_regression)

# Extract coefficients with p-values
coef_table <- summary_model$coefficients

# View the full table with p-values (exact)
print(coef_table, digits = 20)







# ___6) AR analysis ---------------------------------------------------------
# L’idée est donc de refaire l’analyse (régression logistique) CAI ~ LOF et covariables 
# mais en ne considérant que les gènes de la liste Hauck ITP/AHAI qui sont autosomal récessifs (AR). 


# From the 59 genes associated with ITP or AIHA in the dataset, 39 (66%) are AR, 20 are not AR 
# AR genes in dataset:  [1] "RAG1"     "DOCK8"    "LIG4"     "NCF2"     "ATM"      "ITK"      "CARD11"   "LRBA"     "DCLRE1C"  "SMARCAL1" "TFRC"
# [12] "PIK3CG"   "PIK3CD"   "STK4"     "STIM1"    "TPP2"     "RAG2"     "STAT5B"   "ARPC1B"   "ADA"      "NCF1"     "IKBKB"
# [23] "AICDA"    "ARHGEF1"  "DEF6"     "TET2"     "CTNNBL1"  "G6PC3"    "PNP"      "ZAP70"    "LCK"      "FADD"     "CD81"
# [34] "CD3G"     "CYBA"     "NHEJ1"    "LAT"      "ICOS"     "ORAI1"
# 
# Not AR genes in dataset: 
# [1] "IKZF2"     "TNFRSF13B" "NFKB1"     "PIK3R1"    "KDM6A"     "IKZF1"     "STAT3"     "FOXP3"     "SASH3"     "WAS"
# [11] "CTLA4"     "CASP10"    "CYBB"      "KMT2D"     "CD40LG"    "PTEN"      "TBX1"      "FAS"       "TNFSF12"   "SOCS1"



# Supplemental table 2

# Figure 1A


# ______ 6.1 glm cai ~ LoF AR genes only  ---------------------------------------------------------------
### Estelle made this part



library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)
library(OptSurvCutR)
library(ggsurvfit)



### Load data

df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )



df_inheritance <- fread("/home/stn/Documents/GitHub/Doc/analysis/genes_lists/IUIS_inheritance_curated_formatted.tsv")



# Filter AR genes
ar_genes <- df_inheritance %>% 
  filter(Inheritance == "AR")

# Create AR genes list
list_ar_genes <- ar_genes %>% 
  pull(gene_name) %>% 
  strsplit(";") %>% 
  unlist() %>% 
  trimws() # trim white spaces


df_lof_variants <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/002_RareVariants_allscores_bygenes_20260414_no_HC.txt", select = c("Participant.ID", "Gene", "WES_500k_LoF_MAF01"))


hauck_list <- read_xlsx("/home/stn/Documents/GitHub/Doc/analysis/genes_lists/IEI_Hauck_2024_FASLG.xlsx")





### Code

# Create list of genes associated with Hauck (ITP OR AIHA) in the dataset

# Create new columns for further filtering
df_genes_in_CAI <- df_lof_variants %>%
  mutate(itp_associated = "no",
         aiha_associated = "no" ) %>%
  select(Gene, itp_associated, aiha_associated) %>%
  distinct(Gene, .keep_all = T)



hauck_ITP <- hauck_list %>%
  filter(ITP != "") %>%
  select(ITP) %>%
  distinct(.) %>%
  pull() %>%
  strsplit(";\\s*") %>%
  unlist() %>%
  unique() %>%
  trimws()



hauck_AIHA <- hauck_list %>%
  filter(AIHA != "") %>%
  select(AIHA) %>%
  pull() %>%
  strsplit(";\\s*") %>%
  unlist() %>%
  unique() %>%
  trimws()




# Annotate genes associated with ITP
df_genes_in_CAI <- df_genes_in_CAI %>%
  mutate(
    itp_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_ITP), "yes", "no")
    })
  )


# Annotate genes associated with AIHA
df_genes_in_CAI <- df_genes_in_CAI %>%
  mutate(
    aiha_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_AIHA), "yes", "no")
    })
  )


# Filter only ITP or AIHA genes
df_genes_in_CAI <- df_genes_in_CAI %>%
  filter(itp_associated == "yes" | aiha_associated == "yes")


# Annotate AR genes
df_ar_genes <- df_genes_in_CAI %>% 
  mutate(ar_gene = ifelse(Gene %in% list_ar_genes, "yes", "no"))


# Create list of genes CAI AR in the dataset
list_dataset_cai_ar_genes <- df_ar_genes %>%
  filter(ar_gene == "yes") %>%
  pull(Gene)


df_indivs_genes_ar <- df_lof_variants %>% 
  mutate(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR = ifelse(Gene %in% list_dataset_cai_ar_genes, WES_500k_LoF_MAF01, "0")) %>% 
  select(`Participant.ID`, WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR)


# Annotate in the main table
df_main_cai_ar_annotated <- left_join(df, df_indivs_genes_ar, by = c("Participant.ID")) %>% 
  mutate(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR = ifelse(is.na(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR), 0, WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR))






# Subset if needed
df_subset <- df_main_cai_ar_annotated %>% 
  mutate(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR = as.numeric(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR))
nrow(df_subset)



### CAI

# Logistic regression
summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








### ITP

# Logistic regression
summary(glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







### AIHA

# Logistic regression
summary(glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    omop_gender_concept_id +
    age_end_follow_up  + 
    Genetic.principal.components...Array.1 + Genetic.principal.components...Array.2 + 
    Genetic.principal.components...Array.3 + Genetic.principal.components...Array.4 + 
    Genetic.principal.components...Array.5 + Genetic.principal.components...Array.6 + 
    Genetic.principal.components...Array.7 + Genetic.principal.components...Array.8 + 
    Genetic.principal.components...Array.9 + Genetic.principal.components...Array.10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








# _________ AoU replication ---------------------------------------------------------------
### Estelle made this part


required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}





df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    #    first_condition_date,
    #    first_condition_year,
    first_condition_age,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 

# Rename IID col
df <- df |> 
  rename( `Participant.ID` = person_id)




df_inheritance <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/IUIS_inheritance_curated_formatted.tsv")



# Filter AR genes
ar_genes <- df_inheritance %>% 
  filter(Inheritance == "AR")

# Create AR genes list
list_ar_genes <- ar_genes %>% 
  pull(gene_name) %>% 
  strsplit(";") %>% 
  unlist() %>% 
  trimws() # trim white spaces


df_lof_variants <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/002_RareVariants_allscores_bygenes_20260414_no_HC.txt", 
                         select = c("Participant.ID", "Gene", "WES_500k_LoF_MAF01"))


hauck_list <- read_xlsx("/home/rstudio/workspace/workspace-bucket/stennio/000_data/IEI_Hauck_2024_FASLG.xlsx")





### Code

# Create list of genes associated with Hauck (ITP OR AIHA) in df_lof_variants

# Create new columns for further filtering
df_genes_in_CAI <- df_lof_variants %>%
  mutate(itp_associated = "no",
         aiha_associated = "no" ) %>%
  select(Gene, itp_associated, aiha_associated) %>%
  distinct(Gene, .keep_all = T)



hauck_ITP <- hauck_list %>%
  filter(ITP != "") %>%
  select(ITP) %>%
  distinct(.) %>%
  pull() %>%
  strsplit(";\\s*") %>%
  unlist() %>%
  unique() %>%
  trimws()



hauck_AIHA <- hauck_list %>%
  filter(AIHA != "") %>%
  select(AIHA) %>%
  pull() %>%
  strsplit(";\\s*") %>%
  unlist() %>%
  unique() %>%
  trimws()




# Annotate genes associated with ITP
df_genes_in_CAI <- df_genes_in_CAI %>%
  mutate(
    itp_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_ITP), "yes", "no")
    })
  )


# Annotate genes associated with AIHA
df_genes_in_CAI <- df_genes_in_CAI %>%
  mutate(
    aiha_associated = sapply(Gene, function(x) {
      if(is.na(x) || x == "") return("no")
      
      # Divides cell content by ,; or -
      genes <- unlist(strsplit(as.character(x), "[,;-]"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      ifelse(any(genes %in% hauck_AIHA), "yes", "no")
    })
  )


# Filter only ITP or AIHA genes
df_genes_in_CAI <- df_genes_in_CAI %>%
  filter(itp_associated == "yes" | aiha_associated == "yes")


# Annotate AR genes
df_ar_genes <- df_genes_in_CAI %>% 
  mutate(ar_gene = ifelse(Gene %in% list_ar_genes, "yes", "no"))


# Create list of genes CAI AR in the df_lof_variants
list_dataset_cai_ar_genes <- df_ar_genes %>%
  filter(ar_gene == "yes") %>%
  pull(Gene)



### Check genes

list_ar_genes
length(list_ar_genes)


# CAI AR in Hauck
hauck_CAI_AR <- hauck_CAI[hauck_CAI %in% list_ar_genes]
hauck_CAI_AR
length(hauck_CAI_AR)


# CAI AR in df_lof_variants
list_dataset_cai_ar_genes
length(list_dataset_cai_ar_genes)


# Annotate number of variants LOF ITP AIHA AR in df_lof_variants
df_indivs_genes_ar <- df_lof_variants %>% 
  mutate(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR = ifelse(Gene %in% list_dataset_cai_ar_genes, WES_500k_LoF_MAF01, "0")) %>% 
  select(`Participant.ID`, WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR)


# Annotate in the main table
df_main_cai_ar_annotated <- left_join(df, df_indivs_genes_ar, by = c("Participant.ID")) %>% 
  mutate(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR = ifelse(is.na(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR), 0, WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR))



# Subset if needed
df_subset <- df_main_cai_ar_annotated %>% 
  mutate(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR = as.numeric(WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR))
nrow(df_subset)








### CAI

# Logistic regression
summary(glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP_or_AIHA_indivs) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)








### ITP

# Logistic regression
summary(glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_ITP) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)







### AIHA

# Logistic regression
summary(glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)
  
)




model_logistic_regression <- glm(
  as.factor(binary_AIHA) ~  WES_500k_LoF_MAF01_Hauck_ITPnHAI_AR +
    sex_genetic +
    age_end_follow_up  + 
    PC1 +
    PC2 +
    PC3 +
    PC4 +
    PC5 +
    PC6 +
    PC7 +
    PC8 +
    PC9 +
    PC10, 
  data=df_subset, family = binomial)




results <- tidy(model_logistic_regression, conf.int = TRUE, conf.level = 0.95) %>%
  mutate(
    OR = exp(estimate),
    CI_inf = exp(conf.low),
    CI_sup = exp(conf.high)
  ) %>%
  select(term, OR, CI_inf, CI_sup, p.value, statistic)

print(results, digits = 3)










# List of AR and not AR genes in the dataset

# 
# yes_ar <- df_ar_genes %>% 
#   filter(ar_gene == "yes")
# 
# list_yes_ar <- yes_ar %>% 
#   pull(Gene)
# 
# list_yes_ar
# 
# 
# 
# 
# not_ar <- df_ar_genes %>% 
#   filter(ar_gene == "no")
# 
# list_not_ar <- not_ar %>%
#   pull(Gene)
# 
# list_not_ar
# 




# ______ 6.2 AR genes in Hauck ----------------------------------------------------------------



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom", 
                       "ggplot2", 
                       "tidyplots",
                       "lubridate",
                       "survival",
                       "OptSurvCutR",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}




### Load data

df_inheritance <- fread("/home/stn/Documents/GitHub/Doc/analysis/genes_lists/IUIS_inheritance_curated_formatted.tsv")



# Filter AR genes
ar_genes <- df_inheritance %>% 
  filter(Inheritance == "AR")


# Create AR genes list
list_ar_genes <- ar_genes %>% 
  pull(gene_name) %>% 
  strsplit(";") %>% 
  unlist() %>% 
  trimws() # trim white spaces



hauck_list <- read_xlsx("/home/stn/Documents/GitHub/Doc/analysis/genes_lists/IEI_Hauck_2024_FASLG.xlsx")






### Code


hauck_ITP <- hauck_list %>%
  filter(ITP != "") %>%
  select(ITP) %>%
  distinct(.) %>%
  pull() %>%
  strsplit(";\\s*") %>%
  unlist() %>%
  unique() %>%
  trimws()

hauck_AIHA <- hauck_list %>%
  filter(AIHA != "") %>%
  select(AIHA) %>%
  pull() %>%
  strsplit(";\\s*") %>%
  unlist() %>%
  unique() %>%
  trimws()


hauck_SLE <- hauck_list %>%
  filter(`SLE/SLE-like` != "") %>%
  select(`SLE/SLE-like`) %>%
  distinct(.) %>%
  pull() %>%
  strsplit(";\\s*") %>%
  unlist() %>%
  unique() %>%
  trimws()

hauck_other_pheno  <- hauck_list %>%
  filter(other_phenotypes != "") %>%
  select(other_phenotypes) %>%
  distinct(.) %>%
  pull() %>%
  strsplit(";\\s*") %>%
  unlist() %>%
  unique() %>%
  trimws()


hauck_CAI <- c(hauck_ITP, hauck_AIHA) %>% unique()
length(hauck_CAI)

hauck_all_pheno <- c(hauck_ITP, hauck_AIHA, hauck_SLE, hauck_other_pheno)


# all pheno
hauck_all_pheno_AR <- hauck_all_pheno[hauck_all_pheno %in% list_ar_genes]
hauck_all_pheno_AR
length(hauck_all_pheno_AR)

# ITP
hauck_ITP_AR <- hauck_ITP[hauck_ITP %in% list_ar_genes]
hauck_ITP_AR
length(hauck_ITP_AR)

# AIHA
hauck_AIHA_AR <- hauck_AIHA[hauck_AIHA %in% list_ar_genes]
hauck_AIHA_AR
length(hauck_AIHA_AR)

# CAI
hauck_CAI_AR <- hauck_CAI[hauck_CAI %in% list_ar_genes]
hauck_CAI_AR
length(hauck_CAI_AR)

# SLE
hauck_SLE_AR <- hauck_SLE[hauck_SLE %in% list_ar_genes]
hauck_SLE_AR
length(hauck_SLE_AR)

# Other pheno
hauck_other_pheno_AR <- hauck_other_pheno[hauck_other_pheno %in% list_ar_genes]
hauck_other_pheno_AR
length(hauck_other_pheno_AR)







# ___Data for methods ----------------------------------------------------


library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)




df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         MYELMULT,
         
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )



# ____UKB -----------------------------------------------------------------


# _____Mean SD Age at genetic testing and last follow-up -------------------


# Mean age at genetic testing
mean_age_genetic_testing <- df %>% 
  filter(!is.na(`Age.at.recruitment`)) %>% 
  pull(`Age.at.recruitment`) %>% 
  mean()

mean_age_genetic_testing


# SD age at genetic testing
sd_age_genetic_testing <- df %>% 
  filter(!is.na(`Age.at.recruitment`)) %>% 
  pull(`Age.at.recruitment`) %>% 
  sd()

sd_age_genetic_testing


# Mean age at end of follow-up
mean_age_age_end_follow_up <- df %>% 
  filter(!is.na(age_end_follow_up)) %>% 
  pull(age_end_follow_up) %>% 
  mean() / 365.25

mean_age_age_end_follow_up

# SD age at end of follow-up
sd_age_end_follow_up <- df %>% 
  filter(!is.na(age_end_follow_up)) %>% 
  pull(age_end_follow_up) %>% 
  sd() / 365.25

sd_age_end_follow_up




# _____Indivs with ITP AIHA CAI count --------------------------------------


### ITP indivs 
itp_indivs <- df %>% 
  filter(binary_ITP == 1) 

# Count
itp_indivs_count <- itp_indivs %>% 
  nrow()

itp_indivs_count

# Mean age at diagnostic
itp_indivs_mean_age <- itp_indivs %>% 
  filter(!is.na(TPI_lower_age)) %>% 
  pull(TPI_lower_age) %>% 
  mean() / 365.25

itp_indivs_mean_age


# SD at age of diagnostic
itp_indivs_sd <- itp_indivs %>% 
  filter(!is.na(TPI_lower_age)) %>% 
  pull(TPI_lower_age) %>% 
  sd() / 365.25

itp_indivs_sd


### aiha indivs 
aiha_indivs <- df %>% 
  filter(binary_AIHA == 1) 

# Count
aiha_indivs_count <- aiha_indivs %>% 
  nrow()

aiha_indivs_count

# Mean age at diagnostic
aiha_indivs_mean_age <- aiha_indivs %>% 
  filter(!is.na(AIHA_lower_age)) %>% 
  pull(AIHA_lower_age) %>% 
  mean() / 365.25

aiha_indivs_mean_age

# SD at diagnostic
aiha_indivs_sd <- aiha_indivs %>% 
  filter(!is.na(AIHA_lower_age)) %>% 
  pull(AIHA_lower_age) %>% 
  sd() / 365.25

aiha_indivs_sd



### ES indivs 
es_indivs <- df %>% 
  filter(binary_ITP == 1 & binary_AIHA == 1) 

# Count
es_indivs_count <- es_indivs %>% 
  nrow()

es_indivs_count

# Mean age at diagnostic
es_indivs_mean_age <- es_indivs %>% 
  filter(!is.na(CAI_lower_age)) %>% 
  pull(CAI_lower_age) %>% 
  mean() / 365.25

es_indivs_mean_age


# SD at diagnostic
es_indivs_sd <- es_indivs %>% 
  filter(!is.na(CAI_lower_age)) %>% 
  pull(CAI_lower_age) %>% 
  sd() / 365.25

es_indivs_sd



### SLE indivs
sle_indivs <- df %>% 
  filter(binary_SLE == 1)

# Count
sle_indivs_count <- sle_indivs %>% 
  nrow()

sle_indivs_count

# Mean age at diagnostic
sle_indivs_mean_age <- sle_indivs %>% 
  filter(!is.na(SLE_lower_age)) %>% 
  pull(SLE_lower_age) %>% 
  mean() / 365.25

sle_indivs_mean_age

# SD at diagnostic
sle_indivs_sd <- sle_indivs %>% 
  filter(!is.na(SLE_lower_age)) %>% 
  pull(SLE_lower_age) %>% 
  sd() / 365.25

sle_indivs_sd



### LLC indivs
llc_indivs <- df %>% 
  filter(!is.na(LLC_lower_age)) 

# Count
llc_indivs_count <- llc_indivs %>% 
  nrow()

llc_indivs_count

# Mean age at diagnostic
llc_indivs_mean_age <- llc_indivs %>% 
  filter(!is.na(LLC_lower_age)) %>% 
  pull(LLC_lower_age) %>% 
  mean() / 365.25

llc_indivs_mean_age

# SD at diagnostic
llc_indivs_sd <- llc_indivs %>% 
  filter(!is.na(LLC_lower_age)) %>% 
  pull(LLC_lower_age) %>% 
  sd() / 365.25

llc_indivs_sd



### Lymphoma indivs
lymphome_indivs <- df %>% 
  filter(!is.na(LYMPHOME_lower_age))

# Count
lymphome_indivs_count <- lymphome_indivs %>% 
  nrow()

lymphome_indivs_count

# Mean age at diagnostic
lymphome_indivs_mean_age <- lymphome_indivs %>% 
  filter(!is.na(LYMPHOME_lower_age)) %>% 
  pull(LYMPHOME_lower_age) %>% 
  mean() / 365.25

lymphome_indivs_mean_age

# SD at diagnostic
lymphome_indivs_sd <- lymphome_indivs %>% 
  filter(!is.na(LYMPHOME_lower_age)) %>% 
  pull(LYMPHOME_lower_age) %>% 
  sd() / 365.25

lymphome_indivs_sd



### Other leukemia indivs
OL_indivs <- df %>% 
  filter(!is.na(Other_LEUK_lower_age))

# Count
OL_indivs_count <- OL_indivs %>% 
  nrow()

OL_indivs_count

# Mean age at diagnostic
OL_indivs_mean_age <- OL_indivs %>% 
  filter(!is.na(Other_LEUK_lower_age)) %>% 
  pull(Other_LEUK_lower_age) %>% 
  mean() / 365.25

OL_indivs_mean_age

# Mean age at diagnostic
OL_indivs_sd <- OL_indivs %>% 
  filter(!is.na(Other_LEUK_lower_age)) %>% 
  pull(Other_LEUK_lower_age) %>% 
  sd() / 365.25

OL_indivs_sd



### CH indivs
CH_indivs <- df %>% 
  filter(!is.na(CH_Age_sample_lower_age))

# Count
CH_indivs_count <- CH_indivs %>% 
  nrow()

CH_indivs_count

# Mean age at diagnostic
CH_indivs_mean_age <- CH_indivs %>% 
  filter(!is.na(CH_Age_sample_lower_age)) %>% 
  pull(CH_Age_sample_lower_age) %>% 
  mean() / 365.25

CH_indivs_mean_age


# SD at diagnostic
CH_indivs_sd <- CH_indivs %>% 
  filter(!is.na(CH_Age_sample_lower_age)) %>% 
  pull(CH_Age_sample_lower_age) %>% 
  sd() / 365.25

CH_indivs_sd



### CAI indivs
CAI_indivs <- df %>% 
  filter(!is.na(CAI_lower_age))

# Count
CAI_indivs_count <- CAI_indivs %>% 
  nrow()

CAI_indivs_count 


# Mean age at diagnostic
CAI_indivs_mean_age <- CAI_indivs %>% 
  filter(!is.na(CAI_lower_age)) %>% 
  pull(CAI_lower_age) %>% 
  mean() / 365.25

CAI_indivs_mean_age





# _____Num of indivs carrying LoF 0.1%  -------------------------------

df_indivs_1_lof <- df %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI == 1) 

df_indivs_1_lof_CAI <- df_indivs_1_lof %>% 
  filter(!is.na(CAI_lower_age))


df_indivs_2_lof <- df %>% 
  filter(WES_500k_LoF_MAF01__Hauck_ITPnHAI == 2)

df_indivs_2_lof_CAI <- df_indivs_2_lof %>% 
  filter(!is.na(CAI_lower_age))


# Indiv with 2 LoF

df_lof_variants <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/002_RareVariants_allscores_bygenes_20260414_no_HC.txt", select = c("Participant.ID", "Gene", "WES_500k_LoF_MAF01"))

indivs_2_lof <- df_lof_variants %>% 
  filter(Participant.ID == 1671995)




# _____Num of indivs carrying REVEL 0.1%  -------------------------------

## REVEL 0.9
df_CAI_indivs <- df %>% 
  filter(!is.na(CAI_lower_age))


df_CAI_indivs_carrying_1_REVEL_0p9 <- df_CAI_indivs %>% 
  filter(list_ITPnHAI_maf_MAF01_revel_0.9_count == 1)

df_CAI_indivs_carrying_2_REVEL_0p9 <- df_CAI_indivs %>% 
  filter(list_ITPnHAI_maf_MAF01_revel_0.9_count == 2)


## REVEL 0.5
df_CAI_indivs_carrying_1_REVEL_0p5 <- df_CAI_indivs %>% 
  filter(list_ITPnHAI_maf_MAF01_revel_0.5_count == 1)

df_CAI_indivs_carrying_2_REVEL_0p5 <- df_CAI_indivs %>% 
  filter(list_ITPnHAI_maf_MAF01_revel_0.5_count == 2)







# 
# #_____ No indivs CAI LOF/non-LOF -------------------------------------------------
# 
# 
# # 1 - Dans AoU, combien de patients avec CAI ont un LOF (<0.1%) et combien en ont >1?
# 
# 
# 
# library(dplyr)
# library(data.table)
# library(tidyr)
# library(stringr)
# library(readxl)
# library(broom)
# library(tidyplots)
# library(lubridate)
# library(survival)
# 
# 
# 
# 
# df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
#   select(`Participant.ID`, 
#          
#          # Diagnostics
#          SLE,
#          TPI, 
#          HAI,
#          CH_Age_sample,
#          LYMPHOME,
#          Other_LEUK,
#          LLC,
#          MYELMULT,
#          
#          
#          # Binaries
#          binary_CH_Age_sample,
#          binary_ITP_or_AIHA_indivs,
#          binary_AIHA,
#          binary_ITP,
#          binary_SLE,
#          
#          # Lower age
#          AIHA_lower_age,
#          TPI_lower_age,
#          CAI_lower_age,
#          SLE_lower_age,
#          LYMPHOME_lower_age,
#          Other_LEUK_lower_age,
#          LLC_lower_age,
#          CH_Age_sample_lower_age,
#          
#          # Survival analysis
#          Time.blood.sample.collected...Instance.0...Array.0,
#          blood_draw_date,
#          formatted_DATE_BIRTH,
#          age_inclusion,
#          date_end_follow_up,
#          age_end_follow_up,
#          span_age_end_follow_up_age_inclusion,
#          
#          # MAF < 0.1%
#          WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
#          list_ITPnHAI_maf_MAF01_revel_0.5_count,
#          list_ITPnHAI_maf_MAF01_revel_0.9_count,
#          MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
#          
#          # MAF < 1%
#          WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
#          list_ITPnHAI_maf_MAF1_revel_0.5_count,
#          list_ITPnHAI_maf_MAF1_revel_0.9_count,
#          MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
#          
#          # Covariates
#          `Age.at.recruitment`,
#          omop_gender_concept_id,
#          Genetic.principal.components...Array.1,
#          Genetic.principal.components...Array.2,
#          Genetic.principal.components...Array.3,
#          Genetic.principal.components...Array.4,
#          Genetic.principal.components...Array.5,
#          Genetic.principal.components...Array.6,
#          Genetic.principal.components...Array.7,
#          Genetic.principal.components...Array.8,
#          Genetic.principal.components...Array.9,
#          Genetic.principal.components...Array.10,
#          Standard.PRS.for.systemic.lupus.erythematosus..SLE.
#   )
# 
# 
# 
# 
# 
# ### CAI indivs
# CAI_indivs <- df %>% 
#   filter(!is.na(CAI_lower_age))
# 
# 
# # Indivs CAI and 1 LoF
# indivs_cai_1_lof <- CAI_indivs |> 
#   filter(AoU_LoF_MAF01__Hauck_ITPnHAI == 1) |> 
#   nrow()
# 
# 
# # Indivs CAI and > 1 LoF
# indivs_cai_gt_1_lof <- CAI_indivs |> 
#   filter(AoU_LoF_MAF01__Hauck_ITPnHAI > 1) |> 
#   nrow()
# 
# 






# _____ No indivs Hematological malignancies  -----------------------------




library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(readxl)
library(broom)
library(tidyplots)
library(lubridate)
library(survival)




df <- fread("/home/stn/Documents/GitHub/Doc/analysis/ukb_rare_variants/narval/trash/20260501_Summary_age_delays_diag_LOFnotHC_binaries_survival.txt")  %>% 
  select(`Participant.ID`, 
         
         # Diagnostics
         SLE,
         TPI, 
         HAI,
         CH_Age_sample,
         LYMPHOME,
         Other_LEUK,
         LLC,
         MYELMULT,
         
         
         # Binaries
         binary_CH_Age_sample,
         binary_ITP_or_AIHA_indivs,
         binary_AIHA,
         binary_ITP,
         binary_SLE,
         
         # Lower age
         AIHA_lower_age,
         TPI_lower_age,
         CAI_lower_age,
         SLE_lower_age,
         LYMPHOME_lower_age,
         Other_LEUK_lower_age,
         LLC_lower_age,
         CH_Age_sample_lower_age,
         
         # Survival analysis
         Time.blood.sample.collected...Instance.0...Array.0,
         blood_draw_date,
         formatted_DATE_BIRTH,
         age_inclusion,
         date_end_follow_up,
         age_end_follow_up,
         span_age_end_follow_up_age_inclusion,
         
         # MAF < 0.1%
         WES_500k_LoF_MAF01__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF01_revel_0.5_count,
         list_ITPnHAI_maf_MAF01_revel_0.9_count,
         MAF01_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # MAF < 1%
         WES_500k_LoF_MAF1__Hauck_ITPnHAI, # LOF only
         list_ITPnHAI_maf_MAF1_revel_0.5_count,
         list_ITPnHAI_maf_MAF1_revel_0.9_count,
         MAF1_LoF_not_CH_or_revel_0.9_ITPnHAI,
         
         # Covariates
         `Age.at.recruitment`,
         omop_gender_concept_id,
         Genetic.principal.components...Array.1,
         Genetic.principal.components...Array.2,
         Genetic.principal.components...Array.3,
         Genetic.principal.components...Array.4,
         Genetic.principal.components...Array.5,
         Genetic.principal.components...Array.6,
         Genetic.principal.components...Array.7,
         Genetic.principal.components...Array.8,
         Genetic.principal.components...Array.9,
         Genetic.principal.components...Array.10,
         Standard.PRS.for.systemic.lupus.erythematosus..SLE.
  )



# HM indivs
hm_indivs <- df |> 
  filter(!is.na(LLC_lower_age) | 
           !is.na(LYMPHOME_lower_age)| 
           !is.na(Other_LEUK_lower_age))




# Mean age at diagnostic
hm_indivs_mean_age <- hm_indivs %>% 
  mutate(hm_lower_age = pmin(LLC_lower_age,
                             LYMPHOME_lower_age,
                             Other_LEUK_lower_age, 
                             na.rm = T)) |> 
  pull(hm_lower_age) |> 
  mean() / 365.25
    

hm_indivs_mean_age



# SD age at genetic testing
sd_age_genetic_testing <- hm_indivs %>% 
  mutate(hm_lower_age = pmin(LLC_lower_age,
                             LYMPHOME_lower_age,
                             Other_LEUK_lower_age, 
                             na.rm = T)) |> 
  pull(hm_lower_age) |> 
  sd() / 365.25

sd_age_genetic_testing










# ____AoU -----------------------------------------------------------------




df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 


# _____Mean SD Age at genetic testing and last follow-up -------------------


# Mean age at genetic testing
mean_age_genetic_testing <- df %>% 
  filter(!is.na(`Age.at.recruitment`)) %>% 
  pull(`Age.at.recruitment`) %>% 
  mean()

mean_age_genetic_testing


# SD age at genetic testing
sd_age_genetic_testing <- df %>% 
  filter(!is.na(`Age.at.recruitment`)) %>% 
  pull(`Age.at.recruitment`) %>% 
  sd()

sd_age_genetic_testing


# Mean age at end of follow-up
mean_age_age_end_follow_up <- df %>% 
  filter(!is.na(age_end_follow_up)) %>% 
  pull(age_end_follow_up) %>% 
  mean() / 365.25

mean_age_age_end_follow_up

# SD age at end of follow-up
sd_age_end_follow_up <- df %>% 
  filter(!is.na(age_end_follow_up)) %>% 
  pull(age_end_follow_up) %>% 
  sd() / 365.25

sd_age_end_follow_up




# _____Indivs with ITP AIHA CAI count --------------------------------------


### ITP indivs 
itp_indivs <- df %>% 
  filter(binary_ITP == 1)

# Count
itp_indivs_count <- itp_indivs %>% 
  nrow()

itp_indivs_count

# Mean age at diagnostic
itp_indivs_mean_age <- itp_indivs %>% 
  filter(!is.na(TPI_lower_age)) %>% 
  pull(TPI_lower_age) %>% 
  mean() / 365.25

itp_indivs_mean_age


# SD at age of diagnostic
itp_indivs_sd <- itp_indivs %>% 
  filter(!is.na(TPI_lower_age)) %>% 
  pull(TPI_lower_age) %>% 
  sd() / 365.25

itp_indivs_sd


### aiha indivs 
aiha_indivs <- df %>% 
  filter(binary_AIHA == 1) 

# Count
aiha_indivs_count <- aiha_indivs %>% 
  nrow()

aiha_indivs_count

# Mean age at diagnostic
aiha_indivs_mean_age <- aiha_indivs %>% 
  filter(!is.na(AIHA_lower_age)) %>% 
  pull(AIHA_lower_age) %>% 
  mean() / 365.25

aiha_indivs_mean_age

# SD at diagnostic
aiha_indivs_sd <- aiha_indivs %>% 
  filter(!is.na(AIHA_lower_age)) %>% 
  pull(AIHA_lower_age) %>% 
  sd() / 365.25

aiha_indivs_sd



### ES indivs 
es_indivs <- df %>% 
  filter(binary_ITP == 1 & binary_AIHA == 1) 

# Count
es_indivs_count <- es_indivs %>% 
  nrow()

es_indivs_count

# Mean age at diagnostic
es_indivs_mean_age <- es_indivs %>% 
  filter(!is.na(CAI_lower_age)) %>% 
  pull(CAI_lower_age) %>% 
  mean() / 365.25

es_indivs_mean_age


# SD at diagnostic
es_indivs_sd <- es_indivs %>% 
  filter(!is.na(CAI_lower_age)) %>% 
  pull(CAI_lower_age) %>% 
  sd() / 365.25

es_indivs_sd



### SLE indivs
sle_indivs <- df %>% 
  filter(binary_SLE == 1)

# Count
sle_indivs_count <- sle_indivs %>% 
  nrow()

sle_indivs_count

# Mean age at diagnostic
sle_indivs_mean_age <- sle_indivs %>% 
  filter(!is.na(SLE_lower_age)) %>% 
  pull(SLE_lower_age) %>% 
  mean() / 365.25

sle_indivs_mean_age

# SD at diagnostic
sle_indivs_sd <- sle_indivs %>% 
  filter(!is.na(SLE_lower_age)) %>% 
  pull(SLE_lower_age) %>% 
  sd() / 365.25

sle_indivs_sd



### LLC indivs
llc_indivs <- df %>% 
  filter(!is.na(LLC_lower_age)) 

# Count
llc_indivs_count <- llc_indivs %>% 
  nrow()

llc_indivs_count

# Mean age at diagnostic
llc_indivs_mean_age <- llc_indivs %>% 
  filter(!is.na(LLC_lower_age)) %>% 
  pull(LLC_lower_age) %>% 
  mean() / 365.25

llc_indivs_mean_age

# SD at diagnostic
llc_indivs_sd <- llc_indivs %>% 
  filter(!is.na(LLC_lower_age)) %>% 
  pull(LLC_lower_age) %>% 
  sd() / 365.25

llc_indivs_sd



### Lymphoma indivs
lymphome_indivs <- df %>% 
  filter(!is.na(LYMPHOME_lower_age))

# Count
lymphome_indivs_count <- lymphome_indivs %>% 
  nrow()

lymphome_indivs_count

# Mean age at diagnostic
lymphome_indivs_mean_age <- lymphome_indivs %>% 
  filter(!is.na(LYMPHOME_lower_age)) %>% 
  pull(LYMPHOME_lower_age) %>% 
  mean() / 365.25

lymphome_indivs_mean_age

# SD at diagnostic
lymphome_indivs_sd <- lymphome_indivs %>% 
  filter(!is.na(LYMPHOME_lower_age)) %>% 
  pull(LYMPHOME_lower_age) %>% 
  sd() / 365.25

lymphome_indivs_sd



### Other leukemia indivs
OL_indivs <- df %>% 
  filter(!is.na(Other_LEUK_lower_age))

# Count
OL_indivs_count <- OL_indivs %>% 
  nrow()

OL_indivs_count

# Mean age at diagnostic
OL_indivs_mean_age <- OL_indivs %>% 
  filter(!is.na(Other_LEUK_lower_age)) %>% 
  pull(Other_LEUK_lower_age) %>% 
  mean() / 365.25

OL_indivs_mean_age

# Mean age at diagnostic
OL_indivs_sd <- OL_indivs %>% 
  filter(!is.na(Other_LEUK_lower_age)) %>% 
  pull(Other_LEUK_lower_age) %>% 
  sd() / 365.25

OL_indivs_sd



### CH indivs
CH_indivs <- df %>% 
  filter(!is.na(CH_Age_sample_lower_age))

# Count
CH_indivs_count <- CH_indivs %>% 
  nrow()

CH_indivs_count

# Mean age at diagnostic
CH_indivs_mean_age <- CH_indivs %>% 
  filter(!is.na(CH_Age_sample_lower_age)) %>% 
  pull(CH_Age_sample_lower_age) %>% 
  mean() / 365.25

CH_indivs_mean_age


# SD at diagnostic
CH_indivs_sd <- CH_indivs %>% 
  filter(!is.na(CH_Age_sample_lower_age)) %>% 
  pull(CH_Age_sample_lower_age) %>% 
  sd() / 365.25

CH_indivs_sd



### CAI indivs
CAI_indivs <- df %>% 
  filter(!is.na(CAI_lower_age))

# Count
CAI_indivs_count <- CAI_indivs %>% 
  nrow()

CAI_indivs_count 


# Mean age at diagnostic
CAI_indivs_mean_age <- CAI_indivs %>% 
  filter(!is.na(CAI_lower_age)) %>% 
  pull(CAI_lower_age) %>% 
  mean() / 365.25

CAI_indivs_mean_age











#_____ No indivs CAI LOF/non-LOF -------------------------------------------------


# 1 - Dans AoU, combien de patients avec CAI ont un LOF (<0.1%) et combien en ont >1?



required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}



df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 




### CAI indivs
CAI_indivs <- df %>% 
  filter(!is.na(CAI_lower_age))


# Indivs CAI and 1 LoF
indivs_cai_1_lof <- CAI_indivs |> 
  filter(AoU_LoF_MAF01__Hauck_ITPnHAI == 1) |> 
  nrow()


# Indivs CAI and > 1 LoF
indivs_cai_gt_1_lof <- CAI_indivs |> 
  filter(AoU_LoF_MAF01__Hauck_ITPnHAI > 1) |> 
  nrow()








# _____ No indivs Hematological malignancies  -----------------------------




required_packages <- c("data.table", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "readxl", 
                       "broom",
                       "lubridate",
                       "survival",
                       "ggsurvfit")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}



df <- fread("/home/rstudio/workspace/workspace-bucket/stennio/000_data/20260619_AoU_master_table.tsv") %>% 
  select(
    # General data
    person_id,
    date_end_follow_up,
    age_end_follow_up,
    sex_genetic,
    year_of_birth,
    biosample_collection_age,
    span_age_end_follow_up_age_inclusion,
    Age.at.recruitment,
    
    # Diagnostics
    TPI,
    HAI,
    SLE,
    LLC,
    LYMPHOME,
    Other_LEUK,
    CH_Age_sample,
    
    # Binaries 
    binary_ITP_or_AIHA_indivs,
    binary_AIHA,
    binary_ITP,
    binary_SLE,
    binary_LLC,
    binary_LYMPHOME,
    binary_Other_LEUK,
    binary_CH_Age_sample,
    
    # Lower age
    AIHA_lower_age,
    TPI_lower_age,
    CAI_lower_age,
    SLE_lower_age,
    LYMPHOME_lower_age,
    Other_LEUK_lower_age,
    LLC_lower_age,
    CH_Age_sample_lower_age,
    
    # Covariates
    AoU_LoF_MAF01__Hauck_ITPnHAI,
    PGS000196,
    PGS004917,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10
  ) 


# HM indivs
hm_indivs <- df |> 
  filter(!is.na(LLC_lower_age) | 
           !is.na(LYMPHOME_lower_age)| 
           !is.na(Other_LEUK_lower_age))




# Mean age at diagnostic
hm_indivs_mean_age <- hm_indivs %>% 
  mutate(hm_lower_age = pmin(LLC_lower_age,
                             LYMPHOME_lower_age,
                             Other_LEUK_lower_age, 
                             na.rm = T)) |> 
  pull(hm_lower_age) |> 
  mean() / 365.25
    

hm_indivs_mean_age



# SD age at genetic testing
sd_age_genetic_testing <- hm_indivs %>% 
  mutate(hm_lower_age = pmin(LLC_lower_age,
                             LYMPHOME_lower_age,
                             Other_LEUK_lower_age, 
                             na.rm = T)) |> 
  pull(hm_lower_age) |> 
  sd() / 365.25

sd_age_genetic_testing



