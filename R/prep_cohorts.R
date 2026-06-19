# Combine log files -------------------------------------------------------


# Combining WES Data, cohort 2:
setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")
library(dplyr)
library(readr)

# Prep cohort 2
part1 <- read.csv('../backups/wes-mrn-cohort2_part1.csv')
part2 <- read.csv('../backups/wes-mrn-cohort2_part2.csv')

cohort1 <- read.csv('../backups/wes-mrn-cohort1.csv')

# part1 <- part1 %>% mutate(duplicated = mrn %in% mrn[duplicated(mrn)])

part2 <- part2 %>% dplyr::rename(DNA.ID = Collaborator.Participant.ID) %>%
  dplyr::select(Sample.ID,DNA.ID,Root.Sample.s.)

# Checks to ensure the sets match:
sum(part1$DNA.ID %in% part2$DNA.ID)
sum(part2$DNA.ID %in% part1$DNA.ID)

cohort2 <- left_join(part1,part2,by='DNA.ID')


cohort1 <- cohort1 %>% mutate(cohort = "old")
cohort2 <- cohort2 %>% mutate(cohort = "new")


# Combine cohort 2:
# all_cols <- union(names(cohort1), names(cohort2))
# cohort1_full <- cohort1 %>% mutate(across(setdiff(all_cols, names(cohort1)), ~NA))
# cohort2_full <- cohort2 %>% mutate(across(setdiff(all_cols, names(cohort2)), ~NA))
combined <- bind_rows(cohort1, cohort2)


write_csv(cohort2, "../backups/wes-mrn-cohort2.csv")


# read demo file, add record_id

demo <- read.csv('/mmfs1/projects/cardio_darbar_chi/common/data/demographics-0.csv')
demo <- demo %>% select(mrn,record_id) %>% mutate(mrn=as.numeric(mrn))

combined <- left_join(combined,demo,by='mrn')

write_csv(combined, '../wes-mrn-combined.csv')

# sample ID: 
# batch1: 
    # name template: "CCDG_Broad_CVD_AF_Darbar_UIC_Cases-UIC0757"
    # corresponding col: study_id (1, 4, 5 etc)
# batch 2: 
    # name template: "SM-OJRMK"
    # corresponding col: Sample.ID (SM-OJT3Y, SM-OJT3Z etc)

# check set 1 ------------------------------------------------------------------
library(arrow)
library(dplyr)
library(stringr)

table <- read.csv('wes-mrn-combined.csv') # coordinating file
vep_file <- read_parquet('vep-combined.parquet')

vep_file <- vep_file %>%
  mutate(
    study_id = str_extract(sample_id, "\\d{4}$") |> as.integer()
  )

vep_unique <- vep_file %>% filter(batch == 1) %>% select(sample_id,study_id,batch) %>% unique()

table <- table %>% filter(cohort == 'old')

sum(vep_unique$study_id %in% table$study_id) # 101 IDs matching in VEP table, log table 
  # There are unaccounted for study_ids in BOTH datasets
sum(table$study_id %in% vep_unique$study_id )

# ~300 VEP/VCF genetics files
# ~ 300 rows in table/log 


# check set 2 -------------------------------------------------------------
table <- read.csv('wes-mrn-combined.csv') # coordinating file
vep_file <- read_parquet('vep-combined.parquet')

vep_file

vep_unique <- vep_file %>% filter(batch == 2) %>% select(sample_id,batch) %>% unique()

table <- table %>% filter(cohort == 'new')

sum(vep_unique$sample_id %in% table$Sample.ID) # missing only 2 Sample.ID from the table sheet

sum(table$Sample.ID %in% vep_unique$sample_id)



# Check set 1 files -------------------------------------------------------

# Path
path <- "../../../data/genetics/uic_first_batch/vep/"

# List all files in the directory
files <- list.files(path, full.names = TRUE)

# Extract trailing 4 digits before the file extension
ids <- sub(".*UIC(\\d{4}).*", "\\1", files)

# Remove leading zeros
ids_num <- as.integer(ids)

# Create named vector: names = file paths, values = numeric IDs
result <- setNames(ids_num, files)

sum(table$study_id %in% result)
sum(result %in%table$study_id)


# Path
path <- "../../../data/genetics/uic_first_batch/vep_filtered/"

# List all files in the directory
files <- list.files(path, full.names = TRUE)

# Extract trailing 4 digits before the file extension
ids <- sub(".*UIC(\\d{4}).*", "\\1", files)

# Remove leading zeros
ids_num <- as.integer(ids)

# Create named vector: names = file paths, values = numeric IDs
result_filt <- setNames(ids_num, files)




# Path
path <- "../../../data/genetics/uic_first_batch/vcf/"

# List all files in the directory
files <- list.files(path, full.names = TRUE)

# Extract trailing 4 digits before the file extension
ids <- sub(".*UIC(\\d{4}).*", "\\1", files)

# Remove leading zeros
ids_num <- as.integer(ids)

# Create named vector: names = file paths, values = numeric IDs
result_vcf <- setNames(ids_num, files)



# Path
path <- "../../../data/genetics/uic_first_batch/lof/"

# List all files in the directory
files <- list.files(path, full.names = TRUE)

# Extract trailing 4 digits before the file extension
ids <- sub(".*UIC(\\d{4}).*", "\\1", files)

# Remove leading zeros
ids_num <- as.integer(ids)

# Create named vector: names = file paths, values = numeric IDs
result_lof <- setNames(ids_num, files)


# Check set 2 files -------------------------------------------------------
library(arrow)
library(dplyr)
library(readr)
library(stringr)

options(tibble.print_max = Inf)
options(tibble.width = Inf)

# Sample.ID # for batch 2
# study_id # for batch 1


path <- "../../../data/genetics/uic_second_batch/vep/"
files <- list.files(path, full.names = TRUE)

wes <- read.csv('wes-mrn-combined.csv') %>% filter(cohort=='new')



path <- "../../../data/genetics/uic_second_batch/vep/"
files <- list.files(path, full.names = TRUE)
# Extract the prefix before ".hard-filtered.vep"
batch2_filenames <- sub(".*/(.*)\\.hard-filtered\\.vep$", "\\1", files)

sum(wes$Sample.ID %in% batch2_filenames)
sum(batch2_filenames %in% wes$Sample.ID)

length(batch2_filenames)

vep <- read_parquet('vep-combined.parquet') %>% filter(batch == 2) %>% select(sample_id)


# pull patient/ecg data -------------------------------------------------------------------------

# pull patient data for all record_ids in wes-mrn-combined.csv 
# pull ECGs for all mrns/record_ids in wes-mrn-combined.csv 


# build ECG table ---------------------------------------------------------

# Build ECG dataframe
library(EGM)

# Parameters
home_dir <- fs::path('/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes_ml/pt-data/')
input_folder <- fs::path(home_dir,"wfdb")
output_file <- fs::path(home_dir,'ecg_table.parquet')

wfdb_table <- data.frame(ecg_name = tools::file_path_sans_ext(list.files(input_folder,pattern='*.hea')))

wfdb_table$null_ecg <- 0
for (i in 1:nrow(wfdb_table)) {
  sig <- read_signal(record = wfdb_table$ecg_name[i],record_dir = 'wfdb')
  if (all(sig[1,-1] == 0)) {
    wfdb_table$null_ecg[i] <- 1
  }
  if (!i%%500) {
    print(paste(i,'in',nrow(wfdb_table)))
  }
}

wfdb_table$freq_250 <- 0
for (i in 1:nrow(wfdb_table)) {
  hea <- read_header(record = wfdb_table$ecg_name[i],record_dir = 'wfdb')
  if (attributes(hea)$record_line$frequency !=500) {
    wfdb_table$freq_250[i] <- 1
  }
  if (!i%%500) {
    print(paste(i,'in',nrow(wfdb_table)))
  }
}

# Add mrn col:
library(dplyr)
library(stringr)
library(arrow)
library(lubridate)


wfdb_table <- wfdb_table %>%
  mutate(mrn = as.numeric(str_extract(ecg_name, "^[^_]+")))

# Add in record_id
demo <- read_parquet('/mmfs1/projects/cardio_darbar_chi/common/data/demographics-0.parquet')
demo <- demo %>% mutate(mrn = as.numeric(mrn)) %>% select(mrn,record_id)
wfdb_table <- wfdb_table %>% left_join(demo,by='mrn')

# Add in age:
wfdb_table <- wfdb_table %>%
  mutate(
    ecg_date = str_extract(ecg_name, "(?<=_)\\d{8}(?=_)"),
    ecg_date = as.Date(ecg_date, format = "%Y%m%d")
  )


wfdb_table <- wfdb_table %>%
  left_join(demo %>% select(record_id, birth_date), by = "record_id") %>%
  
  # 3. Compute age at ECG
  mutate(
    ecg_age = round(time_length(interval(birth_date, ecg_date), unit = "years"),2)
  )





write_parquet(wfdb_table,output_file)

# (Run verify sinus)


# Find HR for ECG dataframe ------------------------------------------------------------------
library(EGM)
library(arrow)
library(dplyr)

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")

common <- fs::path('/mmfs1/projects/cardio_darbar_chi/common')
source(fs::path(common,'software/cluster_general_software/AF_DM_functions.R'))
source(fs::path(common,'software/cluster_general_software/annotator_prep_functions.R'))
source(fs::path(common,'software/cluster_general_software/pvc_functions_TEMP.R'))



ecg <- read_parquet('../ecg_table.parquet')
wes <- read.csv('../wes-mrn-combined.csv')
ecg$HR <- NA

start <- 1
for (i in start:length(ecg$ecg_name)) {
  wfdb <- read_wfdb(record = ecg$ecg_name[i],record_dir = '../wfdb')
  ann <- wfdb$annotation
  fs <- attributes(wfdb$header)$record_line$frequency
  
  
  filt <- ecg_filter(wfdb$signal[['II']])
  rpeaks <- EGM::detect_QRS(filt,frequency = fs)
  
  time_length <- length(filt) / fs
  rate <- length(rpeaks) / time_length * 60
  
  ecg$HR[i] <- rate
  
  if (!i%%200) {print(paste(i,'of',nrow(ecg)))}
}


# join <- ecg %>% left_join(wes %>% select(mrn,ttn_present,ttn_mutant))
# global_avg_hr <- join %>% dplyr::filter(ttn_present == TRUE & ttn_mutant == TRUE) %>%
#   dplyr::filter(ecg_sinus == TRUE) %>%
#   group_by(record_id) %>%
#   summarise(mean_hr = mean(HR, na.rm = TRUE)) %>%
#   summarise(global_mean_hr = mean(mean_hr, na.rm = TRUE)) %>%
#   pull(global_mean_hr)
# print(global_avg_hr)

write_parquet(ecg,'../ecg_table.parquet')


# Add AF dx to demo -------------------------------------------------------
library(arrow)
library(dplyr)

options(tibble.print_max = Inf)
options(tibble.width = Inf)

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")

common <- fs::path('/mmfs1/projects/cardio_darbar_chi/common')
source(fs::path(common,'software/cluster_general_software/AF_DM_functions.R'))
source(fs::path(common,'software/cluster_general_software/annotator_prep_functions.R'))
source(fs::path(common,'software/cluster_general_software/pvc_functions_TEMP.R'))



ecg <- read_parquet('../ecg_table.parquet')
log <- read.csv('../wes-mrn-combined.csv')
demo <- read_parquet('../demo.parquet')

af_icds <- c('I48')

icdRecid2date_outpt <- function(table,icds,name,option = 'all') {
  
  # Input: table of recids; ICD code/s; name of new column
  # Output: table with column for date of initial dx
  # Note: can handle multiple ICDs for single dx
  
  #   Options: 
    #   'all': diagnoses from ALL providers, including ED
    #   'medicine': diagnoses from IM/FM/cards (non-ED)
    #   'cards': diagnoses from cards providers only (non-ED)
  
  icds <- paste0(icds, collapse = "|") 
  
  if (any(names(table) %in% 'recid')) {
    recid_colname <- 'recid'
  } else {recid_colname <- 'record_id'} 
  
  
 
  all_years <- c()
  for (year in 2010:2025) {
    
    home <- fs::path('/mmfs1/projects/cardio_darbar_chi/common/data/')
    dxfile <- fs::path(home,paste0('Diagnosis/Diagnosis_',year),ext = 'parquet')
    visfile <- fs::path(home,paste0('Visits/Visits_',year),ext = 'parquet')
    
    dx <- arrow::open_dataset(dxfile) |>
      dplyr::filter(grepl(x = icd_code, pattern = icds) & 
                      record_id %in% table[[recid_colname]]) |>
      select(record_id, date,encounter_id) |>
      collect()
    
    vis <- arrow::open_dataset(visfile) |> select(-record_id) |> collect()
    
    if (option == 'cards_only') { # non-ED, cardiologist only
      include_dpt <- 'CARD|HEART'
      exclude_dpt <- 'SURG|PLANNING|CT'
      
      vis <- vis %>% filter(visit_location %in% c('Emergency to Inpatient','Inpatient','Outpatient','Other')) %>%
        filter(grepl(pattern = include_dpt, x = department)) %>%
        filter(!grepl(pattern = exclude_dpt, x = department))
      
    } else if (option == 'medicine') { # non-ED, cardiologist, IM, FM
      include_dpt <- 'CARD|PRIMARY|HEART|FAMILY|INTERNAL\\ MEDICINE|MS-HLTH'
      exclude_dpt <- 'SURG|PLANNING|CT'
      
      vis <- vis %>% filter(visit_location %in% c('Emergency to Inpatient','Inpatient','Outpatient','Other')) %>%
        filter(grepl(pattern = include_dpt, x = department)) %>%
        filter(!grepl(pattern = exclude_dpt, x = department))
      
    } else if (option == 'all') { # any provider, including ED
      exclude_dpt <- 'SURG|PLANNING|CT'
      
      vis <- vis %>% filter(!grepl(pattern = exclude_dpt, x = department))
    }
    
    dx <- inner_join(dx, vis, by = "encounter_id")
    
    # dx <- left_join(dx,vis,by='encounter_id') |>
    #   filter(!grepl("Emergency|Inpatient", visit_location)) |>
    #   filter(grepl(include_dpt, department)) |>
    #   filter(!grepl(exclude_dpt, department))
    
    
    print(paste0('Finished ', year, '. ', nrow(dx)))
    all_years <- rbind(all_years,dx)
    rm(dx,vis)
  }
  
  min_dates <- all_years %>%
    group_by(record_id) %>%
    summarize(date = min(date, na.rm = TRUE))  # Get min date per record_id
  
  if (recid_colname == 'recid') {
    min_dates <- min_dates %>% rename(recid = record_id) # rename
  }
  
  min_dates <- min_dates %>% rename(!!name := date)
  
  final_table <- table %>%
    left_join(min_dates, by = recid_colname)  # Perform a left join on 'rec_id'
  
  return(final_table)
}

demo <- icdRecid2date_outpt(demo,af_icds,'af_dx_any','all')
demo <- icdRecid2date_outpt(demo,af_icds,'af_dx_medicine','medicine')
demo <- icdRecid2date_outpt(demo,af_icds,'af_dx_cards','cards_only')

write_parquet(demo,'../demo.parquet')

demo <- demo %>% select(-age,-ttn_present,-ttn_mutant)

log <- log %>% 
  select(-record_id) %>%
  rename(sex_logdata = sex,
         ethnicity_logdata = ethnicity) %>%
  left_join(demo,by='mrn')

# Remove two duplicated lines
log <- log %>% filter(!duplicated(.))

write_parquet(log,'../wes-mrn-combined.csv')


# Build VEP file --------------------------------------------------------
library(arrow)
library(dplyr)
library(stringr)

options(tibble.print_max = Inf)
options(tibble.width = Inf)

# Recommend 16 GB RAM

vep1 = read.csv('../../../data/genetics/uic_first_batch/vep_filtered/vep_annotations.csv')
vep1 <- vep1 %>% filter(SYMBOL == 'TTN')

vep2 = read.csv('../../../data/genetics/uic_second_batch/vep_filtered/vep_annotations.csv')
vep2 <- vep2 %>% filter(SYMBOL == 'TTN')

vep_combined <- rbind(vep1,vep2)
write_parquet(vep_combined,'vep-combined.parquet')


# wes-mrn-combined: add if TTN is present ---------------------------------------------------------
library(arrow)
library(dplyr)
library(stringr)
library(readr)

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")

options(tibble.print_max = Inf)
options(tibble.width = Inf)

wes_log <- read.csv('../wes-mrn-combined.csv')
vep <- read_parquet('../vep-combined.parquet')

# Sample.ID # for batch 2
# study_id # for batch 1

vep_unique <- vep %>% distinct(sample_id, batch)

vep_unique <- vep_unique %>% 
  filter(!grepl(x=sample_id,pattern='UIC_Cases-AC')) %>% # remove files with AC in the name
  mutate(sample_id = sub(".*UIC(\\d{4}).*", "\\1", sample_id)) %>%
  
  mutate( # mutate the 1st cohort sample_ids, which only contain numbers
    sample_id = if_else(
      grepl("^0*\\d+$", sample_id),   # only digits (with possible leading zeros)
      as.character(as.integer(sample_id)),  # strip leading zeros
      sample_id                       # otherwise leave unchanged
    )
  )


old_ids <- vep_unique %>%
  filter(batch == 1) %>%
  pull(sample_id)

new_ids <- vep_unique %>%
  filter(batch == 2) %>%
  pull(sample_id)

wes_log <- wes_log %>%
  mutate(
    ttn_present = case_when(
      cohort == "old" & study_id %in% old_ids ~ TRUE,
      cohort == "new" & Sample.ID %in% new_ids ~ TRUE,
      TRUE ~ FALSE
    )
  )

write_csv(wes_log, '../wes-mrn-combined.csv')

# Cohort 1: poor overlap between VEP files vs. log file (~300 vs ~300, with ~100 overlap)
# Cohort 2: good overlap between VEP files vs. log file (~240 vs ~250, with 238 overlap)


# Find pathogenic variants --------------------------------------------------------
library(arrow)
library(dplyr)
library(stringr)

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")

options(tibble.print_max = Inf)
options(tibble.width = Inf)

demo <- read_parquet('../demographics-0.parquet')
wes_log <- read.csv('../wes-mrn-combined.csv')
ecg <- read_parquet('../ecg_table.parquet')
vep <- read_parquet('../vep-combined.parquet')

demo <- demo %>% mutate(mrn = as.numeric(mrn))

# vep <- vep %>%
#   mutate(
#     study_id = str_extract(sample_id, "\\d+$") |> as.integer()
#   )


# Label mutant record_ids:
ttn_mutant <-
  vep |>
  filter(CANONICAL == "YES") |>
  filter(
    !is.na(LoF) |
      IMPACT == "HIGH|MODERATE" |
      str_detect(CLIN_SIG, "pathogenic") |
      str_detect(
        Consequence,
        "stop_gained|missense_variant|frameshift_variant|splice_donor_variant|splice_acceptor_variant|start_lost"
      ) |
      str_detect(SIFT, "deleterious") |
      str_detect(PolyPhen, "damaging")
  ) |>
  filter(MAX_AF < 0.01) |>
  filter(!str_detect(CLIN_SIG, "benign"))



vep_unique <- ttn_mutant %>% distinct(sample_id, batch)

vep_unique <- vep_unique %>% 
  filter(!grepl(x=sample_id,pattern='UIC_Cases-AC')) %>% # remove files with AC in the name
  mutate(sample_id = sub(".*UIC(\\d{4}).*", "\\1", sample_id)) %>%
  
  mutate( # mutate the 1st cohort sample_ids, which only contain numbers
    sample_id = if_else(
      grepl("^0*\\d+$", sample_id),   # only digits (with possible leading zeros)
      as.character(as.integer(sample_id)),  # strip leading zeros
      sample_id                       # otherwise leave unchanged
    )
  )


old_ids <- vep_unique %>%
  filter(batch == 1) %>%
  pull(sample_id)

new_ids <- vep_unique %>%
  filter(batch == 2) %>%
  pull(sample_id)

wes_log <- wes_log %>%
  mutate(
    ttn_mutant = case_when(
      cohort == "old" & study_id %in% old_ids ~ TRUE,
      cohort == "new" & Sample.ID %in% new_ids ~ TRUE,
      TRUE ~ FALSE
    )
  )

# write_csv(wes_log, '../wes-mrn-combined.csv')
write_parquet(wes_log,'../wes-mrn-combined.parquet')
# wes_log %>% select(-mrn,-sex,-ethnicity,-age,-Root.Sample.s.)

wes_log <- wes_log %>% select(record_id,ttn_present,ttn_mutant)

demo <- demo %>% left_join(wes_log,by='record_id')

write_parquet(demo,'../demo.parquet')
write_parquet(demo,'../demo.parquet')



# Use only first batch for now 
# batch1 <- vep %>% filter(batch == 1) # all WES
# ttn_batch1 <- ttn_mutant %>% filter(batch == 1) # TTN only WES

# wes_table <- wes_table %>% filter(cohort == 'old') # WES IDs and MRNs
demo <- demo %>% filter(mrn %in% wes_log$mrn) # demographic data
ecg <- ecg %>% filter(mrn %in% wes_log$mrn) # all ECGs

demo$ttn_mutant <- NA

# in ttn variable: convert sample_id col: CCDG_Broad_CVD_AF_Darbar_UIC_Cases-UIC0002 to 2, with name study_id
# then add col to wes_table for 'ttn_mutant' based on study_id
# then add 'ttn_mutant' col to demo based on mrn from wes_table
# can also add 'ttn_mutant' col to ecg



# Add de-identified (ID analgous to MRN)----------------------------------------------------
library(arrow)
library(dplyr)

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")

options(tibble.print_max = Inf)
options(tibble.width = Inf)

demo <- read_parquet('../demo.parquet')
ecg <- read_parquet('../ecg_table.parquet')
log <- read_parquet('../wes-mrn-combined.parquet')

# Add a new numbering system to keep track of unqiue patients (in case user wants to use multiple ECGs per pt)
log <- log %>%
  mutate(pt_number = dense_rank(mrn))
mini_log <- log %>% select(mrn,pt_number)

ecg <- ecg %>% left_join(mini_log,by='mrn')
demo <- demo %>% mutate(mrn = as.numeric(mrn)) %>%
  left_join(mini_log,by='mrn')

write_parquet(demo,'../demo.parquet')
write_parquet(ecg,'../ecg_table.parquet')
write_parquet(log,'../wes-mrn-combined.parquet')

# add unique number for each ECG ------------------------------------------
library(arrow)
library(dplyr)

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")

options(tibble.print_max = Inf)
options(tibble.width = Inf)

ecg <- read_parquet('../ecg_table.parquet')
ecg$ecg_number <- 1:nrow(ecg)

write_parquet(ecg,'../ecg_table.parquet')


# create deidentified ECG folder ------------------------------------------
library(arrow)
library(dplyr)
library(EGM)

options(tibble.print_max = Inf)
options(tibble.width = Inf)

ecg <- read_parquet('ecg_table.parquet')

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn")

deid_wfdb <- function(old_name,old_dir,new_name,new_dir) {
  
  wfdb <- read_wfdb(record = old_name,record_dir = old_dir,annotator = 'ann')
  
  attributes(wfdb$header)$record_line$record_name <- new_name
  wfdb$header$file_name <- rep(paste0(new_name,'.dat'),12)
  
  write_wfdb(data = wfdb,record = new_name,record_dir = new_dir)
  write_annotation(data = wfdb$annotation,record = new_name,annotator = 'ann',record_dir = new_dir)
}

old_dir <- 'wfdb'
new_dir <- 'wfdb-deid'


if (!dir.exists("new_dir")) {
  dir.create("new_dir", recursive = TRUE)
}

start <- 3640
end <- nrow(ecg)
for (i in start:end) {
  
  if (ecg$null_ecg[i] == TRUE) {
    next
  }
  
  old_name <- ecg$ecg_name[i]
  new_name <- as.character(ecg$ecg_number[i])
  
  deid_wfdb(old_name = old_name,old_dir = old_dir,new_name = new_name,new_dir = new_dir)
  if (!i%% 1000) {
    print(paste(i,'of',nrow(ecg)))
  }
}



# Summary of patients / ECGs: ---------------------------------------------
library(arrow)
library(dplyr)

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")

options(tibble.print_max = Inf)
options(tibble.width = Inf)

demo <- read_parquet('../demographics-0.parquet')
ecg <- read_parquet('../ecg_table.parquet')
log <- read_parquet('../wes-mrn-combined.parquet')

ecg <- ecg %>% filter(ttn_present == TRUE)

length(unique(ecg$mrn))
length(unique(ecg %>% filter(ttn_mutant == TRUE) %>% pull(mrn)))
length(unique(ecg %>% filter(ttn_mutant == FALSE) %>% pull(mrn)))

nrow(ecg)
length(ecg %>% filter(ttn_mutant == TRUE) %>% pull(mrn))
length(ecg %>% filter(ttn_mutant == FALSE) %>% pull(mrn))


length(unique(ecg %>% filter()))

# Summary of patients:
# 555 patients were sequenced total across both cohorts (556 rows, 1 duplicate: record_id 146708, batch 1 DNA ID 690/763)
# 548 patients have ECGs on the cluster (14823 ECGs)
  # 14609 ECGs were non-null ECGs, and sampled at 500 Hz (548 patients)

  # 7 patients do NOT have ECGs on the cluster
    # 5 have TTN sequencing data, 0 have TTN mutations

# 5 sequenced patients do not have record_ids, and are not in the patient data sets on cluster
  # 4 of 5 have TTN sequencing data, 1 of 5 is a TTN mutant

# 292  patients have TTN sequencing data.  
  # 72   mutants, 220  wild type
  # 7482 ECGs. 5865 mutants, 1617 wild type

# Find pre-dx ECGs

# Decisions:
#   12 lead vs single lead
#   10 sec vs R to R vs P to P, time averaged etc 

# deid and scp ------------------------------------------------------------
library(arrow)
library(dplyr)

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")

options(tibble.print_max = Inf)
options(tibble.width = Inf)

demo <- read_parquet('../demo.parquet')
ecg <- read_parquet('../ecg_table.parquet')
log <- read_parquet('../wes-mrn-combined.parquet')

demo_deid <- demo %>% select(-mrn,-first_name,-last_name,-ssn,-census_tract)
ecg_deid <- ecg %>% select(-ecg_name,-muse_path,-mrn)
log_deid <- log %>% select(-mrn,-first_name,-last_name,-ssn,-census_tract)

write_parquet(demo_deid,'../deid_demo.parquet')
write_parquet(ecg_deid,'../deid_ecg_table.parquet')
write_parquet(log_deid,'../deid_wes-mrn-combined.parquet')

scp -o MACs=hmac-sha2-256 -r dseaney2@131.193.182.245:/home/dseaney2/cardio_darbar_chi_link/common/cohorts/wes-ml-ttn/deid_demo.parquet 'C:\Users\darre\OneDrive\Documents\UICOM Research\WES ML'
scp -o MACs=hmac-sha2-256 -r dseaney2@131.193.182.245:/home/dseaney2/cardio_darbar_chi_link/common/cohorts/wes-ml-ttn/deid_ecg_table.parquet 'C:\Users\darre\OneDrive\Documents\UICOM Research\WES ML'
scp -o MACs=hmac-sha2-256 -r dseaney2@131.193.182.245:/home/dseaney2/cardio_darbar_chi_link/common/cohorts/wes-ml-ttn/deid_wes-mrn-combined.parquet 'C:\Users\darre\OneDrive\Documents\UICOM Research\WES ML'