# Filter ------------------------------------------------------------------
library(arrow)
library(dplyr)
library(EGM)
library(purrr)
library(lubridate)
library(stringr)

setwd("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/code")

options(tibble.print_max = Inf)
options(tibble.width = Inf)

demo <- read_parquet('../demo.parquet')
log <- read_parquet('../wes-mrn-combined.parquet')
ecg <- read_parquet('../ecg_table.parquet')
vep <- read_parquet('../vep-combined.parquet')

common <- fs::path('/mmfs1/projects/cardio_darbar_chi/common')
source(fs::path(common,'software/cluster_general_software/AF_DM_functions.R'))
source(fs::path(common,'software/cluster_general_software/annotator_prep_functions.R'))
source(fs::path(common,'software/cluster_general_software/pvc_functions_TEMP.R'))


# Tweak names
vep <- vep %>% 
  dplyr::filter(!grepl(x=sample_id,pattern='UIC_Cases-AC')) %>% # remove files with AC in the name
  mutate(sample_id = sub(".*UIC(\\d{4}).*", "\\1", sample_id)) %>%
  
  mutate( # mutate the 1st cohort sample_ids, which only contain numbers
    sample_id = if_else(
      grepl("^0*\\d+$", sample_id),   # only digits (with possible leading zeros)
      as.character(as.integer(sample_id)),  # strip leading zeros
      sample_id                       # otherwise leave unchanged
    )
  )

# Choose method to filter mutations as pathogenic:
choose_method <- 3 # 1 or 2. see below:

if (choose_method == 1) {
  # Method 1 (Original)
  ttn_mutant <-
    vep |>
    dplyr::filter(CANONICAL == "YES") |>
    dplyr::filter(
      !is.na(LoF) |
        IMPACT %in% c("HIGH","MODERATE") |
        str_detect(CLIN_SIG, "pathogenic") |
        str_detect(
          Consequence,
          "stop_gained|missense_variant|frameshift_variant|splice_donor_variant|splice_acceptor_variant|start_lost"
        ) |
        str_detect(SIFT, "deleterious") |
        str_detect(PolyPhen, "damaging")
    ) |>
    dplyr::filter(MAX_AF < 0.01) |>
    dplyr::filter(!str_detect(CLIN_SIG, "benign"))
  
} else if (choose_method == 2) {
  # Method 2 (Claude)
  truncating_csq <- c("stop_gained", "frameshift_variant",
                      "splice_acceptor_variant", "splice_donor_variant")
  
  af_threshold <- 1e-3   # rare; novel variants (NA in MAX_AF) are kept
  #   original: 1e-3
  
  ttn_mutant <- vep |>
    dplyr::filter(SYMBOL == "TTN") |>
    # Prefer the cardiac meta-transcript over CANONICAL == "YES".
    # Confirm the transcript ID in your VEP cache (commonly the TTN
    # meta-transcript, e.g. ENST00000589042 / LRG_391):
    # filter(Feature == "ENST00000589042") |>
    mutate(
      CLIN_SIG   = ifelse(is.na(CLIN_SIG), "", CLIN_SIG),
      truncating = str_detect(Consequence,
                              str_c(truncating_csq, collapse = "|")),
      rare       = is.na(MAX_AF) | MAX_AF < af_threshold
    ) |>
    dplyr::filter(truncating, rare, !str_detect(CLIN_SIG, "benign"))
  # Cardiac-exon restriction (once you have the PSI table from
  # cardiodb.org/titin): join PSI per exon, then keep PSI > 0.9 (hiPSI):
  # |> left_join(ttn_psi, by = "exon_id") |> filter(psi > 0.9)
} else if (choose_method == 3) { # middle ground in terms of strictness
  truncating_csq <- c(
    "stop_gained", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant"
  )
  
  ttn_mutant <- vep |>
    dplyr::filter(SYMBOL == "TTN") |>
    mutate(
      truncating = str_detect(Consequence,
                              str_c(truncating_csq, collapse="|")),
      missense   = str_detect(Consequence, "missense_variant"),
      rare_trunc = is.na(MAX_AF) | MAX_AF < 1e-4,
      rare_miss  = is.na(MAX_AF) | MAX_AF < 1e-5,
      damaging_miss = missense &
        str_detect(SIFT, "deleterious") &
        str_detect(PolyPhen, "damaging")
    ) |>
    dplyr::filter(
      # Keep truncating variants if rare
      (truncating & rare_trunc) |
        # Keep only rare AND damaging missense
        (damaging_miss & rare_miss) |
        # Keep ClinVar pathogenic/likely pathogenic
        str_detect(CLIN_SIG, "pathogenic")
    ) |>
    # Exclude benign ClinVar
    dplyr::filter(!str_detect(CLIN_SIG, "benign"))
  
}

# Combine VEP data into WES log and ECG log
vep_all_unique <- vep %>% distinct(sample_id, batch)
vep_mutant_unique <- ttn_mutant %>% distinct(sample_id, batch)


old_mutant_ids <- vep_mutant_unique %>%
  dplyr::filter(batch == 1) %>%
  pull(sample_id)
new_mutant_ids <- vep_mutant_unique %>%
  dplyr::filter(batch == 2) %>%
  pull(sample_id)
log <- log %>%
  mutate(
    ttn_mutant = case_when(
      cohort == "old" & study_id %in% old_mutant_ids ~ TRUE,
      cohort == "new" & Sample.ID %in% new_mutant_ids ~ TRUE,
      TRUE ~ FALSE
    )
  )


old_all_ids <- vep_all_unique %>%
  dplyr::filter(batch == 1) %>%
  pull(sample_id)
new_all_ids <- vep_all_unique %>%
  dplyr::filter(batch == 2) %>%
  pull(sample_id)
log <- log %>%
  mutate(
    ttn_present = case_when(
      cohort == "old" & study_id %in% old_all_ids ~ TRUE,
      cohort == "new" & Sample.ID %in% new_all_ids ~ TRUE,
      TRUE ~ FALSE
    )
  )

log <- log %>% select(ttn_present,ttn_mutant,pt_number)


ecg <- ecg %>% #select(-ttn_present,-ttn_mutant) %>% 
  left_join(log,by='pt_number')


# select ECGs -------------------------------------------------------------

# Choose method to select ECGs for patients

# Include only sinus ECGs, with TTN sequencing, and good ECGs
ecg <- ecg %>% 
  dplyr::filter(ecg_sinus) %>% # include only sinus ECGs (using muse readout)
  dplyr::filter(ttn_present) %>% # include only ECGs of patients with TTN sequencing available (mutant or wild type)
  dplyr::filter(null_ecg == 0) %>% dplyr::filter(freq_500 == 1) # include only ECGs which are not empty (signal = 0 for whole length), and sampled at 500 hz

# Pick earliest ECG for each patient
# ecg <- ecg %>%
#   group_by(pt_number) %>%
#   slice_min(ecg_date, with_ties = FALSE) %>%
#   ungroup()

# Pick up to the 5 earliest ECGs per patient
# ecg <- ecg %>%
#   group_by(pt_number) %>%
#   slice_min(ecg_date, n = 5, with_ties = FALSE) %>%
#   ungroup()

# Pick ECGs from patients that are at least one month apart, then select up to 5 
library(dplyr)
library(lubridate)

max_ecgs <- 5

pick_ecgs <- function(df) {
  # Ensure ecg_date is Date and sorted
  df <- df %>%
    mutate(ecg_date = as.Date(ecg_date)) %>%
    arrange(ecg_date)
  
  n <- nrow(df)
  
  # If 0 or 1 ECG, just return it
  if (n <= 1) return(df)
  
  keep <- logical(n)
  keep[1] <- TRUE
  last_date <- df$ecg_date[1]
  count <- 1
  
  for (i in 2:n) {
    
    this_date <- df$ecg_date[i]
    
    # safety: skip malformed values
    if (is.na(this_date)) next
    
    # check spacing (≥30 days)
    if (as.numeric(this_date - last_date) >= 30) {
      keep[i] <- TRUE
      last_date <- this_date
      count <- count + 1
      if (count == max_ecgs) break
    }
  }
  
  df[keep, ]
}

ecg <- ecg %>%
  group_by(pt_number) %>%
  group_modify(~ pick_ecgs(.x)) %>%
  ungroup()


# Pull and prep (ECG segments) ----------------------------------------------
library(dplyr)
library(rlang)
library(purrr)
library(tibble)

ecg_dir <- "/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes-ml-ttn/wfdb-deid"

# Parameters: _______________________________________________

set_onset <- 'p_onset'
set_offset <- 't_offset'

lead <- 'I'

max_waves_per_ecg <- 1 # any, or integer to cap

extend_window <- 10 # extends the onset/offset in either direction to increase the wave length. Unlike max_length.
# Ie to include pre/post p-wave electrical activity, or adjust for any annotations which are cut short 
# Units: indices 
max_length <- 500 # total length, achieved by padding zeros to beginning/end. Units: indices (500 Hz)
# Ie 120 -> 240 ms (p-wave), 240 indices for t-wave, 500
exclude_termini <- 50 # exclude onsets and offsets which occur in this window from the signal termini

normalize_to_wholeECG <- FALSE # single normalization for the whole 10 sec ECG
normalize_to_indvWave <- FALSE # normalize each wave to only itself, not to whole 10 sec ECG

# ___________________________________________________________


# Split into wave + point
onset_wave  <- strsplit(set_onset,  "_")[[1]][1]
onset_point <- strsplit(set_onset,  "_")[[1]][2]

offset_wave  <- strsplit(set_offset, "_")[[1]][1]
offset_point <- strsplit(set_offset, "_")[[1]][2]


# Helper functions for creating wave windows
split_wave_point <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  list(wave = parts[1], point = parts[2])
}
make_wave_windows <- function(ann_compact, set_onset, set_offset) {
  
  o  <- split_wave_point(set_onset)
  of <- split_wave_point(set_offset)
  
  ann_by_channel <- split(ann_compact, ann_compact$channel)
  
  map_dfr(ann_by_channel, function(df_ch) {
    
    # Extract onset and offset times
    onset_times <- df_ch %>%
      dplyr::filter(type == o$wave) %>%
      arrange(!!sym(o$point)) %>%
      pull(!!sym(o$point))
    
    offset_times <- df_ch %>%
      dplyr::filter(type == of$wave) %>%
      arrange(!!sym(of$point)) %>%
      pull(!!sym(of$point))
    
    # Greedy matching with onset-skipping rule
    i <- 1
    j <- 1
    res_on  <- integer()
    res_off <- integer()
    
    while (i <= length(onset_times) && j <= length(offset_times)) {
      
      current_onset <- onset_times[i]
      current_offset <- offset_times[j]
      
      # If offset is before onset → move to next offset
      if (current_offset <= current_onset) {
        j <- j + 1
        next
      }
      
      # RULE A: skip onset if another onset occurs before this offset
      if (i < length(onset_times) && onset_times[i + 1] < current_offset) {
        # Skip this onset
        i <- i + 1
        next
      }
      
      # Otherwise, pair them
      res_on  <- c(res_on,  current_onset)
      res_off <- c(res_off, current_offset)
      
      i <- i + 1
      j <- j + 1
    }
    
    tibble(
      channel = unique(df_ch$channel),
      beat    = seq_along(res_on),
      onset   = res_on,
      offset  = res_off,
      length  = res_off - res_on
    )
  })
}


# Helper function for padding:
center_pad <- function(x, target_len) {
  n <- length(x)
  if (n >= target_len) {
    # truncate if too long *** revisit this
    return(x[1:target_len])
  }
  
  total_pad <- target_len - n
  left_pad  <- floor(total_pad / 2)
  right_pad <- ceiling(total_pad / 2)
  
  c(rep(0, left_pad), x, rep(0, right_pad))
}

all_waves <- list()
all_labels <- c()
all_pt_numbers <- c()

n <- nrow(ecg)
for (i in seq_len(n)) {
  
  ecg_name <- as.character(ecg$ecg_number[i])
  
  # read ECG
  wfdb <- read_wfdb(record = ecg_name, record_dir = ecg_dir,annotator = 'ann')
  sig_length <- length(wfdb$signal[[1]])
  lead_order <- setdiff(names(wfdb$signal),'sample')
  internal_lead_number <- which(lead_order == lead) # must handle carefully due to variable AVF/AVL/AVR order across institutions 
  
  ann_compact <- ann_wfdb2compact(wfdb$annotation)
  
  # Filter for lead number
  ann_compact <- ann_compact %>% dplyr::filter(channel == internal_lead_number)
  
  # Filter if waves are close to the beginning/end
  ann_compact <- ann_compact %>% 
    dplyr::filter(!onset < exclude_termini) %>% 
    dplyr::filter(!offset > (exclude_termini + sig_length))
  
  
  # Skip sample if it doesn't have the waves of interest:
  if (nrow(ann_compact %>% dplyr::filter(type == onset_wave)) == 0) {
    print(paste0('Skipping sample ',i,'. No ',onset_wave,' wave.'))
    next
  } 
  if (nrow(ann_compact %>% dplyr::filter(type == offset_wave)) == 0) {
    print(paste0('Skipping sample ',i,'. No ',offset_wave,' wave.'))
    next
  }
  
  # Create dataframe for windows of interest
  windows <- make_wave_windows(
    ann_compact,
    set_onset  = set_onset,
    set_offset = set_offset
  )
  
  # if there are more waves in the ECG than the max cap, remove randomly
  if (is.numeric(max_waves_per_ecg) & nrow(windows) > max_waves_per_ecg) {
    windows <- windows %>% dplyr::slice_sample(n = max_waves_per_ecg)
  }
  
  # Remove sample if window is too long
  if (any(windows$length > max_length)) {
    number_too_long <- sum(windows$length > max_length)
    
    if (number_too_long == nrow(windows)) {
      print(paste('All samples in row',i,'are tooo long.'))
      next
    }
    
    print(paste('Row',i,'has',number_too_long,'samples too long.', 
                'Lengths',paste(windows$length[windows$length > max_length],collapse=','),
                'Removing them.'))
    windows <- windows %>% dplyr::filter(length <= max_length)
  }
  
  if (nrow(windows) == 0) {
    print(paste0('Skipping sample ',i,'. Not enough onset/offset combos'))
    next
  } 
  
  # Adjust window for 'extend_window' parameter:
  windows <- windows %>%
    mutate(
      onset  = pmax(1, onset  - extend_window),
      offset = pmin(sig_length, offset + extend_window),
      length = offset-onset
    )
  
  
  # extract signal
  sig <- wfdb$signal[[lead]]
  # sanity check
  if (is.null(sig)) {
    warning(sprintf("Missing lead for row %d (%s)", i, ecg_name))
    next
  }
  # filter
  filt <- ecg_filter(sig)
  
  # Normalize whole ECG
  if (normalize_to_wholeECG) {
    filt <- (filt - min(filt)) / (max(filt) - min(filt)) * 100
  }
  
  # enforce length 5000
  if (length(filt) != 5000) {
    warning(sprintf("ECG length != 5000 for row %d (%s)", i, ecg_name))
    next
  }
  
  # Normalize and pad segments
  waves <- lapply(1:nrow(windows), function(j) {
    seg <- sig[windows$onset[j]:windows$offset[j]]
    if (normalize_to_indvWave) {
      seg <- (seg - min(seg)) / (max(seg) - min(seg)) * 100
    }
    if (length(seg) > max_length) {
      print(paste('signal is too long for max_length!',
                  'Sample:',i,
                  'Length:', length(seg)))
    }
    center_pad(seg, max_length)
  })
  
  # plot_func2(ecg_filter(wfdb$signal$II),wfdb$annotation %>% dplyr::filter(channel == 2))
  # plot_func2(waves[[1]])
  # result <- rowMeans(do.call(cbind, x)) # mean wave
  
  
  
  # store in list
  all_waves <- c(all_waves, waves)
  all_labels <- c(all_labels, rep(as.numeric(ecg$ttn_mutant[i]), length(waves)))
  all_pt_numbers <- c(all_pt_numbers, rep(ecg$pt_number[i], length(waves)))
  
  if (!i%%20) {
    print(paste(i,'of',n))
  }
}

signal <- array(unlist(all_waves), dim = c(length(all_waves), max_length, length(lead)))
labels <- all_labels
pt_numbers <- all_pt_numbers


# Build model and train (ECG segments) **in progress** ------------------------------------
library(keras)
library(dplyr)

train_size <- 0.7

unique_ids <- unique(pt_numbers)

# Split into training, test and validation (based on patient number, grouping same pt ECGs together)
set.seed(123)
train_ids <- sample(unique_ids, size = train_size * length(unique_ids))
remaining <- setdiff(unique_ids, train_ids)
val_ids <- sample(remaining, size = 0.15 * length(unique_ids))
test_ids <- setdiff(remaining, val_ids)

train_idx <- which(pt_numbers %in% train_ids)
val_idx   <- which(pt_numbers %in% val_ids)
test_idx  <- which(pt_numbers %in% test_ids)

signal_train <- signal[train_idx,,,drop=FALSE]
labels_train <- labels[train_idx]

signal_val <- signal[val_idx,,,drop=FALSE]
labels_val <- labels[val_idx]

signal_test <- signal[test_idx,,,drop=FALSE]
labels_test <- labels[test_idx]


# Build model
input_shape <- c(max_length, length(lead))
model <- keras_model_sequential() %>%
  layer_conv_1d(filters = 32, kernel_size = 7, activation = "relu",
                input_shape = input_shape) %>%
  layer_batch_normalization() %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  
  layer_conv_1d(filters = 64, kernel_size = 5, activation = "relu") %>%
  layer_batch_normalization() %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  
  layer_conv_1d(filters = 128, kernel_size = 5, activation = "relu") %>%
  layer_global_average_pooling_1d() %>%
  
  layer_dense(units = 64, activation = "relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  optimizer = optimizer_adam(learning_rate = 1e-3),
  loss = "binary_crossentropy",
  metrics = c("accuracy", metric_auc())
)


# Train model
history <- model %>% fit(
  signal_train, labels_train,
  validation_data = list(signal_val, labels_val),
  epochs = 10,
  batch_size = 8 # 256
  # callbacks = list(
  #   callback_early_stopping(patience = 5, restore_best_weights = TRUE)
  # )
)

model %>% evaluate(signal_test, labels_test)


# testing 
library(dplyr)

df <- tibble(pt_number = pt_numbers) %>%
  mutate(
    set = case_when(
      pt_number %in% train_ids ~ "train",
      pt_number %in% val_ids   ~ "val",
      pt_number %in% test_ids  ~ "test",
      TRUE ~ NA_character_
    )
  )



# Build model and train (ECG segments) (embedding model) **in progress** ------------------------------------
library(keras)

# Embedding:
max_length <- dim(signal)[2]

embedding_model <- keras_model_sequential() %>%
  layer_conv_1d(filters = 32, kernel_size = 7, activation = "relu",
                input_shape = c(max_length, 1)) %>%
  layer_batch_normalization() %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  
  layer_conv_1d(filters = 64, kernel_size = 5, activation = "relu") %>%
  layer_batch_normalization() %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  
  layer_conv_1d(filters = 128, kernel_size = 5, activation = "relu") %>%
  layer_global_average_pooling_1d() %>%
  
  layer_dense(units = 64, activation = "relu", name = "embedding")

embedding_model %>% compile(
  optimizer = optimizer_adam(1e-3),
  loss = "mse"   # dummy loss; we are not training this model directly
)

# Build model:
input <- layer_input(shape = c(max_length, 1))

embed <- input %>% embedding_model()

output <- embed %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dropout(0.3) %>%
  layer_dense(units = 1, activation = "sigmoid")

full_model <- keras_model(input, output)

full_model %>% compile(
  optimizer = optimizer_adam(1e-3),
  loss = "binary_crossentropy",
  metrics = c("accuracy", metric_auc())
)

# Split
unique_ids <- unique(pt_numbers)

set.seed(123)
train_ids <- sample(unique_ids, size = 0.7 * length(unique_ids))
remaining <- setdiff(unique_ids, train_ids)
val_ids <- sample(remaining, size = 0.15 * length(unique_ids))
test_ids <- setdiff(remaining, val_ids)

train_idx <- which(pt_numbers %in% train_ids)
val_idx   <- which(pt_numbers %in% val_ids)
test_idx  <- which(pt_numbers %in% test_ids)

X_train <- signal[train_idx,,,drop=FALSE]
y_train <- labels[train_idx]

X_val <- signal[val_idx,,,drop=FALSE]
y_val <- labels[val_idx]

# Train:
history <- full_model %>% fit(
  X_train, y_train,
  validation_data = list(X_val, y_val),
  epochs = 30,
  batch_size = 64,
  callbacks = list(
    callback_early_stopping(patience = 5, restore_best_weights = TRUE)
  )
)

# Extract embeddings
embeddings <- embedding_model %>% predict(signal)

# Aggregate embeddings
library(dplyr)

df_embed <- data.frame(
  patient_id = pt_numbers,
  label = labels,
  embeddings
)

patient_embeddings <- df_embed %>%
  group_by(patient_id) %>%
  summarise(
    label = first(label),
    across(starts_with("X"), mean)   # average embedding
  )

# Train a patient level classifier
X_pat <- as.matrix(patient_embeddings %>% select(starts_with("X")))
y_pat <- patient_embeddings$label

# simple classifier
clf <- keras_model_sequential() %>%
  layer_dense(32, activation = "relu", input_shape = ncol(X_pat)) %>%
  layer_dropout(0.3) %>%
  layer_dense(1, activation = "sigmoid")

clf %>% compile(
  optimizer = optimizer_adam(1e-3),
  loss = "binary_crossentropy",
  metrics = c("accuracy", metric_auc())
)

clf %>% fit(
  X_pat, y_pat,
  epochs = 50,
  batch_size = 16,
  validation_split = 0.2,
  callbacks = list(
    callback_early_stopping(patience = 5, restore_best_weights = TRUE)
  )
)

# Evaluate test set
test_pat <- patient_embeddings %>% dplyr::filter(patient_id %in% test_ids)

X_test_pat <- as.matrix(test_pat %>% select(starts_with("X")))
y_test_pat <- test_pat$label

clf %>% evaluate(X_test_pat, y_test_pat)

# Build model and train (ECG segments) (patient level model) ------------------------------------
library(dplyr)
library(purrr)
library(keras)

patients <- unique(pt_numbers)

patient_bags <- map(patients, function(pid) {
  idx <- which(pt_numbers == pid)
  list(
    X = signal[idx,,,drop=FALSE],
    y = labels[idx][1]   # same label for all beats
  )
})
names(patient_bags) <- patients

patient_generator <- function(patient_ids, batch_size = 8) {
  
  function() {
    
    # sample patients for this batch
    batch <- sample(patient_ids, batch_size, replace = TRUE)
    
    X_list <- list()
    y_list <- c()
    
    for (pid in batch) {
      bag <- patient_bags[[as.character(pid)]]
      X_list[[length(X_list) + 1]] <- bag$X
      y_list <- c(y_list, bag$y)
    }
    
    # determine max beats in this batch
    max_beats <- max(sapply(X_list, function(x) dim(x)[1]))
    
    # allocate padded batch array
    X_batch <- array(0, dim = c(batch_size, max_beats, max_length, 1))
    
    # fill in each patient's beats
    for (i in seq_along(X_list)) {
      beats <- X_list[[i]]
      n_beats <- dim(beats)[1]
      X_batch[i, 1:n_beats, , ] <- beats
    }
    
    # return batch (NO yield)
    list(X_batch, y_list)
  }
}


# Build model:
max_length <- dim(signal)[2]

pwave_encoder <- keras_model_sequential() %>%
  layer_conv_1d(32, 7, activation="relu", input_shape=c(max_length,1)) %>%
  layer_max_pooling_1d(2) %>%
  layer_conv_1d(64, 5, activation="relu") %>%
  layer_max_pooling_1d(2) %>%
  layer_conv_1d(128, 5, activation="relu") %>%
  layer_global_average_pooling_1d() %>%
  layer_dense(64, activation="relu")

library(keras)

max_beats <- max(sapply(patient_bags, function(x) dim(x$X)[1]))

input <- layer_input(shape = c(max_beats, max_length, 1))

# TimeDistributed applies the CNN to each beat
embeddings <- input %>%
  time_distributed(pwave_encoder)


patient_embedding <- embeddings %>%
  layer_lambda(function(x) k_mean(x, axis = 2)) 



output <- patient_embedding %>%
  layer_dense(32, activation="relu") %>%
  layer_dropout(0.3) %>%
  layer_dense(1, activation="sigmoid")

mil_model <- keras_model(input, output)

mil_model %>% compile(
  optimizer = optimizer_adam(1e-3),
  loss = "binary_crossentropy",
  metrics = c("accuracy", metric_auc())
)


# Train model
set.seed(123)
train_pat <- sample(patients, 0.7 * length(patients))
remaining <- setdiff(patients, train_pat)
val_pat <- sample(remaining, 0.15 * length(patients))
test_pat <- setdiff(remaining, val_pat)


train_gen <- patient_generator(train_pat, batch_size = 8)
val_gen <- patient_generator(val_pat, batch_size = 8)

mil_model %>% fit(
  train_gen,
  steps_per_epoch = length(train_pat) / 8,
  validation_data = val_gen,
  validation_steps = length(val_pat) / 8,
  epochs = 30,
  callbacks = list(
    callback_early_stopping(patience = 5, restore_best_weights = TRUE)
  )
)

test_gen <- patient_generator(test_pat, batch_size = length(test_pat))

mil_model %>% evaluate(test_gen)
