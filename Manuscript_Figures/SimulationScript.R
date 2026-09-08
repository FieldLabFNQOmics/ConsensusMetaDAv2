suppressPackageStartupMessages({library(SparseDOSSA2); library(MASS)})

################ Stool template parameter ######

## ---- Parameters ------------------------------------------------------------
sample_sizes <- c(100, 200)
spike_props  <- c(0.05, 0.10, 0.25)
n_feature    <- 500
n_metadata   <- 2       # meta_1 binary (tested), meta_2 continuous (nuisance)
eff_lo       <- 2.5
eff_hi       <- 5.0
depth        <- 50000
nrep         <- 100

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) spike_props <- as.numeric(strsplit(args[1], ",")[[1]])

## ---- Helpers ---------------------------------------------------------------
make_metadata <- function(n_sample, n_metadata) {
  M <- MASS::mvrnorm(n_sample, rep(0, n_metadata), diag(n_metadata))
  nb <- floor(n_metadata / 2)
  M[, seq_len(nb)] <- as.numeric(M[, seq_len(nb)] > 0)
  dimnames(M) <- list(paste0("Sample", seq_len(n_sample)),
                      paste0("meta_", seq_len(n_metadata)))
  M
}

## All spikes forced onto meta_1 so the truth-set size is CONSTANT across
## replicates. Random assignment made it swing widely, which would confound
## the sweep with truth-set-size noise.
# old
# make_spike <- function(n_feature, perc, lo, hi) {
#   k <- round(n_feature * perc)
#   f <- sample(paste0("Feature", seq_len(n_feature)), k)
#   do.call(rbind, lapply(f, function(ft) {
#     data.frame(metadata_datum = 1L,
#                feature_spiked = ft,
#                associated_property = c("abundance", "prevalence"),
#                effect_size = runif(2, lo, hi),
#                stringsAsFactors = FALSE)
#   }))
# }

make_spike <- function(n_feature, perc, lo, hi, frac_pos = 0.5) {
  k <- round(n_feature * perc)
  f <- sample(paste0("Feature", seq_len(n_feature)), k)
  signs <- ifelse(runif(k) < frac_pos, 1, -1)          # per-feature direction
  do.call(rbind, Map(function(ft, s) {
    data.frame(metadata_datum = 1L,
               feature_spiked = ft,
               associated_property = c("abundance", "prevalence"),
               effect_size = s * runif(2, lo, hi),      # signed magnitude
               stringsAsFactors = FALSE)
  }, f, signs))
}

run_one <- function(seed, n_sample, perc, outdir) {
  set.seed(seed)
  meta <- make_metadata(n_sample, n_metadata)
  sp   <- make_spike(n_feature, perc, eff_lo, eff_hi)
  
  sim <- SparseDOSSA2(template = "Stool", n_sample = n_sample,
                      new_features = TRUE, n_feature = n_feature,
                      spike_metadata = sp, metadata_matrix = meta,
                      median_read_depth = depth, verbose = FALSE)
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cnt   <- sim$simulated_data
  truth <- unique(sp$feature_spiked)
  
  write.table(cnt,  file.path(outdir, "abundance.tsv"), sep = "\t",
              quote = FALSE, col.names = NA)
  write.table(meta, file.path(outdir, "metadata.tsv"),  sep = "\t",
              quote = FALSE, col.names = NA)
  writeLines(truth, file.path(outdir, "truth_features.txt"))
  write.table(sp, file.path(outdir, "truth_full.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  
  c(n_samp = ncol(cnt), n_feat = nrow(cnt),
    n_true = length(truth), sparsity = round(mean(cnt == 0), 3))
}

## ---- Generate --------------------------------------------------------------
rows <- list()

for (pp in spike_props) {
  for (ns in sample_sizes) {
    
    cell <- file.path("sim_data", sprintf("n%d", ns), sprintf("spike%.2f", pp))
    cat(sprintf("\n=== %s ===\n", cell))
    st <- Sys.time()
    
    res <- t(sapply(seq_len(nrep), function(i)
      run_one(seed   = ns * 10000 + round(pp * 100) * 10 + i,
              n_sample = ns, perc = pp,
              outdir = file.path(cell, paste0("Replicate_", i)))))
    
    print(res)
    cat(sprintf("  %.1f min\n",
                as.numeric(difftime(Sys.time(), st, units = "mins"))))
    
    rows[[length(rows)+1]] <- data.frame(
      n_sample = ns, spike = pp,
      n_true = mean(res[, "n_true"]),
      sd_true = sd(res[, "n_true"]),
      sparsity = round(mean(res[, "sparsity"]), 3))
  }
}

cat("\n===== SUMMARY =====\n")
print(do.call(rbind, rows), row.names = FALSE)
cat("\nCheck: n_true should be exactly 20 / 100 / 150 / 200 for\n",
    "spike 0.05 / 0.10 / 0.25, with sd_true = 0 everywhere.\n")

################# IBD template ######


## ---- Parameters ------------------------------------------------------------
sample_sizes <- c(100, 200)
spike_props  <- c(0.05, 0.10, 0.25, 0.50)
n_feature    <- 2000
n_metadata   <- 2       # meta_1 binary (tested), meta_2 continuous (nuisance)
eff_lo       <- 2.5
eff_hi       <- 5.0
depth        <- 50000
nrep         <- 100

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) spike_props <- as.numeric(strsplit(args[1], ",")[[1]])

## ---- Helpers ---------------------------------------------------------------
make_metadata <- function(n_sample, n_metadata) {
  M <- MASS::mvrnorm(n_sample, rep(0, n_metadata), diag(n_metadata))
  nb <- floor(n_metadata / 2)
  M[, seq_len(nb)] <- as.numeric(M[, seq_len(nb)] > 0)
  dimnames(M) <- list(paste0("Sample", seq_len(n_sample)),
                      paste0("meta_", seq_len(n_metadata)))
  M
}

## All spikes forced onto meta_1 so the truth-set size is CONSTANT across
## replicates. Random assignment made it swing widely, which would confound
## the sweep with truth-set-size noise.
# old
# make_spike <- function(n_feature, perc, lo, hi) {
#   k <- round(n_feature * perc)
#   f <- sample(paste0("Feature", seq_len(n_feature)), k)
#   do.call(rbind, lapply(f, function(ft) {
#     data.frame(metadata_datum = 1L,
#                feature_spiked = ft,
#                associated_property = c("abundance", "prevalence"),
#                effect_size = runif(2, lo, hi),
#                stringsAsFactors = FALSE)
#   }))
# }

make_spike <- function(n_feature, perc, lo, hi, frac_pos = 0.5) {
  k <- round(n_feature * perc)
  f <- sample(paste0("Feature", seq_len(n_feature)), k)
  signs <- ifelse(runif(k) < frac_pos, 1, -1)          # per-feature direction
  do.call(rbind, Map(function(ft, s) {
    data.frame(metadata_datum = 1L,
               feature_spiked = ft,
               associated_property = c("abundance", "prevalence"),
               effect_size = s * runif(2, lo, hi),      # signed magnitude
               stringsAsFactors = FALSE)
  }, f, signs))
}
run_one <- function(seed, n_sample, perc, outdir) {
  set.seed(seed)
  meta <- make_metadata(n_sample, n_metadata)
  sp   <- make_spike(n_feature, perc, eff_lo, eff_hi)
  
  sim <- SparseDOSSA2(template = "IBD", n_sample = n_sample,
                      new_features = TRUE, n_feature = n_feature,
                      spike_metadata = sp, metadata_matrix = meta,
                      median_read_depth = depth, verbose = FALSE)
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cnt   <- sim$simulated_data
  truth <- unique(sp$feature_spiked)
  
  write.table(cnt,  file.path(outdir, "abundance.tsv"), sep = "\t",
              quote = FALSE, col.names = NA)
  write.table(meta, file.path(outdir, "metadata.tsv"),  sep = "\t",
              quote = FALSE, col.names = NA)
  writeLines(truth, file.path(outdir, "truth_features.txt"))
  write.table(sp, file.path(outdir, "truth_full.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  
  c(n_samp = ncol(cnt), n_feat = nrow(cnt),
    n_true = length(truth), sparsity = round(mean(cnt == 0), 3))
}

## ---- Generate --------------------------------------------------------------
rows <- list()

for (pp in spike_props) {
  for (ns in sample_sizes) {
    
    cell <- file.path("sim_data", sprintf("n%d", ns), sprintf("spike%.2f", pp))
    cat(sprintf("\n=== %s ===\n", cell))
    st <- Sys.time()
    
    res <- t(sapply(seq_len(nrep), function(i)
      run_one(seed   = ns * 10000 + round(pp * 100) * 10 + i,
              n_sample = ns, perc = pp,
              outdir = file.path(cell, paste0("Replicate_", i)))))
    
    print(res)
    cat(sprintf("  %.1f min\n",
                as.numeric(difftime(Sys.time(), st, units = "mins"))))
    
    rows[[length(rows)+1]] <- data.frame(
      n_sample = ns, spike = pp,
      n_true = mean(res[, "n_true"]),
      sd_true = sd(res[, "n_true"]),
      sparsity = round(mean(res[, "sparsity"]), 3))
  }
}

cat("\n===== SUMMARY =====\n")
print(do.call(rbind, rows), row.names = FALSE)
