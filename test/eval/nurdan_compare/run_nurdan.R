#!/usr/bin/env Rscript
# run_nurdan.R — run MKLMM (ARG-LMM) on one Nurdan simulation case
#
# Usage:
#   Rscript run_nurdan.R <sim_dir> <trait> <pheno_col> <out_json>
#
# sim_dir   : path to Nurdan simulation directory (contains data/, pheno/, kernels/)
# trait     : e.g. vc_balanced, vc_highfix_lowrand, vc_lowfix_highrand
# pheno_col : replicate index 1..5 (matches truth.csv 'phenotype' column)
# out_json  : path to write JSON results

suppressPackageStartupMessages({ library(jsonlite) })

this_dir <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])))
source(file.path(this_dir, "mklmm_core.R"))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) stop("Usage: Rscript run_nurdan.R <sim_dir> <trait> <pheno_col> <out_json>")

sim_dir   <- args[1]
trait     <- args[2]
pheno_col <- as.integer(args[3])
out_json  <- args[4]

cat(sprintf("=== MKLMM (ARG-LMM) on %s / %s rep%d ===\n", sim_dir, trait, pheno_col))

t_start <- proc.time()["elapsed"]

# ---------- load genotypes ----------
geno_file <- file.path(sim_dir, "data", "all.genotype.csv")
cat("Loading genotypes...\n")
genotype <- read.csv(geno_file)
genome   <- as.matrix(genotype[, 7:ncol(genotype)])
genome   <- genome[, !duplicated(t(genome)), drop=FALSE]
n_indiv  <- nrow(genome)
n_snps   <- ncol(genome)
cat(sprintf("  %d individuals x %d SNPs\n", n_indiv, n_snps))

# ---------- standardize genotypes (all individuals) ----------
sd_all <- apply(genome, 2, sd)
X_all  <- genome[, sd_all > 0, drop=FALSE]
X_all  <- scale(X_all, center=colMeans(X_all), scale=apply(X_all,2,sd))

# ---------- load phenotype ----------
pheno_file <- file.path(sim_dir, "pheno", paste0("all.", trait, ".pheno"))
pheno <- read.table(pheno_file)
Y_all <- as.numeric(pheno[, pheno_col + 2])  # cols 1-2 are IDs
Y_all <- Y_all - mean(Y_all)

# ---------- load truth ----------
truth_file <- file.path(sim_dir, "pheno", paste0("all.", trait, ".truth.csv"))
truth <- read.csv(truth_file)
truth_row <- truth[truth$phenotype == pheno_col, ]

# ---------- load kernels ----------
cat("Loading ARG kernels...\n")
rdata_file <- file.path(sim_dir, "kernels", "Cmatrix_unique.RData")
load(rdata_file)   # loads Cmatrix_full_list
n_kernels <- length(Cmatrix_full_list)
cat(sprintf("  %d unique kernels\n", n_kernels))

# normalize: tr(K)/n = 1
K_all <- lapply(Cmatrix_full_list, function(Ki) {
  Ki <- (Ki + t(Ki)) / 2
  Ki / mean(diag(Ki))
})

# ---------- fit ARG-LMM ----------
cat("Fitting ARG-LMM (PLMEfit)...\n")
t_fit_start <- proc.time()["elapsed"]
fit <- PLMEfit(K=K_all, Y=Y_all, X=X_all, lambdarange=c(0,1), minheri=1e-4)
t_fit <- proc.time()["elapsed"] - t_fit_start

# ---------- heritability ----------
her <- compute_heritability(
  beta=fit$fixed$fixed, dr=fit$random$dr, sigma=fit$random$sigma,
  X=X_all, G_list=K_all
)

t_total <- proc.time()["elapsed"] - t_start

cat(sprintf("\n=== Results ===\n"))
cat(sprintf("  h2_SNP   = %.4f  (truth %.4f)\n", her$h2_SNP,   truth_row$h2_fix))
cat(sprintf("  h2_ARG   = %.4f  (truth %.4f)\n", her$h2_ARG,   truth_row$h2_rand))
cat(sprintf("  h2_total = %.4f  (truth %.4f)\n", her$h2_total, truth_row$h2_total))
cat(sprintf("  sigma2   = %.4f  (truth %.4f)\n", fit$random$sigma^2, truth_row$sigma2))
cat(sprintf("  n_active_kernels = %d  (truth %d)\n",
            sum(fit$random$dr > 0), truth_row$n_active_kernels))
cat(sprintf("  wall_fit = %.1f s  wall_total = %.1f s\n", t_fit, t_total))

# ---------- write JSON ----------
result <- list(
  method         = "MKLMM_ARG-LMM",
  sim_dir        = sim_dir,
  trait          = trait,
  pheno_col      = pheno_col,
  n_indiv        = n_indiv,
  n_snps         = n_snps,
  n_kernels      = n_kernels,
  h2_snp_est     = her$h2_SNP,
  h2_arg_est     = her$h2_ARG,
  h2_total_est   = her$h2_total,
  sigma2_est     = fit$random$sigma^2,
  n_active_kernels_est   = sum(fit$random$dr > 0),
  h2_snp_truth   = truth_row$h2_fix,
  h2_arg_truth   = truth_row$h2_rand,
  h2_total_truth = truth_row$h2_total,
  sigma2_truth   = truth_row$sigma2,
  n_active_kernels_truth = truth_row$n_active_kernels,
  wall_fit_s     = t_fit,
  wall_total_s   = t_total,
  n_backfit_iters= fit$iter
)
write(toJSON(result, auto_unbox=TRUE, digits=6), out_json)
cat(sprintf("\nWrote results to %s\n", out_json))
