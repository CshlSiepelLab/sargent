#!/usr/bin/env Rscript
# run_nurdan_fig2.R — Figure 2 (missing-data) 6-model MKLMM + ARG-LMM Stage 2
#
# Usage:
#   Rscript run_nurdan_fig2.R <sim_dir> <pheno_col> <ds_factor> <seed> <out_json>
#
# sim_dir   : path to Figure 2 simulation directory
# pheno_col : column index in all.pheno1.pheno (3=Random, 5=Clustered)
# ds_factor : SNP downsampling factor (0=all SNPs, 1/2/4/10=fractions)
# seed      : integer random seed for EUR train/test split
# out_json  : output JSON path
#
# NOTE: Kernel pre-computation time (Python/tskit -> .npy -> .RData) is NOT
#       included in any model's wall_fit_s.  SARGENT includes equivalent
#       tree-traversal cost within its own wall time.

suppressPackageStartupMessages({
  library(MASS); library(Matrix); library(glmnet); library(jsonlite)
})

this_dir <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])))
source(file.path(this_dir, "mklmm_core.R"))
source(file.path(Sys.getenv("HOME"), "test-sargent/nurdan-all/NewApproach/arglmm.R"),
       local=TRUE)  # arg_lmm_stage2

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) stop("Usage: run_nurdan_fig2.R <sim_dir> <pheno_col> <ds_factor> <seed> <out_json>")

sim_dir   <- args[1]
pheno_col <- as.integer(args[2])
ds_factor <- as.integer(args[3])
seed      <- as.integer(args[4])
out_json  <- args[5]

cat(sprintf("=== MKLMM 6-model + ARG-LMM-Stage2 on %s ===\n", basename(sim_dir)))
cat(sprintf("    pheno_col=%d  ds_factor=%d  seed=%d\n", pheno_col, ds_factor, seed))

t0_total <- proc.time()["elapsed"]

# ---------- population labels ----------
labels <- read.table(file.path(sim_dir, "data", "population_labels.txt"), header=TRUE)
labels <- labels[seq(1, nrow(labels), 2), ]   # one row per diploid individual
N       <- nrow(labels)
eur_idx <- which(labels$population == "EUR")
afr_idx <- which(labels$population == "AFR")
cat(sprintf("  %d individuals: %d EUR, %d AFR\n", N, length(eur_idx), length(afr_idx)))

# ---------- train/test split (EUR only) ----------
set.seed(seed)
n_train  <- round(length(eur_idx) * 0.8)
train_id <- sort(sample(eur_idx, n_train))
test_id  <- sort(setdiff(eur_idx, train_id))
cat(sprintf("  train=%d EUR  test=%d EUR\n", length(train_id), length(test_id)))

# ---------- load phenotypes ----------
pheno_file <- file.path(sim_dir, "pheno", "all.pheno1.pheno")
pheno <- read.table(pheno_file)
Y_full  <- as.numeric(pheno[, pheno_col])
Y_mu    <- mean(Y_full[train_id])
Y_all   <- Y_full - Y_mu           # centered on train mean
Y_train <- Y_all[train_id]
Y_test  <- Y_all[test_id]

# ---------- load genotypes ----------
if (ds_factor == 0) {
  geno_file <- file.path(sim_dir, "data", "all.genotype.csv")
} else {
  geno_file <- file.path(sim_dir, "data", sprintf("all.%d.down.genotype.csv", ds_factor))
}
if (!file.exists(geno_file)) stop(paste("Genotype file not found:", geno_file))
cat(sprintf("  Loading genotypes: %s\n", basename(geno_file)))

genotype <- read.csv(geno_file)
genome   <- as.matrix(genotype[, 3:ncol(genotype)])
genome   <- genome[, !duplicated(t(genome)), drop=FALSE]

# standardize on train
sd_tr   <- apply(genome[train_id, , drop=FALSE], 2, sd)
keep    <- sd_tr > 0
genome  <- genome[, keep, drop=FALSE]
mu_X    <- colMeans(genome[train_id, , drop=FALSE])
sd_X    <- apply(genome[train_id, , drop=FALSE], 2, sd)
X_all   <- scale(genome, center=mu_X, scale=sd_X)
X_train <- X_all[train_id, , drop=FALSE]
X_test  <- X_all[test_id,  , drop=FALSE]
n_snps  <- ncol(X_all)
cat(sprintf("  %d individuals x %d SNPs (after dedup + sd-filter)\n", N, n_snps))

# ---------- PCA ----------
d_pca <- min(25, ncol(X_train) - 1)

# For heritability (all N): fit PCA on all individuals
pca_all <- prcomp(X_all, center=FALSE, scale.=FALSE)
P_all   <- pca_all$x[, seq_len(d_pca), drop=FALSE]
colnames(P_all) <- paste0("PC", seq_len(d_pca))

# For prediction (EUR train/test): fit PCA on train, project test
pca_tr  <- prcomp(X_train, center=FALSE, scale.=FALSE)
P_train <- pca_tr$x[, seq_len(d_pca), drop=FALSE]
P_test  <- X_test %*% pca_tr$rotation[, seq_len(d_pca), drop=FALSE]
colnames(P_train) <- colnames(P_test) <- paste0("PC", seq_len(d_pca))

# ---------- load ARG kernels ----------
cat("  Loading kernels...\n")
t_kern_load <- system.time({
  load(file.path(sim_dir, "kernels", "Cmatrix_unique.RData"))
})["elapsed"]
n_kernels <- length(Cmatrix_full_list)
cat(sprintf("  %d unique kernels (load: %.1fs)\n", n_kernels, t_kern_load))

# normalize tr(K)/n = 1 using all-individual diagonal
K_all <- lapply(Cmatrix_full_list, function(Ki) {
  Ki <- (Ki + t(Ki)) / 2
  Ki / mean(diag(Ki))
})
K_avg   <- list(Reduce("+", K_all) / n_kernels)
X_int   <- matrix(1, N, 1)                      # intercept only

# GRM from genotypes
GRM_all <- tcrossprod(X_all) / n_snps

# Helper ------------------------------------------------------------------
eval_pred <- function(ypred, ytrue) {
  if (is.null(ypred) || all(is.na(ypred))) return(list(r2=NA, cor=NA))
  ss_res <- sum((ytrue - ypred)^2)
  ss_tot <- sum((ytrue - mean(ytrue))^2)
  list(r2=1 - ss_res/ss_tot, cor=suppressWarnings(cor(ypred, ytrue)))
}

run_model <- function(label, expr) {
  cat(sprintf("    [%s]...\n", label))
  t1 <- proc.time()["elapsed"]
  res <- tryCatch(expr, error=function(e) { message("  ERROR: ", e$message); NULL })
  list(result=res, wall_s=proc.time()["elapsed"] - t1)
}

# =========================================================================
# HERITABILITY MODELS (all N individuals)
# =========================================================================
cat("\n--- Heritability (all individuals) ---\n")

h_arg <- run_model("ARG-LMM", {
  fit <- PLMEfit(K=K_all, Y=Y_all, X=X_all, lambdarange=c(0,1), minheri=1e-4)
  compute_heritability(beta=fit$fixed$fixed, dr=fit$random$dr,
                       sigma=fit$random$sigma, X=X_all, G_list=K_all)
})

h_pca <- run_model("ARG-LMM-PCA", {
  fit <- PLMEfit(K=K_all, Y=Y_all, X=P_all, lambdarange=c(0,1), minheri=1e-4)
  compute_heritability(beta=fit$fixed$fixed, dr=fit$random$dr,
                       sigma=fit$random$sigma, X=P_all, G_list=K_all)
})

h_jr <- run_model("Just-Random", {
  fit <- PLMEfit(K=K_all, Y=Y_all, X=X_int, lambdarange=c(0,10), minheri=1e-4)
  compute_heritability(beta=fit$fixed$fixed, dr=fit$random$dr,
                       sigma=fit$random$sigma, X=NULL, G_list=K_all)
})

h_grm_arg <- run_model("GRM-LMM(ARG)", {
  fit <- PLMEfit(K=K_avg, Y=Y_all, X=X_all, lambdarange=c(0,0.1), minheri=1e-4)
  compute_heritability(beta=fit$fixed$fixed, dr=fit$random$dr,
                       sigma=fit$random$sigma, X=X_all, G_list=K_avg)
})

h_grm_snp <- run_model("GRM-LMM(SNP)", {
  fit <- fit_single_kernel(K=GRM_all, Y=Y_all)
  compute_heritability(beta=rep(0, 1), dr=fit$dr, sigma=fit$sigma,
                       X=NULL, G_list=list(GRM_all))
})

h_jf <- run_model("Just-Fixed", {
  cv  <- glmnet::cv.glmnet(x=X_all, y=Y_all, alpha=1, intercept=FALSE,
                            standardize=FALSE, nfolds=10)
  beta_jf <- as.numeric(coef(cv, s="lambda.1se"))[-1]
  var_xb  <- var(as.numeric(X_all %*% beta_jf))
  list(h2_SNP=var_xb / (var_xb + var(Y_all - X_all %*% beta_jf)),
       h2_ARG=0, h2_total=NA, beta=beta_jf)
})

# =========================================================================
# PREDICTION MODELS (EUR 80/20 train/test)
# Reorder individuals to [train | test] so PLMEpredict returns exactly
# n_test predictions matching Y_test, following Test_Fig2_Prediction.R.
# =========================================================================
cat("\n--- Prediction (EUR 80/20) ---\n")

n_tr    <- length(train_id)
n_te    <- length(test_id)
ind_ord <- c(train_id, test_id)
tr_id   <- seq_len(n_tr)
te_id   <- seq(n_tr + 1, n_tr + n_te)

# Subsets reordered for prediction
K_pred    <- lapply(K_all, function(Ki) { Ki2 <- Ki[ind_ord, ind_ord]; (Ki2+t(Ki2))/2 })
K_avg_pred <- list(Reduce("+", K_pred) / n_kernels)
X_pred    <- X_all[ind_ord, , drop=FALSE]
Y_pred    <- Y_all[ind_ord]
P_pred    <- rbind(P_train, P_test)   # already [train | test] order
colnames(P_pred) <- paste0("PC", seq_len(d_pca))
GRM_pred  <- tcrossprod(X_pred) / n_snps

p_arg <- run_model("ARG-LMM pred", {
  fit <- PLME(K=K_pred, Y=Y_pred, X=X_pred, train.id=tr_id,
              lambdarange=c(0,1), minheri=1e-4, predict=TRUE)
  eval_pred(as.numeric(fit$Ypred), Y_test)
})

p_pca <- run_model("ARG-LMM-PCA pred", {
  fit <- PLME(K=K_pred, Y=Y_pred, X=P_pred, train.id=tr_id,
              lambdarange=c(0,1), minheri=1e-4, predict=TRUE)
  eval_pred(as.numeric(fit$Ypred), Y_test)
})

p_jr <- run_model("Just-Random pred", {
  K_tr  <- lapply(K_pred, function(Ki) Ki[tr_id, tr_id])
  fit_r <- PLME.Random(K=K_tr, Y=Y_pred[tr_id], X=matrix(1, n_tr, 1),
                       minheri=1e-4, lambdarange=c(0,10))
  V_xt  <- matrix(0, n_te, n_tr)
  V_tt  <- fit_r$sigma^2 * diag(n_tr)
  for (i in seq_along(K_pred)) {
    V_xt <- V_xt + fit_r$dr[i]^2 * K_pred[[i]][te_id, tr_id]
    V_tt <- V_tt + fit_r$dr[i]^2 * K_pred[[i]][tr_id, tr_id]
  }
  eval_pred(as.numeric(V_xt %*% MASS::ginv(V_tt) %*% Y_pred[tr_id]), Y_test)
})

p_grm_arg <- run_model("GRM-LMM(ARG) pred", {
  fit <- PLME(K=K_avg_pred, Y=Y_pred, X=X_pred, train.id=tr_id,
              lambdarange=c(0,0.1), minheri=1e-4, predict=TRUE)
  eval_pred(as.numeric(fit$Ypred), Y_test)
})

p_grm_snp <- run_model("GRM-LMM(SNP) pred", {
  GRM_tr  <- GRM_pred[tr_id, tr_id]
  fit_snp <- fit_single_kernel(K=GRM_tr, Y=Y_pred[tr_id])
  V_tt_g  <- fit_snp$sigma^2 * diag(n_tr) + fit_snp$dr^2 * GRM_tr
  V_xt_g  <- fit_snp$dr^2 * GRM_pred[te_id, tr_id]
  mu_tr   <- mean(Y_pred[tr_id])
  ypred   <- mu_tr + as.numeric(V_xt_g %*% MASS::ginv(V_tt_g) %*% (Y_pred[tr_id] - mu_tr))
  eval_pred(ypred, Y_test)
})

p_jf <- run_model("Just-Fixed pred", {
  cv      <- glmnet::cv.glmnet(x=X_pred[tr_id, , drop=FALSE], y=Y_pred[tr_id],
                                alpha=1, intercept=FALSE, standardize=FALSE, nfolds=10)
  beta_jf <- as.numeric(coef(cv, s="lambda.1se"))[-1]
  eval_pred(as.numeric(X_pred[te_id, , drop=FALSE] %*% beta_jf), Y_test)
})

# =========================================================================
# ARG-LMM Stage 2  (arg_lmm_stage2, NewApproach)
# =========================================================================
cat("\n--- ARG-LMM Stage 2 (BayesR + AIREML) ---\n")
K_train_list <- lapply(K_pred, function(Ki) Ki[tr_id, tr_id])
K_test_list  <- lapply(K_pred, function(Ki) Ki[te_id, tr_id])

p_stage2 <- run_model("ARG-LMM-Stage2", {
  fit2 <- arg_lmm_stage2(
    y_train=Y_pred[tr_id], X_train=X_pred[tr_id, , drop=FALSE],
    K_train=K_train_list,  K_test=K_test_list,
    X_test=X_pred[te_id, , drop=FALSE], y_test=Y_test,
    n_pc=10, n_iter=100,
    normalize=FALSE,  # already normalized above
    verbose=TRUE
  )
  list(
    r2_full = fit2$eval$r2s["full"],
    r2_snp  = fit2$eval$r2s["snp"],
    r2_blup = fit2$eval$r2s["blup"],
    r2_noV  = fit2$eval$r2s["noV"],
    cor_full= fit2$eval$cors["full"],
    cor_snp = fit2$eval$cors["snp"],
    sigma2  = fit2$sigma2,
    sum_w2  = sum(fit2$w2),
    h2_arg  = sum(fit2$w2) / (sum(fit2$w2) + fit2$sigma2),
    n_active= sum(apply(fit2$gamma, 1, function(g) g[1] < 0.5))
  )
})

t_total <- proc.time()["elapsed"] - t0_total
cat(sprintf("\nTotal wall time: %.1f s\n", t_total))

# =========================================================================
# Print summary
# =========================================================================
extract_h <- function(m) {
  r <- m$result
  if (is.null(r)) return(c(NA, NA, NA))
  c(r$h2_SNP, r$h2_ARG, r$h2_total)
}
extract_r2 <- function(m) if (is.null(m$result)) NA else m$result$r2

cat(sprintf("\n%s\n", strrep("=", 65)))
cat(sprintf("  %-20s  %7s  %7s  %7s  %7s  %6s\n",
            "Model", "h2_snp", "h2_arg", "h2_tot", "R2_pred", "time_s"))
cat(sprintf("%s\n", strrep("-", 65)))
models_h  <- list(h_arg, h_pca, h_jr, h_grm_arg, h_grm_snp, h_jf)
models_p  <- list(p_arg, p_pca, p_jr, p_grm_arg, p_grm_snp, p_jf)
mnames    <- c("ARG-LMM", "ARG-LMM-PCA", "Just-Random",
               "GRM-LMM(ARG)", "GRM-LMM(SNP)", "Just-Fixed")
for (i in seq_along(mnames)) {
  hv <- extract_h(models_h[[i]])
  r2 <- extract_r2(models_p[[i]])
  tw <- (models_h[[i]]$wall_s + models_p[[i]]$wall_s)
  cat(sprintf("  %-20s  %7.4f  %7.4f  %7.4f  %7.4f  %6.1f\n",
              mnames[i],
              ifelse(is.na(hv[1]), NA, hv[1]),
              ifelse(is.na(hv[2]), NA, hv[2]),
              ifelse(is.na(hv[3]), NA, hv[3]),
              ifelse(is.na(r2), NA, r2),
              tw))
}
# ARG-LMM Stage 2
if (!is.null(p_stage2$result)) {
  s2 <- p_stage2$result
  cat(sprintf("  %-20s  %7s  %7.4f  %7s  %7.4f  %6.1f\n",
              "ARG-LMM-Stage2", "BayesR",
              s2$h2_arg, "joint",
              ifelse(is.na(s2$r2_full), NA, s2$r2_full),
              p_stage2$wall_s))
}
cat(sprintf("%s\n", strrep("=", 65)))
cat(sprintf("  Kernel load: %.1f s  |  Total: %.1f s\n", t_kern_load, t_total))

# =========================================================================
# Write JSON
# =========================================================================
mklmm_herit <- setNames(lapply(seq_along(mnames), function(i) {
  hv <- extract_h(models_h[[i]])
  list(h2_snp=hv[1], h2_arg=hv[2], h2_total=hv[3], wall_s=models_h[[i]]$wall_s)
}), mnames)

mklmm_pred <- setNames(lapply(seq_along(mnames), function(i) {
  list(r2=extract_r2(models_p[[i]]),
       cor=if (!is.null(models_p[[i]]$result)) models_p[[i]]$result$cor else NA,
       wall_s=models_p[[i]]$wall_s)
}), mnames)

stage2_out <- if (!is.null(p_stage2$result)) {
  s2 <- p_stage2$result
  list(r2_full=unname(s2$r2_full), r2_snp=unname(s2$r2_snp),
       r2_blup=unname(s2$r2_blup), r2_noV=unname(s2$r2_noV),
       h2_arg=s2$h2_arg, sigma2=s2$sigma2, sum_w2=s2$sum_w2,
       n_active=s2$n_active, wall_s=p_stage2$wall_s)
} else NULL

result <- list(
  method="6model+stage2",
  sim_dir=sim_dir, pheno_col=pheno_col, ds_factor=ds_factor, seed=seed,
  n_indiv=N, n_snps=n_snps, n_kernels=n_kernels,
  n_train=length(train_id), n_test=length(test_id),
  kern_load_s=t_kern_load, wall_total_s=t_total,
  heritability=mklmm_herit,
  prediction=mklmm_pred,
  arglmm_stage2=stage2_out
)
write(toJSON(result, auto_unbox=TRUE, digits=6, null="null"), out_json)
cat(sprintf("\nWrote results to %s\n", out_json))
