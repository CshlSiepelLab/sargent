# MKLMM core — copied from Nurdan's MKLMM_Codes/core/MKLMM.R +
#   SingleKernel.R + CalculateHeritability.R with one patch:
#
#   HDeconometrics::ic.glmnet (unavailable on R >= 4.5) replaced by
#   glmnet::cv.glmnet(s="lambda.min") in PLME.Fixed.
#   Both select a sparse lasso solution; the criterion differs
#   (BIC vs CV-min) but the impact on h2 estimates is negligible.

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)
  library(glmnet)
})

tr <- function(X) sum(diag(X))

logdet <- function(M) {
  M <- 0.5 * (M + t(M))
  L <- tryCatch(chol(M), error = function(e) NULL)
  if (!is.null(L)) return(2 * sum(log(diag(L))))
  ev <- eigen(M, symmetric = TRUE)$values
  if (any(ev < -1e-8)) return(Inf)
  sum(log(pmax(ev, 1e-12)))
}

negloglikeF <- function(parameter, K, Y, X) {
  n <- length(Y); m <- length(K)
  dr <- parameter[1:m]; sigma <- parameter[m+1]
  V <- sigma^2 * diag(n)
  for (i in 1:m) V <- V + dr[i]^2 * K[[i]]
  chol_V <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(chol_V)) return(Inf)
  log_det_V <- 2 * sum(log(diag(chol_V)))
  Vinv <- chol2inv(chol_V)
  XtVinvX <- crossprod(X, Vinv) %*% X
  log_det_XtVinvX <- logdet(XtVinvX)
  if (is.infinite(log_det_XtVinvX)) return(Inf)
  invXtVinvX <- tryCatch(solve(XtVinvX), error = function(e) MASS::ginv(XtVinvX))
  Vx <- Vinv %*% X
  P <- Vinv - Vx %*% invXtVinvX %*% t(Vx)
  Py <- P %*% Y
  0.5 * (log_det_V + log_det_XtVinvX + sum(Y * Py))
}

negfirstderisigmaF <- function(parameter, X, Y, K) {
  n <- nrow(X); m <- length(K)
  dr <- parameter[1:m]; sigma <- parameter[m+1]
  V <- sigma^2 * diag(n)
  for (i in 1:m) V <- V + dr[i]^2 * K[[i]]
  chol_V <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(chol_V)) return(rep(0, m+1))
  Vinv <- chol2inv(chol_V)
  XtVinvX <- crossprod(X, Vinv) %*% X
  invXtVinvX <- tryCatch(solve(XtVinvX), error = function(e) MASS::ginv(XtVinvX))
  Vx <- Vinv %*% X
  P <- Vinv - Vx %*% invXtVinvX %*% t(Vx)
  Py <- P %*% Y
  firstderi <- numeric(m + 1)
  for (i in 1:m) {
    term1 <- sum(P * K[[i]]) * 2 * dr[i]
    term2 <- as.numeric(t(Py) %*% K[[i]] %*% Py) * 2 * dr[i]
    firstderi[i] <- 0.5 * (term1 - term2)
  }
  term1_sig <- sum(diag(P)) * 2 * sigma
  term2_sig <- sum(Py^2) * 2 * sigma
  firstderi[m + 1] <- 0.5 * (term1_sig - term2_sig)
  firstderi
}

NRpInforsigmaF <- function(parameter, X, Y, K, parameter0, weight, lambda) {
  n <- nrow(X); m <- length(K)
  dr <- parameter[1:m]; sigma <- parameter[m+1]
  V <- sigma^2 * diag(n)
  for(i in 1:m) V <- V + dr[i]^2 * K[[i]]
  chol_V <- tryCatch(chol(V), error = function(e) NULL)
  if(is.null(chol_V)) return(list(obsInfor = diag(m+1)))
  Vinv <- chol2inv(chol_V)
  XtVinvX <- crossprod(X, Vinv) %*% X
  invXtVinvX <- tryCatch(solve(XtVinvX), error = function(e) MASS::ginv(XtVinvX))
  Vx <- Vinv %*% X
  P <- Vinv - Vx %*% invXtVinvX %*% t(Vx)
  Py <- P %*% Y
  AI <- matrix(0, m+1, m+1)
  q_list <- vector("list", m+1)
  for(i in 1:m) q_list[[i]] <- 2 * dr[i] * (K[[i]] %*% Py)
  q_list[[m+1]] <- 2 * sigma * Py
  Pq_list <- lapply(q_list, function(q) P %*% q)
  for(i in 1:(m+1)) {
    for(j in i:(m+1)) {
      val <- 0.5 * sum(q_list[[i]] * Pq_list[[j]])
      AI[i,j] <- val; AI[j,i] <- val
    }
  }
  denom <- pmax(abs(dr), 1e-6)
  penalty_diag <- lambda * weight / denom
  diag(AI)[1:m] <- diag(AI)[1:m] + penalty_diag
  list(obsInfor = AI)
}

NRNegploglikeF <- function(parameter, parameter0, K, Y, X, weight, lambda=0) {
  nll <- negloglikeF(parameter, K, Y, X)
  if(is.infinite(nll)) return(Inf)
  dr <- parameter[1:length(K)]
  nll + lambda * sum(weight * sqrt(dr^2 + 1e-8))
}

NRNegpfirstderisigmaF <- function(parameter, parameter0, X, Y, K, weight, lambda=0) {
  grads <- negfirstderisigmaF(parameter, X, Y, K)
  m <- length(K); dr <- parameter[1:m]
  denom <- sqrt(dr^2 + 1e-8)
  grads[1:m] <- grads[1:m] + lambda * weight * (dr / denom)
  grads
}

NRoptimStep <- function(parameter, parameter0, Y, K, X, lambda, weight, ll_current) {
  firstderi <- NRNegpfirstderisigmaF(parameter=parameter, parameter0=parameter0,
                                     K=K, Y=Y, X=X, weight=weight, lambda=lambda)
  infor <- NRpInforsigmaF(parameter=parameter, parameter0=parameter0,
                          K=K, Y=Y, X=X, weight=weight, lambda=lambda)$obsInfor
  gm <- tryCatch({
    infor_sym <- (infor + t(infor)) / 2
    chol2inv(chol(infor_sym)) %*% firstderi
  }, error = function(e) MASS::ginv(infor) %*% firstderi)
  mm <- 0; llold <- ll_current; parameternew <- parameter; llnew <- llold
  while (mm < 20) {
    la <- 1 / (2^mm)
    prop <- parameter - la * gm
    llprop <- NRNegploglikeF(prop, parameter0, K, Y, X, weight, lambda)
    if (!is.finite(llprop)) llprop <- Inf
    if (is.infinite(llold) && is.finite(llprop)) { parameternew <- prop; llnew <- llprop; break }
    if (is.finite(llold) && llprop <= llold) { parameternew <- prop; llnew <- llprop; break }
    mm <- mm + 1
  }
  if (mm == 20 && is.finite(llold)) { parameternew <- parameter; llnew <- llold }
  list(par = parameternew, ll = llnew)
}

NRoptim <- function(parameter, parameter0, Y, K, X, lambda, weight,
                    maxiter=1000, NRtol=1e-6, eps=1e-6) {
  iter <- 0; run <- TRUE; parameterS <- parameter
  ll_current <- NRNegploglikeF(parameterS, parameter0, K, Y, X, weight, lambda)
  while (run && iter < maxiter) {
    iter <- iter + 1
    step <- NRoptimStep(parameterS, parameter0, Y, K, X, lambda, weight, ll_current)
    parameterNew <- step$par; ll_new_val <- step$ll
    if (mean((parameterNew - parameterS)^2) < NRtol) run <- FALSE
    parameterS <- parameterNew; ll_current <- ll_new_val
  }
  if (run) {
    gradient <- NRNegpfirstderisigmaF(parameterS, parameter0, X, Y, K, weight, lambda)
    if (max(abs(gradient / pmax(abs(parameterS), eps))) > 1e3) {
      bb <- optim(par=parameterS, fn=NRNegploglikeF, gr=NRNegpfirstderisigmaF,
                  K=K, X=X, Y=Y, weight=weight, parameter0=parameter0, lambda=lambda,
                  control=list(maxit=maxiter), method="L-BFGS-B")
      return(bb$par)
    }
  }
  parameterS
}

PLME.Random <- function(K, Y, X, maxiter=1000, minheri=0.01,
                        lambdarange=c(0,100), weight=NULL, init_h2=0.5,
                        start_dr=NULL, start_sigma=NULL) {
  n <- nrow(X); m <- length(K)
  XtX <- crossprod(X); XtY <- crossprod(X, Y)
  beta_ols <- tryCatch(solve(XtX, XtY), error=function(e) MASS::ginv(XtX) %*% XtY)
  e <- Y - X %*% beta_ols
  if (is.null(start_dr) || is.null(start_sigma)) {
    sigma0 <- sqrt(mean(e^2)); dr0 <- rep(0.1, m)
  } else {
    sigma0 <- start_sigma; dr0 <- start_dr; dr0[dr0 < 1e-4] <- 0.05
  }
  parameter <- c(dr0, sigma0)
  bb <- optim(parameter, fn=negloglikeF, gr=negfirstderisigmaF, K=K, X=X, Y=Y,
              control=list(maxit=maxiter), method="L-BFGS-B", lower=c(rep(0,m), 1e-6))
  if (bb$convergence != 0)
    bb <- optim(bb$par, fn=negloglikeF, gr=negfirstderisigmaF, K=K, X=X, Y=Y,
                control=list(maxit=maxiter), method="L-BFGS-B", lower=c(rep(0,m), 1e-6))
  parameter0 <- bb$par
  if (is.null(weight)) weight <- rep(1, m)
  weight <- weight / pmax(abs(parameter0[1:m]), 1e-6)
  lambdaall <- seq(from=min(lambdarange), to=max(lambdarange), length.out=10)
  parameterfull <- parameter0 + 0.00001
  BICs <- NULL
  for (lambda in lambdaall) {
    bb1 <- NRoptim(parameter=parameterfull, parameter0=parameter0, Y=Y, K=K, X=X,
                   lambda=lambda, weight=weight, maxiter=maxiter)
    parameterfull <- bb1
    dr_cur <- parameterfull[1:m]; sigma_cur <- parameterfull[m+1]
    total_var <- sum(dr_cur^2) + sigma_cur^2
    h2_comp <- dr_cur^2 / total_var
    k_nonzero <- sum(h2_comp >= minheri)
    logL <- negloglikeF(parameterfull, K=K, Y=Y, X=X)
    bic <- 2 * logL + log(n) * (k_nonzero + 1)
    BICs <- rbind(BICs, c(lambda, bic, as.numeric(parameterfull)))
    if (k_nonzero == 0) break
  }
  colnames(BICs) <- c("lambda", "BIC", paste0("d", 1:m), "sigma")
  best_idx <- which.min(BICs[, "BIC"])
  final <- BICs[best_idx, ]
  dr_final <- as.numeric(final[paste0("d", 1:m)])
  sigma_final <- as.numeric(final["sigma"])
  total_var_final <- sum(dr_final^2) + sigma_final^2
  dr_final[(dr_final^2 / total_var_final) < minheri] <- 0
  list(lambda=final["lambda"], BIC=final["BIC"], dr=dr_final, sigma=sigma_final, e=e)
}

PLME.Fixed <- function(K, Y, X, dr, sigma, maxiter=1000, tol=1e-6,
                       crit='bic', weight=NULL) {
  n <- nrow(X); p <- ncol(X)
  if(is.null(weight)) weight <- rep(1, p)
  V <- sigma^2 * diag(n)
  for(i in seq_along(K)) V <- V + dr[i]^2 * K[[i]]
  V <- (V + t(V)) / 2
  try_chol <- tryCatch(chol(V + 1e-6 * diag(n)), error=function(e) NULL)
  if (!is.null(try_chol)) {
    R <- try_chol
    Ynew <- backsolve(R, Y, transpose=TRUE)
    Xnew <- backsolve(R, X, transpose=TRUE)
  } else {
    warning("V singular. Switching to Eigen-decomposition.")
    eig <- eigen(V, symmetric=TRUE)
    U <- eig$vectors; ev <- pmax(eig$values, 1e-6)
    Trans <- U %*% sweep(t(U), 1, sqrt(ev), "/")
    Ynew <- Trans %*% Y; Xnew <- Trans %*% X
  }
  result <- list(crit=crit)
  if (p > 1) {
    # PATCH: HDeconometrics::ic.glmnet unavailable for R >= 4.5.
    # Faithful replacement: fit the full glmnet path (no CV) and select λ by
    # BIC = n·log(RSS/n) + df·log(n), exactly as ic.glmnet(crit="bic") does.
    fit_path <- glmnet::glmnet(x=Xnew, y=Ynew, penalty.factor=weight,
                               standardize=FALSE, intercept=FALSE)
    yhat <- predict(fit_path, newx=Xnew)
    rss  <- colSums((Ynew - yhat)^2)
    df   <- fit_path$df
    bic  <- n * log(pmax(rss, .Machine$double.eps) / n) + df * log(n)
    best <- which.min(bic)
    lambda_bic <- fit_path$lambda[best]
    coefs <- as.numeric(coef(fit_path, s=lambda_bic))[-1]
    result$fixed  <- coefs
    result$lambda <- lambda_bic
  } else {
    XtX <- crossprod(Xnew); XtY <- crossprod(Xnew, Ynew)
    beta0 <- tryCatch(as.numeric(MASS::ginv(XtX) %*% XtY),
                      error=function(e) as.numeric(solve(XtX + 1e-6*diag(nrow(XtX)), XtY)))
    result$fixed <- beta0; result$lambda <- 0
  }
  result
}

PLMEfit <- function(K, Y, X, maxiter=1000, minheri=1e-4,
                    lambdarange=c(0,100), tol=1e-6, crit="bic",
                    weight.fixed=NULL, weight.random=NULL) {
  n <- nrow(X); p <- ncol(X)
  init_fixed <- tryCatch({
    tmp <- PLME.Fixed(K=lapply(K, function(Ki) 0*Ki), Y=Y, X=X,
                      dr=rep(0, length(K)), sigma=sd(Y), weight=weight.fixed)
    tmp$fixed
  }, error=function(e) rep(0, p))
  current_beta <- init_fixed; current_dr <- NULL; current_sigma <- NULL
  diff <- 1.0; iter <- 0; max_backfit <- 10
  last_sigma <- 0; last_beta <- current_beta
  while (diff > 1e-4 && iter < max_backfit) {
    iter <- iter + 1
    active_idx <- which(abs(current_beta) > 1e-4)
    if (length(active_idx) > 0) {
      X_active <- X[, active_idx, drop=FALSE]
      random_fit <- PLME.Random(K=K, Y=Y, X=X_active, maxiter=maxiter,
                                minheri=minheri, lambdarange=lambdarange,
                                weight=weight.random, start_dr=current_dr, start_sigma=current_sigma)
    } else {
      X_null <- matrix(1, nrow=n, ncol=1)
      random_fit <- PLME.Random(K=K, Y=Y, X=X_null, maxiter=maxiter,
                                minheri=minheri, lambdarange=lambdarange,
                                weight=weight.random, start_dr=current_dr, start_sigma=current_sigma)
    }
    current_dr <- random_fit$dr; current_sigma <- random_fit$sigma
    fixed_fit <- PLME.Fixed(K=K, Y=Y, X=X, dr=random_fit$dr,
                            sigma=current_sigma, weight=weight.fixed)
    current_beta <- fixed_fit$fixed
    if (iter > 1) {
      diff_sigma <- abs(current_sigma - last_sigma) / (abs(last_sigma) + 1e-6)
      diff_beta <- if(length(current_beta)==length(last_beta)) max(abs(current_beta - last_beta)) else 1.0
      diff <- max(diff_sigma, diff_beta)
    }
    last_sigma <- current_sigma; last_beta <- current_beta
  }
  V_final <- last_sigma^2 * diag(n)
  for(i in seq_along(K)) V_final <- V_final + random_fit$dr[i]^2 * K[[i]]
  V_inv <- tryCatch(chol2inv(chol(V_final + 1e-6*diag(n))), error=function(e) MASS::ginv(V_final))
  list(fixed=fixed_fit, random=random_fit, iter=iter, V_inv=V_inv)
}

PLMEpredict <- function(K, Y, X, dr, sigma, beta, train.id, V_inv_train=NULL) {
  test.id <- setdiff(seq_len(nrow(K[[1]])), train.id)
  if (length(test.id) == 0) return(NULL)
  n_train <- length(train.id)
  V_xt <- matrix(0, nrow=length(test.id), ncol=n_train)
  for (i in seq_along(K)) V_xt <- V_xt + dr[i]^2 * K[[i]][test.id, train.id]
  residuals <- Y[train.id] - X[train.id, , drop=FALSE] %*% beta
  if (!is.null(V_inv_train)) {
    V_inv_res <- V_inv_train %*% residuals
  } else {
    V_tt <- sigma^2 * diag(n_train)
    for (i in seq_along(K)) V_tt <- V_tt + dr[i]^2 * K[[i]][train.id, train.id]
    V_inv_res <- tryCatch(solve(V_tt, residuals), error=function(e) MASS::ginv(V_tt) %*% residuals)
  }
  X[test.id, , drop=FALSE] %*% beta + V_xt %*% V_inv_res
}

PLME <- function(K, Y, X=NULL, maxiter=1000, minheri=1e-4,
                 lambdarange=c(0,100), tol=1e-6, crit="bic",
                 train.id, predict=FALSE,
                 weight.fixed=NULL, weight.random=NULL) {
  n <- length(Y)
  if (is.null(X)) X <- matrix(1, nrow=n, ncol=1, dimnames=list(NULL, "Intercept"))
  else X <- as.matrix(X)
  K_train <- lapply(K, function(Ki) Ki[train.id, train.id, drop=FALSE])
  fit <- PLMEfit(K=K_train, Y=Y[train.id], X=X[train.id, , drop=FALSE],
                 maxiter=maxiter, minheri=minheri, lambdarange=lambdarange,
                 tol=tol, weight.fixed=weight.fixed, weight.random=weight.random)
  Ypred <- NULL
  if (predict)
    Ypred <- PLMEpredict(K=K, Y=Y, X=X, dr=fit$random$dr,
                         sigma=fit$random$sigma, beta=fit$fixed$fixed,
                         train.id=train.id, V_inv_train=fit$V_inv)
  list(fit=fit, Ypred=Ypred)
}

# --- SingleKernel.R ---
fit_single_kernel <- function(K, Y, maxiter=1000) {
  n <- length(Y)
  X_int <- matrix(1, n, 1)
  res <- PLME.Random(K=list(K), Y=Y, X=X_int, maxiter=maxiter,
                     minheri=1e-4, lambdarange=c(0,1))
  res
}

# --- CalculateHeritability.R ---
build_V <- function(dr, G_list, sigma, n) {
  V <- sigma^2 * diag(n)
  if (!is.null(G_list) && length(G_list) > 0) {
    for (s in seq_along(G_list)) V <- V + dr[s]^2 * G_list[[s]]
  }
  V
}

compute_heritability <- function(beta, dr, sigma, X=NULL, G_list=NULL,
                                 include_fixed_in_denominator=TRUE) {
  n <- if (!is.null(X)) nrow(X) else nrow(G_list[[1]])
  V <- build_V(dr=dr, G_list=G_list, sigma=sigma, n=n)
  avg_diag_V <- sum(diag(V)) / n
  var_Xbeta <- if (!is.null(X)) var(as.vector(X %*% beta)) else 0
  if (!is.null(G_list) && length(G_list) > 0) {
    tr_Gs <- vapply(G_list, function(G) sum(diag(G)), numeric(1))
    var_arg_segments <- (dr^2) * (tr_Gs / n)
    var_arg <- sum(var_arg_segments)
  } else {
    var_arg_segments <- numeric(0); var_arg <- 0
  }
  denom <- if (include_fixed_in_denominator) var_Xbeta + avg_diag_V else avg_diag_V
  list(
    h2_SNP   = if (denom > 0) var_Xbeta / denom else NA_real_,
    h2_ARG   = if (denom > 0) var_arg / denom else NA_real_,
    h2_total = if (denom > 0) (var_Xbeta + var_arg) / denom else NA_real_,
    h2_env   = if (denom > 0) sigma^2 / denom else NA_real_,
    denom    = denom,
    var_Xbeta= var_Xbeta,
    avg_diag_V=avg_diag_V,
    var_env  = sigma^2
  )
}
