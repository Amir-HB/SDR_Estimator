# Import libraries for Iterative reweighting
library(Rmosek)
#mosek_attachbuilder("~/mosek/9.3/tools/platform/linux64x86/bin", pos=2L, name="Rmosek:builder", warn.conflicts = TRUE)
#install.rmosek()
library(gsubfn)
library(ggplot2)
#library(rlist)
library(reticulate)
suppressWarnings(library(CVXR, warn.conflicts=FALSE))
installed_solvers() # Check if Mosek is installed

updateWeights <- function(eps, M, Sigma) {
  p <- NCOL(Sigma)
  n <- dim(M)[1] 
  w <- Variable(n, nonneg=TRUE)
  
  M_w <- t(M) %*% w
  M_w <- reshape_expr(M_w, c(p, p))
  
  objective <- Minimize(pos(lambda_max(M_w)))
  constr1 <- list(max_entries(w) <= 1./(n - n * eps))
  constr2 <- list(sum_entries(w) == 1.)
  ptm <- proc.time()
  prob <- Problem(objective, c(constr1, constr2))
  result <- solve(prob, solver="MOSEK", warm_start=TRUE)
  # feastol=sqrt(eps)/10)
  #print(c("Optimized in ", proc.time() - ptm))
  
  #print(result$status)
  #print(result$value)
  return(list(result$getValue(w), result$value))
}

mu_IR <- function (X, eps, Sigma) {
  # X is a matrix of size n times p. 
  # eps is a contamination parameter
  # Sigma is a covariance structure of each X_i. 
  
  n <- NROW(X)
  p <- NCOL(X)
  history <- list()
  
  gm <- geo_median(X)
  mu_0 <- gm$p
  
  if(eps > 0.25 || eps == 0) {
    K <- 10
  } else {
    # K <- 2 + floor(2 * log(18 * sqrt(p + 1), base=exp(1)) / log(1/ (eps * 3), base=exp(1)))
    # nom <- 2 * log(6 * sqrt(p + 1)) + 2 * log(1 / eps)
    nom <- log(4 * tr(Sigma) / norm(Sigma, type="2")) - 2 * log(eps * (1 - 2 * eps))
    denom <- 2 * log(1 - 2 * eps) - log(eps) - log(1 - eps)
    K <- floor(nom / denom)
  }
  
  mu <- mu_0
  
  weight_iters <- lapply(1:K, matrix, data=NA, nrow=1, ncol=n)
  # M <- lapply(1:n, matrix, data=NA, nrow=p, ncol=p)
  
  p2 <- p * p
  nas <- rep(NA, n*p2)
  M <- array(nas, c(n, p2))
  
  for(k in 1:K) {
    print(c("Iteration: ", k))
    ptm <- proc.time()
    for(i in 1:n) {
      M[i, ] <- as.vector(outer(X[i, ] - mu, X[i, ] - mu))
      # M[[i]] <- outer(X[i, ] - mu, X[i, ] - mu)
    }
    #print(proc.time() - ptm)
    
    val <- 0
    list[w, val] <- updateWeights(eps, M, Sigma)
    history <- c(history, val)
    
    #print(c("Current weight: ", w))
    # change to weighted mean
    #write.csv(file=paste("weights/iter", k), w)
    weight_iters[[k]] <- w
    mu <- rep(0, p)
    for(i in 1:n){
      mu <- mu + w[i] * X[i, ]
    }
    # print(c("Current mu: ", mu))
  }
  
  # write.csv(weight_iters, "data/weights_30.csv", row.names=FALSE)
  #print(history)
  
  return(mu)
}

