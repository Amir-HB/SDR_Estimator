# Import libraries needed for SDR
library(Gmedian)
library(psych)
library(pracma)
library(Matrix)
library(MASS)

# This function determines the dimension reduction regime
# Our choice is to divide dimension by e at each iteration 
# The output is a vector indicating the dimension for each iteration
p_l <-function(p){
	# p is the initial dimension of data
  v = c()
	while(p>1){
		v = c(v,p)
		p = ceiling(p/exp(1))
	}
	v = c(v,p)
	return(v)
}

# SDR estimator as described in the paper
mu_SDR <- function(X,t,Sigma){
  # X is the design matrix (n*p) of the contaminated data
  # t is the threshold's value used for filtering data
  # The varaibles are named as in the algorithm presented in the paper
  n <- NROW(X) # Number of data points
  p <- NCOL(X) # Dimension of data points
  pl = p_l(p) # Dimension reduction regime
  L = length(pl) # Number of iterations
  V = diag(p) # New subspace with reduced dimension
  mu_hat = zeros(p,1) # Our final estimation of mean
  if (p>1){
    for(l in 1:(L-1)){
      cat(" **  Step ",  l, " : Dim = ",  pl[l], "\n")
      # Projection of the point over subspace V
      XV = X %*% V
      # Geometric median of the projected points
      bar_mu_V <- geo_median(XV,tol=1e00,maxiter=15)$p
      #bar_mu_V <- geo_median(XV)$p
      # Distance of the points from the geometric median
      Z = apply((t(XV)-bar_mu_V)^2,2,sum)
      # Remove points with large distance
      X_filtered = XV[Z<(t^2*pl[l]),]
      # Mean of the filtered projected points
      Xbar_filtered = apply(X_filtered,2,mean)
      
      # Covariance matrix of the filtered projected points
      S_filtered = cov(X_filtered);
      # Eigen values and vectors of the difference between
      # the emprical cov matrix and the theoretical cov matrix over subspace V
      eigs <- eigen(S_filtered - t(V) %*%Sigma%*% V, symmetric=TRUE)
      
      # Subspace spanned by the bottom pricipal components
      U = eigs$vectors[,(pl[l+1]+1):pl[l]]
      
      # Mean of the filtered projected points over subspace U
      mu_l = V %*% U %*% t(U) %*% t(t(Xbar_filtered))
      # Adding the mean estimation for subspace U to the final estimation
      mu_hat = mu_hat + mu_l
      # Subspace spanned by the top pricipal components
      U_ortho = eigs$vectors[,1:pl[l+1]]
      # The new subspace to project the points over and work on
      V = V %*% U_ortho
      #cat(" **    mu_hat = [", round(mu_hat[1:5],2) ,"]\n")
    }
  }
  # print(pl[L])
  # Projection of the points over subspace V
  XV = X %*% V
  # The last subspace of the dimension reduction is not necessarily of dimension one
  if(pl[L]>1)
  {
    bar_mu_V = geo_median(XV,tol = 1e-00,maxiter=15)$p
    Z = apply((t(XV)-bar_mu_V)^2,2,sum)
    XV_filtered = XV[Z<t^2*pl[L],]
    Xbar_filtered = apply(XV_filtered, 2, mean)
  }
  cat(" **  Step ",  L, " : Dim = ",  pl[L], "\n")
  if(pl[L]==1)
  {
    #bar_mu_V = median(XV)
    #Z = (XV-bar_mu_V)^2
    #XV_filtered = XV[Z<t^2*pl[L]]
    #Xbar_filtered = mean(XV_filtered)
    Xbar_filtered = median(XV)
  }
  mu_L = V%*%t(t(Xbar_filtered))
  mu_hat = mu_hat + mu_L
  #cat(" **    mu_hat[1:5] = [", round(mu_hat[1:5],2) ,"]\n")
  return(mu_hat)
}     
