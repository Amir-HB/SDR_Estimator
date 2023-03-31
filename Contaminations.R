# Replaces the most correlated observations with the empricial cov matrix
# by vectors proportional to the least pricipal eigenvector
# Output: contaminated Gaussian data in form of a matrix
flipping_contam <- function(n, p, eps, Sigma) {
  # n: sample size
  # p: dimension of samples
  # eps: contamination fraction
  # Sigma: the cov matrix of the Gaussian distribution of the inliers
  mu <- rep(0, p)
  o <- as.numeric(floor(n * eps))
  inliers <- mvrnorm(n=n, mu=mu,Sigma=Sigma)
  
  if (o == 0) {
    return(inliers)
  }
  
  scov = cov(inliers) 
  
  eigs <- eigen(scov, symmetric=TRUE)
  v1 <- eigs$vectors[,p]
  
  cors = abs(inliers%*%v1/p - apply(inliers,1,mean)*mean(v1))/(apply(inliers,1,std)*std(v1))
  
  topmost_ind <- head(sort(cors, index.return=TRUE, decreasing=TRUE)$ix, o)
  inliers[topmost_ind,] = matrix(rep(sqrt(p) * v1, o), ncol=p, nrow=o, byrow=TRUE)
  X = rbind(inliers[-topmost_ind,],inliers[topmost_ind,])
  
  return(X)
}


# This generates outliers from Gaussian distributions with uniform r.v. as center
generate_X <- function(n, p, eps, Sigma, outlier_type="large") {
  # Generate a matrix of size p times n.
  mu <- rep(0, p)
  o <- as.numeric(floor(n * eps))
  inliers <- mvrnorm(n=n - o, mu=mu, Sigma=Sigma)
  
  if(o == 0) {
    return(inliers)
  }
  
  if (outlier_type == "large") {
    outliers <- matrix(0L, nrow=p, ncol=o)
    for (i in 1:o) {
      theta <- runif(n=p, min=0, max=3)
      outliers[1:p, i] <- mvrnorm(n=1, mu=mu, Sigma=Sigma) + theta
    }
    #print(norm2(rowMeans(outliers)))
    X <- rbind(inliers, t(outliers))
    return(X)
  } else if (outlier_type == "const") {
    outliers <- 0.5 * ones(p, o)
    X <- rbind(inliers, t(outliers))
    return(X)
  }
  
}

#Generates data from a mixture of two Gaussians: one for inliers and one for outliers
mixture <- function(n,p,eps,Sigma,shift){
  # shift is the distance between the centers of the two Gaussians
  mu <- rep(0, p)
  o <- as.numeric(floor(n * eps))
  inliers <- mvrnorm(n=n - o, mu=mu, Sigma=Sigma)
  if(o == 0) return(inliers)
  outliers = mvrnorm(n=o, mu=(mu+shift), Sigma=Sigma)
  return(rbind(inliers, outliers))
}

scattered_contam <- function(n,p,eps,Sigma,d){
	mu <- rep(0, p)
	o <- as.numeric(floor(n * eps))
	data <- mvrnorm(n=n - o, mu=mu, Sigma=Sigma)
	if(o == 0) return(data)
	n_d = floor(o/d)
	shift = 20*p
	for(i in 1:d){	
		#u = zeros(n_d,p)
		u = ones(n_d,p)
		u[,i] = ((2*(runif(n_d)>.5))-1)*shift
		#u[,i] = ones(n_d,1)*sqrt(p/d)
  		data = rbind(data,u)
  	}
  	return(data)
}


