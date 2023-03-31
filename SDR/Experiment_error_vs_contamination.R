source("SDR.R")
source("IRM.R")
source("Contaminations.R")
source("Experiment_materiels.R")

n = 100
p = 10
#n = 10*floor(p/eps^2)
eps = .2
Sigma = diag(p)
parameter = "eps" # choices: n, p, eps
range = c(seq(.04,.48,.06))
#range = c(seq(5,65,10))
#range = c(seq(50,300,50))
contamination = "flip" # existing contaminations: "mixt", "unif", "mixt"
estimators= c("SDR"
				#,"Iterative filtering"
				,"Geometric median"
				#,"Oracle"
				#,"Sample mean"
				,"Componentwise median"
				#,"Iterative reweighting"
				#,"Oracle"
				#,"Tukey's median"
				#,"New iterative reweighting"
				) 
#existing estimators: SDR", "Iterative filtering","Geometric median","Oracle","Sample mean", "Iterative reweighting","Tukey's median"

rep_num = 50
loss = variation(n,p,eps,Sigma,parameter,range,rep_num,contamination,estimators)

estimators2 = c(#"SDR"
				#,"Iterative filtering"
				#,"Geometric median"
				#"Oracle"
				#,"Sample mean"
				"Iterative reweighting"
				#"Tukey's median"
				) 
rep_num = 10
loss_2 = variation(n,p,eps,Sigma,parameter,range,rep_num,contamination,estimators2)
loss = rbind(loss,loss_2)
estimators = c(estimators,estimators2)

estimators2 = c(#"SDR"
				#,"Iterative filtering"
				#,"Geometric median"
				"Oracle"
				#,"Sample mean"
				#"Iterative reweighting"
				#"Tukey's median"
				) 
rep_num = 50
loss_2 = variation(n,p,eps,Sigma,parameter,range,rep_num,contamination,estimators2)
loss = rbind(loss,loss_2)
estimators = c(estimators,estimators2)



#write.table(rbind(range,loss), "err_eps_flip", row.names = FALSE, col.names = FALSE)
#Tabl = read.table("err_eps_tuk")
#M = as.matrix(Tabl)
#range = M[1,]
#loss = M[-1,]


rows <- nrow(loss)
# extracting odd rows
odd_cols <- seq_len(rows) %% 3
q1 <- loss[odd_cols == 1,]
q2 <- loss[odd_cols == 2,]
q3 <- loss[odd_cols == 0,]

rows <- nrow(q1) 
odd_cols <- seq_len(rows) %% 2
l_q1 <- matrix(q1[odd_cols == 1,],ncol=ncol(loss))
t_q1 <- matrix(q1[odd_cols == 0,],ncol=ncol(loss))

rows <- nrow(q2) 
odd_cols <- seq_len(rows) %% 2
l_q2 <- matrix(q2[odd_cols == 1,],ncol=ncol(loss))
t_q2 <- matrix(q2[odd_cols == 0,],ncol=ncol(loss))

rows <- nrow(q3) 
odd_cols <- seq_len(rows) %% 2
l_q3 <- matrix(q3[odd_cols == 1,],ncol=ncol(loss))
t_q3 <- matrix(q3[odd_cols == 0,],ncol=ncol(loss))


#pdf(file='err_eps.pdf') 
par(mar=c(5, 5, 5, 2) + 0.1)
cl <- brewer.pal(2*nrow(l_q2), "Set1")
cl[6]=cl[7]
pc <- c(0,1,2,5,6,11)
plot(0.6, 0.5, xlab=expression(epsilon),
     ylab="error",
     xlim = c(min(range),0.5), ylim = c(min(l_q2),max(l_q2)+2),
     type = "n", main=bquote(atop(bold("Error vs. contamination rate"),atop("n=100   p=10","Contamination by smallest eigen vector"))), 
     cex.main=2.2, cex.lab=2.2, cex.axis=2)	
	for (i in 1:nrow(l_q2)){
		lines(x= range, y=l_q2[i,], type='o', col=cl[i], pch=pc[i], lwd = 2.5, cex=2)
		arrows(x0=range, y0=l_q1[i,], x1=range, y1=l_q3[i,], code=3, col=cl[i], lwd=2, length=0)
	}
legend("topleft", legend = estimators, 
              col=c(cl, "black"),
              pch=pc, cex=2, pt.cex = 2, lwd=2) # optional legend
#dev.off()

