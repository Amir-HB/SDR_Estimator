source("SDR.R")
source("IRM.R")
source("Contaminations.R")
source("Experiment_materiels.R")

n = 100
p = 10
#n = 10*floor(p/eps^2)
eps = .2
Sigma = diag(p)
parameter = "p" # choices: n, p, eps
#range = c(seq(.04,.48,.04))
range = c(seq(10,60,10))
#range = c(seq(50,300,50))
contamination = "flip" # existing contaminations: "mixt", "unif", "mixt"
estimators= c("SDR"
				#,"Iterative filtering"
				,"Geometric median"
				#,"Oracle"
				#,"Sample mean"
				#,"Iterative reweighting"
				#,"Tukey's median"
				#,"New iterative reweighting"
				,"Componentwise median"
				) 
#existing estimators: SDR", "Iterative filtering","Geometric median","Oracle","Sample mean", "Iterative reweighting","Tukey's median"

rep_num = 50
loss = variation(n,p,eps,Sigma,parameter,range,rep_num,contamination,estimators)

estimators2 = c(#"SDR"
				#,"Iterative filtering"
				#,"Geometric median"
				#,"Oracle"
				#,"Sample mean"
				"Iterative reweighting"
				#,"Tukey's median"
				) 
rep_num = 10
loss_2 = variation(n,p,eps,Sigma,parameter,range,rep_num,contamination,estimators2)
loss = rbind(loss,loss_2)
estimators = c(estimators,estimators2)

#write.table(rbind(range,loss), "intro2", row.names = FALSE, col.names = FALSE)
#Tabl = read.table("intro2")
#loss = as.matrix(Tabl)
#range = loss[1,]
#loss = loss[-1,]

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

#pdf(file='intro2.pdf') ; par(mar=c(4, 5, 4, 2) + 0.1)
cl <- brewer.pal(2*nrow(l_q2), "Set1")
cl[6]=cl[7]
pc <- c(0,1,2,5,6,11)
plot(0.6, 0.5, xlab=p,
     ylab="error",
     xlim = c(min(range),0.5), ylim = c(min(l_q2)-1 ,max(l_q2)+3),
     type = "n", main=paste("Error vs dimension \n n=100", expression(epsilon),"=0.2"), 
     cex.main=2.2, cex.lab=2.2, cex.axis=2)	
	for (i in 1:nrow(l_q2)){
		lines(x= range, y=l_q2[i,], type='o', col=cl[i], pch=pc[i-1], lwd = 2.5, cex=2)
		arrows(x0=range, y0=l_q1[i,], x1=range, y1=l_q3[i,], code=3, col=cl[i], lwd=2, length=0)
	}
legend("topleft", legend = estimators, 
              col=c(cl, "black"),
              pch=pc, cex=2, pt.cex = 2, lwd=2) # optional legend
#dev.off()
