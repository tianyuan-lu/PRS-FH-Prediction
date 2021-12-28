args = commandArgs(trailingOnly = T)
h2_PRS <- as.numeric(args[1])
h2_midp <- as.numeric(args[2])
file <- read.table(args[3],header = T)

alpha = sqrt(h2_PRS)
beta = sqrt(sqrt(2 * h2_midp) - alpha^2)

Sigma_geno <- matrix(c(1,0,1/2,0,1,1/2,1/2,1/2,1),3,3)
Sigma_PG <- matrix(rep(0,9),3,3)
Sigma_PY <- matrix(c(alpha,0,alpha/2,0,alpha,alpha/2,alpha/2,alpha/2,alpha),3,3)
Sigma_GY <- matrix(c(beta,0,beta/2,0,beta,beta/2,beta/2,beta/2,beta),3,3)
Sigma_pheno <- matrix(c(1,0,(alpha^2 + beta^2)/2,0,1,(alpha^2 + beta^2)/2,(alpha^2 + beta^2)/2,(alpha^2 + beta^2)/2,1),3,3)

Sigma <- rbind(cbind(Sigma_geno, Sigma_PG, Sigma_PY),
               cbind(Sigma_PG, Sigma_geno, Sigma_GY),
               cbind(Sigma_PY, Sigma_GY, Sigma_pheno))
rownames(Sigma) <- colnames(Sigma) <- c("P_1","P_2","P_3","G_1","G_2","G_3","Y_1","Y_2","Y_3")

directFun <- function(Y_1, Y_2, P_3) {
  targetidx <- c(1,2,4,5,6,9)
  mean = Sigma[targetidx, -targetidx] %*% solve(Sigma[-targetidx, -targetidx]) %*% as.matrix(c(P_3,Y_1,Y_2))
  return(mean)
}

directFunParent <- function(Y_1, P_3) {
  targetidx <- c(1,2,4,5,6,8,9)
  mean = Sigma[targetidx, -targetidx] %*% solve(Sigma[-targetidx, -targetidx]) %*% as.matrix(c(P_3,Y_1))
  return(mean)
}

file$Joint <- NA
if ("Mother" %in% colnames(file) & "Father" %in% colnames(file)) {
  for (ind in 1:nrow(file)) {
    if (ind %% 100 == 0) {
      cat("processed",ind,"samples...",collapse="\n")
    }
    file$Joint[ind] <- directFun(Y_1 = file$Mother[ind],
                                 Y_2 = file$Father[ind],
                                 P_3 = file$PRS[ind])[6]
  }
}
if ("Mother" %in% colnames(file) & !"Father" %in% colnames(file)) {
  for (ind in 1:nrow(file)) {
    if (ind %% 100 == 0) {
      cat("processed",ind,"samples...",collapse="\n")
    }
    file$Joint[ind] <- directFunParent(Y_1 = file$Mother[ind],
                                       P_3 = file$PRS[ind])[7]
  }
}
if (!"Mother" %in% colnames(file) & "Father" %in% colnames(file)) {
  for (ind in 1:nrow(file)) {
    if (ind %% 100 == 0) {
      cat("processed",ind,"samples...",collapse="\n")
    }
    file$Joint[ind] <- directFunParent(Y_1 = file$Father[ind],
                                       P_3 = file$PRS[ind])[7]
  }
}
save(file,file = paste0(args[3],".joint.pred.continuous.RData"))


