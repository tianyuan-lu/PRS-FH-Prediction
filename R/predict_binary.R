args = commandArgs(trailingOnly = T)
hatmu_0 <- as.numeric(args[1])
hatalpha <- as.numeric(args[2])
hatbeta <- as.numeric(args[3])
hatsex <- as.numeric(args[4])
file <- read.table(args[5],header = T)

library(pracma)
hatb <- Inf
upperbound <- 0.9
while(hatb >= upperbound) {
  upperbound = upperbound + 0.1
  jointfun <- function(x,y,mu_0,z1,z3,beta) {
    p1 = (exp(x + mu_0) / (1 + exp(x + mu_0))) ^ z1 * (1 - exp(x + mu_0) / (1 + exp(x + mu_0))) ^ (1 - z1)
    p2 = (exp(y + mu_0) / (1 + exp(y + mu_0))) ^ z3 * (1 - exp(y + mu_0) / (1 + exp(y + mu_0))) ^ (1 - z3)
    p3 = 0.5 / pi / sqrt((1 - (2^-1)^2) * ( beta^2)^2)
    p4 = exp(- (x^2 - x * y / (2^(1 - 1)) + y^2) / (2 * ((1 - (2^-1)^2) * ( beta^2))))
    return(p1 * p2 * p3 * as.numeric(p4))
  }
  betafun <- function(beta) {
    log(integral2(function(x,y) {jointfun(x,y,z1 = 1,z3 = 1,beta = beta,mu_0 = hatmu_0)}, xmin = -10, xmax = 10, ymin = -10, ymax = 10)$Q *
          integral2(function(x,y) {jointfun(x,y,z1 = 0,z3 = 0,beta = beta,mu_0 = hatmu_0)}, xmin = -10, xmax = 10, ymin = -10, ymax = 10)$Q /
          integral2(function(x,y) {jointfun(x,y,z1 = 1,z3 = 0,beta = beta,mu_0 = hatmu_0)}, xmin = -10, xmax = 10, ymin = -10, ymax = 10)$Q /
          integral2(function(x,y) {jointfun(x,y,z1 = 0,z3 = 1,beta = beta,mu_0 = hatmu_0)}, xmin = -10, xmax = 10, ymin = -10, ymax = 10)$Q)
  }
  initialseq <- 50
  finalseq <- 1000
  searchfun <- function(targetbeta) {
    maxbeta <- seq(1e-3,upperbound,length.out = initialseq)[which.max(sapply(seq(1e-3,upperbound,length.out = initialseq),betafun))]
    beta <- seq(1e-3,maxbeta,length.out = finalseq)[which.min(abs(targetbeta - sapply(seq(1e-3,maxbeta,length.out = finalseq),betafun)))]
    return(beta)
  }
  hatb <- searchfun(hatbeta)
}
set.seed(1)

denomfun <- function(x) {
  (exp(hatmu_0 + x) / (1 + exp(hatmu_0 + x))) * exp(-x^2 / 2 / (hatalpha^2 + hatb^2))
}
denom <- integral(denomfun,-10,10)
pdffun1 <- function(x) {
  (exp(hatmu_0 + x) / (1 + exp(hatmu_0 + x))) * exp(-x^2 / 2 / (hatalpha^2 + hatb^2)) / denom
}
pdffun0 <- function(x) {
  (1 / (1 + exp(hatmu_0 + x))) * exp(-x^2 / 2 / (hatalpha^2 + hatb^2)) / denom
}

denomfunfather <- function(x) {
  (exp(hatmu_0 + x + hatsex) / (1 + exp(hatmu_0 + x + hatsex))) * exp(-x^2 / 2 / (hatalpha^2 + hatb^2))
}
denomfather <- integral(denomfunfather,-10,10)
pdffun1father <- function(x) {
  (exp(hatmu_0 + x + hatsex) / (1 + exp(hatmu_0 + x + hatsex))) * exp(-x^2 / 2 / (hatalpha^2 + hatb^2)) / denomfather
}
pdffun0father <- function(x) {
  (1 / (1 + exp(hatmu_0 + x + hatsex))) * exp(-x^2 / 2 / (hatalpha^2 + hatb^2)) / denomfather
}

support <- seq(-10,10,by = 0.001)
prob1 <- pdffun1(support)
prob1 <- prob1 / sum(prob1)
prob0 <- pdffun0(support)
prob0 <- prob0 / sum(prob0)

prob1father <- pdffun1father(support)
prob1father <- prob1father / sum(prob1father)
prob0father <- pdffun0father(support)
prob0father <- prob0father / sum(prob0father)

samp_fun1 <- function(n) sample(
  support, 
  n, 
  TRUE,  
  prob1
)

samp_fun0 <- function(n) sample(
  support, 
  n, 
  TRUE,  
  prob0
)

samp_fun1father <- function(n) sample(
  support, 
  n, 
  TRUE,  
  prob1father
)

samp_fun0father <- function(n) sample(
  support, 
  n, 
  TRUE,  
  prob0father
)

file$Joint_Binary <- NA
if ("Mother_Disease" %in% colnames(file) & "Father_Disease" %in% colnames(file)) {
  constant <- matrix(c((hatalpha^2 + hatb^2)/2,(hatalpha^2 + hatb^2)/2,hatalpha),1,3) %*% solve(matrix(c(hatalpha^2 + hatb^2, 0, hatalpha/2, 0, hatalpha^2 + hatb^2, hatalpha/2, hatalpha/2, hatalpha/2, 1),3,3))
  distn1father <- samp_fun1father(1e6)
  distn0father <- samp_fun0father(1e6)
  distn1mother <- samp_fun1(1e6)
  distn0mother <- samp_fun0(1e6)
  for (ind in 1:nrow(file)) {
    if (ind %% 100 == 0) {
      cat("processed",ind,"samples...",collapse="\n")
    }
    if (file$Father_Disease[ind] == 1) {
      simulatedfather = distn1father
    }
    else {
      simulatedfather = distn0father
    }
    if (file$Mother_Disease[ind] == 1) {
      simulatedmother = distn1mother
    }
    else {
      simulatedmother = distn0mother
    }
    scr <- rbind(simulatedfather, simulatedmother, rep(file$PRS[ind],1e6))
    file$Joint_Binary[ind] = mean(constant %*% scr)
  }
}
if ("Mother_Disease" %in% colnames(file) & !"Father_Disease" %in% colnames(file)) {
  constant <- matrix(c((hatalpha^2 + hatb^2)/2,hatalpha),1,2) %*% solve(matrix(c(hatalpha^2 + hatb^2,hatalpha/2,hatalpha/2,1),2,2))
  distn1 <- samp_fun1(1e6)
  distn0 <- samp_fun0(1e6)
  for (ind in 1:nrow(file)) {
    if (ind %% 100 == 0) {
      cat("processed",ind,"samples...",collapse="\n")
    }
    if (file$Mother_Disease[ind] == 1) {
      simulatedmother = distn1
    }
    else {
      simulatedmother = distn0
    }
    scr <- rbind(simulatedmother, rep(file$PRS[ind],1e6))
    file$Joint_Binary[ind] = mean(constant %*% scr)
  }
}
if (!"Mother_Disease" %in% colnames(file) & "Father_Disease" %in% colnames(file)) {
  constant <- matrix(c((hatalpha^2 + hatb^2)/2,hatalpha),1,2) %*% solve(matrix(c(hatalpha^2 + hatb^2,hatalpha/2,hatalpha/2,1),2,2))
  distn1 <- samp_fun1father(1e6)
  distn0 <- samp_fun0father(1e6)
  for (ind in 1:nrow(file)) {
    if (ind %% 100 == 0) {
      cat("processed",ind,"samples...",collapse="\n")
    }
    if (file$Father_Disease[ind] == 1) {
      simulatedfather = distn1
    }
    else {
      simulatedfather = distn0
    }
    scr <- rbind(simulatedfather, rep(file$PRS[ind],1e6))
    file$Joint_Binary[ind] = mean(constant %*% scr)
  }
}
save(file,file = paste0(args[5],".joint.pred.binary.RData"))





