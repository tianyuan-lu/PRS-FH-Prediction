library(pracma) ## please install this package

estimate_h2 <- function(hatmu_0, hatbeta, degree = 1, upperbound = 2) {
  ilogit <- function(x) {
    exp(x) / (exp(x) + 1)
  }
  jointfun <- function(x,y,mu_0,z1,z3,beta) {
    p1 = (exp(x + mu_0) / (1 + exp(x + mu_0))) ^ z1 * (1 - exp(x + mu_0) / (1 + exp(x + mu_0))) ^ (1 - z1)
    p2 = (exp(y + mu_0) / (1 + exp(y + mu_0))) ^ z3 * (1 - exp(y + mu_0) / (1 + exp(y + mu_0))) ^ (1 - z3)
    p3 = 0.5 / pi / sqrt((1 - (2^-degree)^2) * ( beta^2)^2)
    p4 = exp(- (x^2 - x * y / (2^(degree - 1)) + y^2) / (2 * ((1 - (2^-degree)^2) * ( beta^2))))
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
  Flogit <- function(x) {
    1 / (1 + exp(-x * pi / sqrt(3)))
  }
  trans.b.tau <- function(b,b0) {
    a <- qnorm(Flogit(b + b0)) - qnorm(Flogit(b0))
    out <- a/sqrt(1 + (a^2))
    return(out)
  }
  hatb <- searchfun(hatbeta)
  hattau <- trans.b.tau(hatb, hatmu_0)
  cat("Estimated effect size for genetic component (std. normal distn.) - logistic model: ",hatb,collapse="\n")
  cat("Estimated effect size for genetic component (std. normal distn.) - liability scale: ",hattau,collapse="\n")
  cat("Estimated heritability: ",hattau^2 * 100,"%",collapse="\n")
}
