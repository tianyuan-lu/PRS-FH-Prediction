set.seed(1)
X <- rnorm(1000)
E <- rnorm(1000)
Y <- exp(0.5 * X + E - 3) / (1 + exp(0.5 * X + E - 3))
Z <- sapply(1:length(X), function(i) rbinom(1,1,Y[i]))
mod <- glm(Z~X + E, family = binomial)



