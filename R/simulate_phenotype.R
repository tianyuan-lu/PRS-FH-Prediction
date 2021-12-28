set.seed(1)
library(mvtnorm)

alpha <- 0.3
beta <- 0.8

P <- rmvnorm(10000, sigma = matrix(c(1,0,0.5,0,1,0.5,0.5,0.5,1),3,3))
P1 <- P[,1]
P2 <- P[,2]
P3 <- P[,3]

G <- rmvnorm(10000, sigma = matrix(c(1,0,0.5,0,1,0.5,0.5,0.5,1),3,3))
G1 <- G[,1]
G2 <- G[,2]
G3 <- G[,3]

Y1 <- alpha * P1 + beta * G1 + rnorm(10000, sd = sqrt(1 - alpha ^ 2 - beta ^ 2))
Y2 <- alpha * P2 + beta * G2 + rnorm(10000, sd = sqrt(1 - alpha ^ 2 - beta ^ 2))
Y3 <- alpha * P3 + beta * G3 + rnorm(10000, sd = sqrt(1 - alpha ^ 2 - beta ^ 2))
Y_midp <- (Y1 + Y2) / 2

Z1 <- sapply(1:length(Y1), function(i) rbinom(1,1,exp(-1 + Y1[i]) / (1 + exp(-1 + Y1[i]))))
Z2 <- sapply(1:length(Y2), function(i) rbinom(1,1,exp(-1 + Y2[i]) / (1 + exp(-1 + Y2[i]))))
Z3 <- sapply(1:length(Y3), function(i) rbinom(1,1,exp(-1 + Y3[i]) / (1 + exp(-1 + Y3[i]))))

cor(Y3,P3)^2 
cor(Y3,Y_midp)^2

summary(glm(Z3 ~ P3, family = binomial))
summary(glm(Z3 ~ Z1 + Z2, family = binomial))


