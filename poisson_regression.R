source("Poisson_regression_functions.R")
eta_start <- log(rowMeans(data)+1)
mu_start <- mean(eta_start)
lambda <- 0.001
rep <- 1e5

#### IS samples ####

isbark <- mybarker(eta_start, mu_start, lambda, sigma_eta, iter = rep, delta = 0.003, data)
output_chain_bark <- list(isbark[[1]])

samps <- output_chain_bark[[1]]
wts <- isbark[[2]]

####### True samples #####

foo <- barker(eta_start, mu_start, sigma_eta, iter = rep, delta = 0.0012, data)

##### Estimating the marginal density using a univariate-Gaussian kernel

# j <- 51 ## for mu
j <- 5

h <- 0.08 # can tune
y <- seq(-5, 5, length = 1e3)


marg_samps <- samps[, j]


est <- numeric(length(y))


for(i in 1:length(y))
{
  pts <- y[i] - marg_samps
  kh <- (dnorm(pts/h, mean = 0, sd = 1))/h
  est[i] <- (kh %*% wts)/(sum(wts))
}




################# TRUE SAMPLES #############################


samps2 <- foo[[1]]
marg_samps_2 <- samps2[, 51]
est_2 <- numeric(length(y))
n <- length(marg_samps_2)

for(i in 1:length(y))
{
  pts_2 <- y[i] - marg_samps_2
  kh <- (dnorm(pts_2/h, mean = 0, sd = 1))/(n*h)
  est_2[i] <- sum(kh)
}

plot(y, est, type = "l", lwd = 2, col = "blue", ylim = c(0, 1.5),  ylab = "Density", xlab = expression(mu))
lines(y, est_2, col = "red", lwd = 2)
# lines(density(marg_samps_2), col = "green", lwd = 2)
legend("topright", legend = c("IS", "True"), col = c("blue", "red"), lwd = 2)

