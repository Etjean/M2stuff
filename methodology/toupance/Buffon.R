

#------------------------------------------------------------------------------
# Estimate pi using Buffon's needle method 
#  -n- needle drops
#  -r- ratio between length of needle L and width of board D: r=L/D
#  returns -one- estimation of pi using Buffon's needle estimator (n, r)

estim_pi <- function(n, r) {
	x1 <- runif(n, min=0, max=1)
	theta <- runif(n, min=-pi/2, max=pi/2)
	x2 <- x1 + sin(theta)*r
	ns <- sum(x2 > 1 | x2 < 0)
	return(2*n/ns*r)
}




#------------------------------------------------------------------------------
# Compute -one- estimation of pi using n=100 and r=0.3

n <- 100
r <- 0.3
estim_pi(n, r)




#------------------------------------------------------------------------------
# Compute -nrep- independant estimations of pi
#  to characterize Buffon's needle estimator (n=100 and r=0.3)
#  repeat estimation nrep=100000 times
#  store results in -res-
#  draw histogram 
#  compute mean, variance and mse 

nrep <- 100000
res <- rep(NA, nrep)
for (i in 1:nrep) {
	res[i] <- estim_pi(n, r)
}

require(MASS)
truehist(res, col="gray50", xlab="pi_hat", ylab="density")
abline(v=pi, col="red", lwd=2)

m <- mean(res)
s2 <- var(res)
mse <- sum((res-pi)^2)/nrep
print(m)
print(s2)
print(mse)

abline(v=m, col="blue", lwd=2, lty=2)

legend("topright", lwd=2, col=c("red", "blue"), lty=c(1, 2), legend=c("pi", "mean(pi_hat)"))




#------------------------------------------------------------------------------
# Compute -mse- of Buffon's needle estimator (n, r)
#  using -nrep- independant simulations
#  store results in -res-
#  draw histogram 
#  compute mean, variance and mse 

mse_pi <- function(nrep, n, r) {
	res <- rep(NA, nrep)
	for (i in 1:nrep) {
		res[i] <- estim_pi(n, r)
	}
	return(sum((res-pi)^2)/nrep)
}



#------------------------------------------------------------------------------
# Compute -mse- of Buffon's needle estimator (n=100, r)
#  using nrep=10000 independant simulations
#  for r ratio varying from 0.1 to 2.0

n <- 100
nrep <- 100000

r_vec <- seq(from=0.1, to=2.0, by=0.05)
r_vec <- seq(from=0.95, to=1.1, len=20)
m <- length(r_vec)
mse_vec <- rep(NA, m)
for (i in 1:m) {
	mse_vec[i] <- mse_pi(nrep, n, r_vec[i])
}

mse_vec[mse_vec == Inf] <- NA
print(mse_vec)
plot(r_vec, mse_vec, xlab="r", ylab="mse")



