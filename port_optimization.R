my_files=list.files(pattern="^[0-9]{1,6}.csv")
n=length(my_files)
all_price=c()
for(f in my_files)
{
data=read.csv(f)
price=data[,"Close.Price"]
all_price=cbind(all_price,price)
}
colnames(all_price)=substr(my_files,1,6)
d=nrow(all_price)
daily_return=(all_price[1:(d-1),]-all_price[2:d,])/all_price[2:d,]
C=cov(daily_return)

w = rep(1/n,n)
risk <- 100*sqrt(t(w) %*% C %*% w)

w <- c(0.1, 0.125, 0.225, 0.3, 0.25)
risk

w <- c(0.25, 0.15, 0.10, 0.30, 0.20)

#solving system of linear equation

A <- cbind(C, rep(1,n))
A <- rbind(A, c(rep(1,n),0))
round((A),4)
b <- c(rep(0,n),1)
w <- solve(A,b)[1:n] #ignoring lambda
risk <- 100*sqrt(t(w) %*% C %*% w)
risk

mu=(all_price[1,]-all_price[d,])/all_price[d,]
mu


#markowitz curve

w1 <- seq(-1, 2, length = 5000)
w2 <- 1-w1
all_w <- rbind(w1,w2)
sigma_square <- rep(0,5000)
r <- rep(0,5000)
for (i in 1:5000)
{
  w = all_w[,i]
  sigma_square[i] <- t(w)%*%C%*%w
  r[i] <- sum(mu*w)
}

plot(sigma_square, r, type = "l")

#sigma_square <- t(w)%*%C%*%w
#r <- colSums(mu*w)
#C
#t(w)
#dim(sigma_square)

sigma1 <- sqrt(C[1,1])
sigma2 <- sqrt(C[2,2])
rho <- C[1,2] / (sigma1*sigma2)
w1_star <- (sigma2^2 - sigma1*sigma2*rho)/(sigma1^2+sigma2^2-2*sigma1*sigma2*rho)
w2_star <- 1-w1_star
w_star <- c(w1_star, w2_star)
r_star <- sum(mu*w_star)
sigma_square_star <- t(w_star)%*%C%*%w_star
points(sigma_square_star, r_star, col = "red", lwd = 3)
idx <- which.max(r/sqrt(sigma_square))
w_sr_opt <- all_w[,idx] #highest sharpe ratio
r_sr_opt <- sum(mu*w_sr_opt)
sigma_square_sr_opt <- t(w_sr_opt) %*% C %*% w_sr_opt
points(sigma_square_sr_opt, r_sr_opt, col = "blue", lwd = 3)




install.packages("quadprog")
require(quadprog)

#optimization with more constraints

#1. historical


mu=(all_price[1,]-all_price[d,])/all_price[d,]


#D <- C
d <- rep(0,n)
A <- matrix(c(rep(-1,n),mu, diag(n)), nc=n, byrow=T)
rf <- 0.08
b <- c(-1, rf, rep(0,n))
solve.QP(C, d, t(A), b, meq = 0)

sqrt(2*ax$value)*100 #this is 1/2 * t(W) * C *W




