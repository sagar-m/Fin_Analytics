#Monte-carlo simulation
#1. expectation of random variable
N=100000
z=runif(N,min=0,max=2)
y=function(z) pmax(z,rep(1,N))
mean(y(z))

#2. Integration
N=1000000
z=runif(N,min=0,max=3)
y=function(z) z^2
mean(y(z))*(3-0)

N=1000000
z=runif(N,min=0,max=3)
y=function(z) exp(z^2)
mean(y(z))*(3-0)

#3. Value of pi

N=1000000
x=runif(N,min=0,max=2)
y=runif(N,min=0,max=2)
d=sqrt((x-1)^2+(y-1)^2)
N1=sum(d<1)
4*N1/N

#4. Expected Stock price simulation (independent case)
data=read.csv("532667.csv")
price=data[,"Close.Price"]
d=length(price)
r=(price[1:(d-1)]-price[2:d])/price[2:d]
N=1000000
price_y=data[1:249,"Close.Price"]
d=length(price_y)
r=(price_y[1:(d-1)]-price_y[2:d])/price_y[2:d]
mu=mean(r)*d
sigma=sd(r)
t=1/12
S0=price_y[1]
z=rnorm(N)
St=S0*exp((mu+(sigma^2)/2)*t+sigma*sqrt(t)*z)
mean(St)


#5. Stock price simulation (correlated case)
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
B=chol(C)
S0=all_price[1,]
mu=colMeans(daily_return)*d
sigma=sqrt(diag(C))
z=rnorm(n)
St=S0*exp((mu+(sigma^2)/2)*t+sqrt(t)*B%*%z)
simulations=c()
N=1000000

t1=proc.time()
for(i in 1:N)
{
z=rnorm(n)
St=S0*exp((mu+(sigma^2)/2)*t+sqrt(t)*B%*%z)
simulations=cbind(simulations,St)
}
t2=proc.time()
t2-t1

t1=proc.time()

z=matrix(rnorm(n*N),nr=n,nc=N)
St=S0*exp((mu+(sigma^2)/2)*t+sqrt(t)*B%*%z)
t2=proc.time()
t2-t1
rowMeans(St)


#6. VaR using MC
theta=rep(1/n,n)
W0=sum(S0*theta)
Wt=colSums(St*theta)
probability=sum(Wt/W0>1.02)/N



#7. Option price using MC

t=25/365
z=matrix(rnorm(n*N),nr=n,nc=N)
St=S0*exp((mu+(sigma^2)/2)*t+sqrt(t)*B%*%z)
ongc=St[2,]
K=172.50
mean(pmax(ongc-K,0))

#testing using BSM
d1=	(log(S0[2]/K)+(mu[2]+(sigma[2]^2)/2)*t)/(sigma[2]*sqrt(t))
d2=d1-sigma[2]*sqrt(t)
S0[2]*pnorm(d1)-pnorm(d2)*K*exp(-mu[2]*t)


#8. Greeks
#delta
ds=.1
ps=pmax(ongc-K,0)
psds=pmax(ongc+ds-K,0)
delta=mean((psds-ps)/ds)
#gamma
ds=.1
ps=pmax(ongc-K,0)
psds=pmax(ongc+ds-K,0)
ps2ds=pmax(ongc+2*ds-K,0)

gamma=mean((ps2ds+ps-2*psds)/(ds^2))

#vega
t=25/365
z=matrix(rnorm(n*N),nr=n,nc=N)
dv=.001
St=S0*exp((mu+(sigma^2)/2)*t+sqrt(t)*B%*%z)
Stdv=S0*exp((mu+((sigma+dv)^2)/2)*t+sqrt(t)*(B+dv)%*%z)

ongc=St[2,]
ongcdv=Stdv[2,]
K=172.50
dv=.1
pv=pmax(ongc-K,0)
pvdv=pmax(ongcdv-K,0)
mean((pvdv-pv)/dv



