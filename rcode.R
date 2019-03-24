
#--------------RUN-THIS-FIRST--------------------------------
##### Simulations:
# This function will generate our data for us base 
# on n (sample length), theta (true parameter).
generate_dat <- function(n, theta){
  # generating a random (Cauchy) sample
  X = rcauchy(n, location = theta, scale = 1) 
  return(X)
}

# Define a constant value for theta
theta = 2

### R-Code for MSE
# The MSE takes in a vector then return a single value for MSE
MSE<-function(X,theta){
  MSE1 = mean((X-theta)^2)
  return(MSE1)
}

### R-Code for Coverage Probability

con_prob<-function(X,theta,epsilon){
  return(sum(abs(X-theta)<epsilon)/length(X))
}

#------------------------------------------------------------

## Sample Mean 

### Dynamic Plot for Estimator vs Sample Size

set.seed(21)
table_size=2 
plot_mean<-function(X,theta,table_size){
  n = length(X)
  # computes the sample mean of subsamples (X[1],...X[k]) for k=1,...,n
  submean = cumsum(X)/seq(1:n)
  plot(submean, xlab="k", ylab="Mean", type="l", main = n)  # plot of M1
  abline(h=theta)  # to compare with the true value 
}
par(mfrow=c(table_size,table_size))
plot_mean(generate_dat(100,theta),theta,table_size)
plot_mean(generate_dat(1000,theta),theta,table_size)
plot_mean(generate_dat(10000,theta),theta,table_size)
plot_mean(generate_dat(100000,theta),theta,table_size)

### MSE vs Sample Size

set.seed(20)
mean_MSE<-function(n){
  X = generate_dat(100,theta)
  submean = cumsum(X)/seq(1:100)
  MSE_vec = MSE(submean,theta)
  for(i in seq(110,n,10)){
    X = generate_dat(i,theta)
    submean = cumsum(X)/seq(1:i)
    MSE_vec = c(MSE_vec,MSE(submean,theta))
  }
  plot(seq(100,n,10),MSE_vec,
       col="blue",
       xlab = "sample size",
       ylab = "MSE", 
       type = "l", 
       main = "MSE vs sample size")
}


### Coverage vs Sample Size

set.seed(20)
mean_Con_Prob<-function(n){
  X = generate_dat(100,theta)
  submean = cumsum(X)/seq(1:100)
  co_prob_vec = con_prob(submean,theta,0.1)
  for(i in seq(110,n,10)){
    X = generate_dat(i,theta)
    submean = cumsum(X)/seq(1:i)
    co_prob_vec = c(co_prob_vec,con_prob(submean,theta,0.1))
  }
  plot(seq(100,n,10),co_prob_vec,
       col="blue",
       xlab = "sample size",
       ylab = "Coverage Probability", 
       type = "l", 
       main = "Coverage Probability vs sample size")
  abline(h=theta)
}

## Plotting the MSE and coverage of sample mean
set.seed(20)
par(mfrow=c(1,2))
mean_MSE(1000)
mean_Con_Prob(1000)

## Median

### Dynamic Plot for Estimator vs Sample Size

set.seed(20)
plot_median<-function(X,theta){
  n = length(X)
  # This function produce a vector os median sub-samples.
  Median = function(X) {z=vector(length=n); for(k in 1:n) z[k]=median(X[1:k]); z}
  M2=Median(X)
  plot(M2, xlab="k", ylab="Median", type="l", main = n)
  abline(h=theta)
}
par(mfrow=c(table_size,table_size))
plot_median(generate_dat(100,theta),theta)
plot_median(generate_dat(1000,theta),theta)
plot_median(generate_dat(10000,theta),theta)
plot_median(generate_dat(100000,theta),theta)

### MSE vs Sample Size

set.seed(123)
median_MSE<-function(n){
  MSE_vec = rep(0,n)
  for(i in seq(1,n)){
    temp_vec = rep(0,100)
    for(j in seq(1,100)){
      X = generate_dat(i,theta)
      temp_vec[j]=median(X)
    }
    MSE_vec[i]=MSE(temp_vec,theta)
  }
  y = rep(0,n)
  for(i in seq(1,n)){
    y[i]=2/i
  }
  plot(seq(1,n),MSE_vec,
       ylim=c(0,2.5),
       xlab = "sample size",
       ylab = "MSE", 
       type = "l", 
       main = "MSE vs sample size",
       col='blue')
  lines(seq(1,n),y,col='red')
  lines(seq(1,n),rep(0,n),col='green')
  legend('topright', 
         legend=c("MSE of Estimator", "n/I(theta) {n* 1/Information}","Zero Line"),
         col=c("blue", "red",'green'), 
         lty=1:2, 
         cex=0.8)
}

### Coverage vs Sample Size

set.seed(200)
# n > 110
median_Cov<-function(m){
  cov_vec=rep(0,m)
  for(i in seq(1,m)){
    temp_vec = rep(1,100)
    for(j in seq(1,100)){
      X = generate_dat(i,theta)
      temp_vec[j]=median(X)
    }
    cov_vec[i]=con_prob(temp_vec,theta,0.1)
  }
  plot(seq(1,m),cov_vec,xlab = "Sample Size",
       ylab = "MSE", type = "l", main = "Coverage Probability vs Sample Size",col='blue')
  lines(seq(1,m),rep(1,m),col='red')
  legend('bottomright', 
         legend=c("Coverage Probability", "1 Line"),
         col=c("blue", "red"), 
         lty=1:2, 
         cex=0.8)
}

### Plotting MSE and coverage vs sample size
par(mfrow=c(1,2))
set.seed(123)
median_MSE(100)
set.seed(200)
median_Cov(1000)


## Modified Sample Median

### Dynamic Plot vs Sample Size
set.seed(20)
plot_adj_median<-function(n,theta){
  # This function produce a vector of median sub-samples.
  M2 = rep(0,n)
  for(i in seq(1,n)){
    X = generate_dat(i,theta)
    M2[i] = median(X)+(4/i)*(sum((X-median(X))/(1+(X-median(X))^2)))
  }
  plot(M2, xlab="k", ylab="Median", type="l", main = n, col="blue")
  abline(h=theta, col="red")
}
plot_adj_median(10000,theta)

### MSE vs Sample Size
set.seed(5)
theta = 2
# n > 110
plot_Adj_MSE<-function(m){
  MSE_vec = rep(0,m)
  for(i in seq(1,m)){
    temp_vec = rep(0,100)
    for(j in seq(1,100)){
      X = generate_dat(i,theta)
      temp_vec[j]=median(X) + (4/i)*sum((X-median(X))/(1+(X-median(X))^2))
    }
    MSE_vec[i]=MSE(temp_vec,theta)
  }
  y = rep(0,m)
  for(i in seq(1,m)){
    y[i]=2/i
  }
  plot(seq(1,m),MSE_vec,
       ylim=c(0,2.5),
       xlab = "sample size",
       ylab = "MSE", type = "l", 
       main = "MSE vs sample size", 
       col='blue')
  lines(seq(1,m),y,col='red')
  lines(seq(1,m),rep(0,m),col='green')
  legend('topright', 
         legend=c("MSE of Estimator", "n/I(theta) {n* 1/Information}","Zero Line"),
         col=c("blue", "red",'green'), lty=1:2, cex=0.8)
}

### Coverage Probability vs Sample Size
set.seed(20)
theta = 2
coverage_mod_med<-function(m){
  co_prob_vec = rep(0,m)
  for(i in seq(1,m)){
    temp_vec = rep(0,100)
    for(j in seq(1,100)){
      Y = generate_dat(i,theta)
      temp_vec[j]=median(Y) + (4/i)*sum((Y-median(Y))/(1+(Y-median(Y))^2))
    }
    co_prob_vec[i]=con_prob(temp_vec,theta,0.1)
  }
  plot(seq(1,m),co_prob_vec,xlab = "sample size",
       ylab = "coverage probability", type = "l", 
       main = "Coverage Probability vs sample size",col='blue')
  lines(seq(1,m),rep(1,m),col='red')
}

### Plotting MSE and coverage vs sample size
par(mfrow=c(1,2))
set.seed(5)
plot_Adj_MSE(100)
set.seed(20)
coverage_mod_med(1000)


## Maximum Likelihood Estimator

### Uniqueness of Roots
X = c(1,6)
logL <- function(theta_t){sum(-log1p(1+(X-theta_t)^2))}
# Need to vectorize our logL else it wouldn't work.
y = Vectorize(logL)
plot(y, xlab="theta", ylab="log-likelihood",
     from=0, 
     to=7, 
     type="l")
abline(v=0.5*(7 + sqrt(21)))
abline(v=0.5*(7 - sqrt(21)))

### Dynamic Plot vs Sample Size
set.seed(20)
plot_MSE<-function(n){
  X = generate_dat(n,theta)
  Y = sample(X,1)
  logL <- function(theta_t){-sum(-log1p(1+(Y-theta_t)^2))}
  a = optim(1,logL, method = "Brent",lower = 0, upper = 5)
  max_vec = a$par
  for (i in 2:n){
    Y = sample(X,i)
    logL <- function(theta_t){-sum(-log1p(1+(Y-theta_t)^2))}
    a = optim(1,logL, method = "Brent", lower = 0, upper = 5)
    max_vec = c(max_vec,a$par)
  }
  plot(max_vec, xlab="k", ylab="MLE", type="l", main = n)
  abline(h=theta)
}
plot_MSE(10000)

### MSE vs Sample Size
set.seed(20)
MLE_MSE<-function(n){
  MSE_vec = rep(0,n)
  for(i in seq(1,n)){
    temp_vec = rep(0,100)
    for(j in seq(1,100)){
      X = generate_dat(i,theta)
      logL <- function(theta_t){-sum(-log1p(1+(X-theta_t)^2))}
      a = optim(1,logL, method = "Brent",lower = 0, upper = 5)
      temp_vec[j]=a$par
    }
    MSE_vec[i]=MSE(temp_vec,theta)
  }
  y = rep(0,n)
  for (i in seq(1,n)){
    y[i]=2/i
  }
  plot(seq(1,n),MSE_vec,
       ylim=c(0,2.5),
       xlab = "sample size",
       ylab = "MSE", 
       type = "l", 
       main = "MSE vs sample size", 
       col='blue')
  lines(seq(1,n),y, col = 'red')
  lines(seq(1,n),rep(0,n),col='green')
  legend('topright', legend=c("MSE of Estimator", "n/I(theta)","Zero Line"),
         col=c("blue", "red","green"), lty=1:2, cex=0.8)
}

### Coverage probability for MLE
set.seed(20)
# n > 110
Cov_MSE<-function(n){
  cov_vec = rep(0,n)
  for(i in seq(1,n)){
    temp_vec = rep(0,100)
    for(j in seq(1,100)){
      X = generate_dat(i,theta)
      logL <- function(theta_t){-sum(-log1p(1+(X-theta_t)^2))}
      a = optim(1,logL, method = "Brent",lower = 0, upper = 5)
      temp_vec[j]=a$par
    }
    cov_vec[i]= con_prob(temp_vec,theta,0.1)
  }
  plot(seq(1,n),cov_vec,xlab = "sample size",
       ylab = "Coverage Probability", type = "l", 
       main = "Coverage Probability vs sample size", col='blue')
  lines(seq(1,n),rep(1,n),col='red')
  legend('bottomright', 
         legend=c("Coverage Probability of Estimator", "1 Line"),
         col=c("blue","red"), 
         lty=1:2, 
         cex=0.8)
}

### Plotting MSE and Coverage probability
par(mfrow=c(1,2))
set.seed(20)
MLE_MSE(100)
set.seed(20)
Cov_MSE(1000)

