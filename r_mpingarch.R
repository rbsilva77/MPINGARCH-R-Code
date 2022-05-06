#------------------------------------------
# R code for generating MPINGARCH processes
#------------------------------------------

# Instructions:

# Set the theoretical values for the parameters

# Set the sample size

r.pigingarch <- function(size, order, par){
  
  pacman::p_load(tseries, gamlss.dist) # R libraries
  
  # Generating count responses from PIG-INGARCH process
  
  r.pigingarch1 <- function(size, par){  #PIGINGARCH(1,0) 
    
    alpha0 <- par[1]
    alpha1 <- par[2]
    phi <- par[3]
    
    n <- size
    
    x <- NULL
    mu <- NULL
    
    mu[1] <- alpha0 / (1 - alpha1)
    x[1] <- rPIG(n = 1, mu = mu[1], sigma = 1 / phi)
    
    for (i in 1:(n-1)){
      mu[i+1] <- alpha0 + alpha1 * x[i]
      x[i+1] <- rPIG(n = 1, mu = mu[i+1], sigma = 1 / phi)
    }
    
    return(x)
  }
  
  r.pigingarch2 <- function(size , par){ #PIGINGARCH(2,0)
    
    alpha0 <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    phi <- par[4]
    
    n <- size
    
    x <- NULL
    mu <- NULL
    
    mu[1] <- alpha0 / (1 - alpha1 - alpha2)
    x[1] <- rPIG(n = 1, mu = mu[1], sigma = 1 / phi)
    
    mu[2] <- alpha0 + alpha1 * x[1] + alpha2 * rPIG(n = 1, mu = mu[1], sigma = 1 / phi)
    x[2] <- rPIG(n = 1, mu = mu[2], sigma = 1/phi)
    
    for (i in 2:(n-1)){
      mu[i+1] <- alpha0 + alpha1 * x[i] + alpha2 * x[i-1]
      x[i+1] <- rPIG(n = 1, mu = mu[i+1], sigma = 1/phi)
    }
    
    return(x)
  }
  
  r.pigingarch11 <- function(size, par){ #PIGINGARCH(1,1)
    
    alpha0 <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    beta1 <- par[4]
    phi <- par[5]
    
    n <- size
    
    x <- NULL
    mu <- NULL
    
    mu[1] <- alpha0 / (1 - alpha1 - beta1)
    x[1] <- rPIG(n = 1, mu = mu[1], sigma = 1 / phi, max.value = 10000)
    
    for (i in 1:(n-1)){
      mu[i+1] <- alpha0 + alpha1 * x[i] + beta1 * mu[i]
      x[i+1] <- rPIG(n = 1, mu = mu[i+1], sigma = 1 / phi, max.value = 10000)
    }
    
    return(x)
  }
  
  switch(order,
         "1" = r.pigingarch1(size, par),
         "2" = r.pigingarch2(size, par),
         "11" = r.pigingarch11(size, par)
  )
  
}

# Negative binomial case

r.nbingarch <- function(size, order, par){
  
  pacman::p_load(tseries, gamlss.dist) # R libraries
  
  # Generating count responses from NB-INGARCH process
  
  r.nbingarch1 <- function(size, par){  #NBINGARCH(1,0) 
    
    alpha0 <- par[1]
    alpha1 <- par[2]
    phi <- par[3]
    
    n <- size
    
    x <- NULL
    mu <- NULL
    
    mu[1] <- alpha0 / (1 - alpha1)
    x[1] <- rnbinom(n = 1, size = phi, prob = phi / (mu[1] + phi))
    
    for (i in 1:(n-1)){
      mu[i+1] <- alpha0 + alpha1 * x[i]
      x[i+1] <- rnbinom(n = 1, size = phi, prob = phi / (mu[i+1] + phi))
    }
    
    return(x)
  }
  
  r.nbingarch2 <- function(size , par){ #NBINGARCH(2,0)
    
    alpha0 <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    phi <- par[4]
    
    n <- size
    
    x <- NULL
    mu <- NULL
    
    mu[1] <- alpha0 / (1 - alpha1 - alpha2)
    x[1] <- rnbinom(n = 1, size = phi, prob = phi / (mu[1] + phi))
    
    mu[2] <- alpha0 + alpha1 * x[1] + alpha2 * rnbinom(n = 1, size = phi, prob = phi / (mu[1] + phi))
    x[2] <- rnbinom(n = 1, size = phi, prob = phi/(mu[2] + phi))
    
    for (i in 2:(n-1)){
      mu[i+1] <- alpha0 + alpha1 * x[i] + alpha2 * x[i-1]
      x[i+1] <- rnbinom(n = 1, size = phi, prob = phi / (mu[i+1] + phi))
    }
    
    return(x)
  }
  
  r.nbingarch11 <- function(size, par){ #NBINGARCH(1,1)
    
    alpha0 <- par[1]
    alpha1 <- par[2]
    alpha2 <- par[3]
    beta1 <- par[4]
    phi <- par[5]
    
    n <- size
    
    x <- NULL
    mu <- NULL
    mu[1] <- alpha0 / (1 - alpha1 - beta1)
    x[1] <- rnbinom(n = 1, size = phi, prob = phi / (mu[1] + phi))
    for (i in 1:(n-1)){
      mu[i+1] <- alpha0 + alpha1 * x[i] + beta1 * mu[i]
      x[i+1] <- rnbinom(n = 1, size = phi, prob = phi / (mu[i+1] + phi))
    }
    
    return(x)
  }
  
  switch(order,
         "1" = r.nbingarch1(size, par),
         "2" = r.nbingarch2(size, par),
         "11" = r.nbingarch11(size, par)
  )
  
}





