#---------------------------------------
# R code for fitting MPINGARCH processes
#---------------------------------------

# Arguments of the fit.mpingarh function:

# x: count time series data
# type: PIG or NB distributions
# order: model order
# method: optimization method
# tol: tolerance level in the EM algorithm


fit.mpingarch <- function(x, type, order, method, tol){
 
  pacman::p_load(tseries, gamlss.dist) # R libraries
  
  # PIGINGARCH(1,0) case
  
  pigingarch1 <- function(x, method, tol){
    
    n <- length(x)
    
    # Auxiliary functions
    qsi0 <- -1/2
    B <- function(x) {-sqrt((-2 * x))}
    d2B <- function(x) {(2 * sqrt(2) * (-x)^(3 / 2))^(-1)}
    T <- function(x) {x}
    D <- function(x) {log(x) / 2}
    dD <- function(x) {1 / (2 * x)}
    H <- function(x) {-log(2 * pi * x^3) / 2}
    G <- function(x) {-1 / (2 * x)}
    
    # Conditional expectations
    eta_r <- function(x){
      
      y <- x[1:n]
      mu <- x[(n+1):(2*n)]
      phi <- x[2*n+1]
      
      (y + 1) * dPIG(y + 1, mu , 1 / phi) / (mu * dPIG(y, mu, 1 / phi))
    }
    
    kappa_r <- function(x){
      
      y<-x[1:n]
      mu<-x[(n+1):(2*n)]
      phi<-x[(2*n+1)]
      
      res <- (y == 0) * ((1 + sqrt(phi * (2 * mu + phi))) / phi) + (y > 0) * (mu * dPIG((y > 0) * (y - 1), mu, 1 / phi) /
                                                                                (((y == 0) + y) * dPIG(y, mu, 1 / phi))) 
      - .5 * res
    }
    
    # Initial guess for the parameters through Yule-Walker estimation
    fit_yw <- ar.yw(x, aic = FALSE, order.max = 1)
    
    alpha1old <- fit_yw$ar
    alpha0old <- (1 - alpha1old) * mean(x)
    phiold <- (mean(x)^2 * d2B(qsi0)) / (var(x) - mean(x))
    
    alpha0yw <- alpha0old
    alpha1yw <- alpha1old
    phiyw <- phiold
    
    # Begin of the estimation process
    # Setting the tolerance - stopping criteria
    tol <- tol
    tolerance <- 1
    
    while(tolerance > tol){
      
      # previous mu
      muold <- NULL
      muold[1] <- alpha0old / (1 - alpha1old)  
      for (i in 1:(n-1)){
        muold[i+1] <- alpha0old + alpha1old * x[i]
      }  
      
      # expected likelihood function    
      q.function <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        
        # Current mu
        mu <- NULL
        mu[1] <- alpha0 / (1 - alpha1)
        
        for (i in 1:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i]
        }
        
        # Expected likelihood function
        sum(x*log(mu) - mu * eta_r(c(x, muold, phiold)) + D(phiold) + phiold*(eta_r(c(x, muold, phiold)) * qsi0 - 
                                                                                B(qsi0) + kappa_r(c(x, muold, phiold))))
      }
      
      # Gradient vector
      grad <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        
        # Current mu
        mu <- NULL
        mu[1] <- alpha0 / (1 - alpha1)
        
        for (i in 1:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i]
        }
        
        # grad functions
        grad1 <- function(par){
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum(x[1:n] / mu[1:n] - eta_r(c(x, muold, phiold))[1:n])
        }
        
        grad2 <- function(par){
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum(x[2:n] * x[1:(n-1)] / mu[2:n] - x[1:(n-1)] * eta_r(c(x, muold, phiold))[2:n])
        }
        
        c(grad1(c(x, mu, muold)), grad2(c(x, mu, muold)))
      }
      
      # Formulation of the constraints: ui %*% theta >= ci,
      # where theta = (alpha0, alpha1)
      ui <- diag(1,2)
      ci <- c(0, 0)
      fit <- try(constrOptim(theta = c(alpha0old, alpha1old), f = q.function, grad = grad, ui = ui, ci = ci,
                             method = method, control = list(fnscale = -1)))
      
      if(class(fit) != "try-error"){
        
        alpha0new <- fit$par[1]
        alpha1new <- fit$par[2]
        
        # Updating phi        
        phinew <- n / (2 * sum(eta_r(c(x, muold, phiold))/2 - 1 - kappa_r(c(x, muold, phiold))))
        
        # Tolerance criterion
        tolerance <- max(abs(c(alpha0new, alpha1new) - c(alpha0old, alpha1old)))
        
        # updating parameters values
        alpha0old <- alpha0new
        alpha1old <- alpha1new
        phiold <- phinew
      }
      
      if(class(fit) == "try-error"){
        tolerance = tol
        stop("Error in the estimation process!")
      }  
    }
    # end of the estimation process
    
    im <- function(z){
      
      x <- z[1:n]
      alpha0 <- z[n+1]
      alpha1 <- z[n+2]
      phi <- z[n+3]
      
      n <- length(x)
      
      # Auxiliary functions
      qsi0 <- -1/2
      B <- function(x) {-sqrt((-2 * x))}
      d2B <- function(x) {(2 * sqrt(2) * (-x)^(3 / 2))^(-1)}
      T <- function(x) {x}
      D <- function(x) {log(x) / 2}
      dD <- function(x) {1 / (2 * x)}
      d2D <- function(x) {-1 / (2 * x^2)}
      H <- function(x) {-log(2 * pi * x^3) / 2}
      G <- function(x) {-1/(2 * x)}
      
      #Conditional expectations: E(Z|F_(t-1))
      eta_t <- function(x){
        
        y <- x[1:n]
        mu <- x[(n+1):(2*n)]
        phi <- x[2*n+1]
        
        (y + 1) * dPIG(y + 1, mu, 1 / phi) / (mu * dPIG(y, mu, 1 / phi))
      }
      
      # E(Z^2|F_(t-1))
      gamma_t <- function(x){
        
        y <- x[1:n]
        mu <- x[(n+1):(2*n)]
        phi <- x[2*n+1]
        
        (y + 1) * (y + 2) * dPIG(y + 2, mu, 1 / phi) / (mu^2 * dPIG(y, mu, 1 / phi))
      }
      
      # E[g(Z)|F_(t-1)]
      kappa_t <- function(x){
        
        y <- x[1:n]
        mu <- x[(n+1):(2*n)]
        phi <- x[(2*n+1)]
        
        res <- (y == 0) * ((1 + sqrt(phi * (2 * mu + phi))) / phi) +
          (y > 0) * (mu * dPIG((y > 0) * (y - 1), mu, 1 / phi) /
                       (((y == 0) + y) * dPIG(y, mu, 1 / phi))) 
        - .5 * res
      }
      
      # estimated mu
      mu <- NULL
      mu[1] <- alpha0 / (1 - alpha1)
      
      for (i in 1:(n-1)){
        mu[i+1] <- alpha0 + alpha1 * x[i]
      }
      
      # First derivatives of mu_t
      # d_mu/d_ao
      dmut.a0 <- function(z){
        
        a1 <- z[1]
        a2 <- z[2]
        b1 <- z[3]
        
        dmut <- NULL
        dmut[1] <- 1 / (1 - a1 - a2 - b1)
        
        for(i in 2:n){
          dmut[i] <- 1 + b1 * dmut[i-1] 
        }
        dmut
      }
      
      # d_mu/d_a1
      dmut.a1 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        a2 <- z[3]
        b1 <- z[4]
        
        dmut <- NULL
        dmut[1] <- a0 / (1 - a1 - a2 - b1)^2
        
        for(i in 2:n){
          dmut[i] <- x[i-1] + b1 * dmut[i-1]
        }
        dmut
      }
      
      E_d2l_a0 <- sum((x/mu^2) * dmut.a0(c(alpha1, 0, 0))^2)
      E_d2l_a1 <- sum((x/mu^2) * dmut.a1(c(alpha0, alpha1, 0, 0))^2)
      E_d2l_a0_a1 <- sum((x/mu^2) * dmut.a0(c(alpha1, 0, 0)) * dmut.a1(c(alpha0, alpha1, 0, 0)))
      E_d2l_a0_phi <- 0
      E_d2l_a1_phi <- 0
      E_d2l_phi <- - n * d2D(phi)
      
      z1 <- dD(phi) + eta_t(c(x, mu, phi)) * qsi0 - B(qsi0) + kappa_t(c(x, mu, phi))
      z2_a0 <- (x/mu - eta_t(c(x, mu, phi))) * dmut.a0(c(alpha1, 0, 0))
      z2_a1 <- (x/mu - eta_t(c(x, mu, phi))) * dmut.a1(c(alpha0, alpha1, 0, 0))
      z3 <- (x/mu)^2 - 2 * (x/mu) * eta_t(c(x, mu, phi)) + gamma_t(c(x, mu, phi))
      
      m_a0 <- outer(as.numeric(z2_a0), as.numeric(z2_a0))
      m_a1 <- outer(as.numeric(z2_a1), as.numeric(z2_a1))
      m_a0_a1 <- outer(as.numeric(z2_a0), as.numeric(z2_a1))
      m_a0_phi <- outer(as.numeric(z2_a0), as.numeric(z1))
      m_a1_phi <- outer(as.numeric(z2_a1), as.numeric(z1))
      m_phi_phi <- outer(as.numeric(z1), as.numeric(z1))
      
      E_dl_a0_a0 <- sum(z3 * dmut.a0(c(alpha1, 0, 0))^2) + sum(m_a0) - sum(diag(m_a0))
      E_dl_a0_a1 <- sum(z3 * dmut.a0(c(alpha1, 0, 0)) * dmut.a1(c(alpha0, alpha1, 0, 0))) + sum(m_a0_a1) - sum(diag(m_a0_a1))
      E_dl_a1_a1 <- sum(z3 * dmut.a1(c(alpha0, alpha1, 0, 0))^2) + sum(m_a1) - sum(diag(m_a1))
      E_dl_a0_phi <- sum(m_a0_phi)
      E_dl_a1_phi <- sum(m_a1_phi)
      E_dl_phi_phi <- sum(m_phi_phi)
      
      I_a0_a0 <- E_d2l_a0 - E_dl_a0_a0
      I_a0_a1 <- E_d2l_a0_a1 - E_dl_a0_a1
      I_a1_a1 <- E_d2l_a1 - E_dl_a1_a1
      I_a0_phi <- E_d2l_a0_phi - E_dl_a0_phi
      I_a1_phi <- E_d2l_a1_phi - E_dl_a1_phi
      I_phi_phi <- E_d2l_phi - E_dl_phi_phi
      
      I <- matrix(
        c(I_a0_a0, I_a0_a1, I_a0_phi, 
          I_a0_a1, I_a1_a1, I_a1_phi,
          I_a0_phi, I_a1_phi, I_phi_phi),
        nrow = 3, ncol = 3)
      
      se <- round(sqrt(diag(solve(I))), 3)
      return(list("se" = se))
    }
    
    #Summary
    thetahat <- round(c(alpha0new, alpha1new, phinew), 3)
    se <- im(c(x, alpha0new, alpha1new, phinew))$se
    
    # Estimated mu
    muhat <- NULL
    muhat[1] <- alpha0new / (1 - alpha1new)
    
    for (i in 1:(n-1)){
      muhat[i+1] <- alpha0new + alpha1new * x[i]
    }
    
    loglik <- sum(log(dPIG(x, muhat, 1 / phinew)))
    # AIC: -2 * l(thetahat) + 2 * (p+q+1)
    aic <- -2 * loglik + (1 + 0 + 1) * 2
    # BIC: -2 * l(thetahat) + (p+q+1) * log(n)
    bic <- -2 * loglik + (1 + 0 + 1) * log(n)
    
    # Output
    return(list("Parameters" = thetahat,
                "SEs" = se,
                "Likelihood" = loglik,
                "AIC" = aic,
                "BIC" = bic))
  }# end of the pigingarch1 function
  
  # PIGINGARCH(2,0) case
  
  pigingarch2 <- function(x, method, tol){
    
    n <- length(x)
    
    # Auxiliary functions
    qsi0 <- -1/2
    B <- function(x) {-sqrt((-2 * x))}
    d2B <- function(x) {(2 * sqrt(2) * (-x)^(3/2))^(-1)}
    T <- function(x) {x}
    D <- function(x) {log(x) / 2}
    dD <- function(x) {1/(2 * x)}
    H <- function(x) {-log(2 * pi * x^3) / 2}
    G <- function(x) {-1/(2 * x)}
    
    # Conditional expectations
    eta_r <- function(x){
      
      y <- x[1:n]
      mu <- x[(n+1):(2*n)]
      phi <- x[2*n+1]
      
      (y + 1) * dPIG(y + 1, mu, 1 / phi) / (mu * dPIG(y, mu, 1 / phi))
    }
    
    kappa_r <- function(x){
      
      y<-x[1:n]
      mu<-x[(n+1):(2*n)]
      phi<-x[2*n+1]
      
      res <- (y == 0) * ((1 + sqrt(phi * (2 * mu + phi))) / phi) +
        (y > 0) * (mu * dPIG((y > 0) * (y - 1), mu, 1 / phi) /
                     (((y == 0) + y) * dPIG(y, mu, 1 / phi))) 
      -.5 * res
    }
    
    # Initial guess for the parameters through Yule-Walker estimation
    fit_yw <- ar.yw(x, aic = FALSE, order.max = 2)
    
    alpha1old <- fit_yw$ar[1]
    alpha2old <- fit_yw$ar[2]
    
    alpha0old <- (1 - alpha1old - alpha2old) * mean(x)
    phiold <- (mean(x)^2 * d2B(qsi0)) / (var(x) - mean(x))
    
    # Begin of the estimation process
    tolerance <- 1
    while(tolerance > tol){
      
      # previous mu
      muold <- NULL
      muold[1] <- alpha0old / (1 - alpha1old - alpha2old)
      muold[2] <- alpha0old + alpha1old * x[1]
      
      for (i in 2:(n-1)){
        muold[i+1] <- alpha0old + alpha1old * x[i] + alpha2old * x[i-1]
      }  
      
      # expected likelihood function    
      q.function <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        alpha2 <- z[3]
        
        # Current mu
        mu <- NULL
        mu[1] <- alpha0 / (1 - alpha1 - alpha2)
        mu[2] <- alpha0 + alpha1 * x[1] 
        
        for (i in 2:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i] + alpha2 * x[i-1]
        }
        
        # Expected likelihood function
        sum(x * log(mu) - mu * eta_r(c(x, muold, phiold)) + 
              D(phiold) + phiold * (eta_r(c(x, muold, phiold)) * qsi0 - 
                                      B(qsi0) + kappa_r(c(x, muold, phiold))))
      }
      
      # Gradient vector
      grad <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        alpha2 <- z[3]
        
        # Current mu
        mu <- NULL
        mu[1] <- alpha0 / (1 - alpha1 - alpha2)
        mu[2] <- alpha0 + alpha1 * x[1]
        
        for (i in 2:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i] + alpha2 * x[i-1]
        }
        
        # grad functions
        grad1 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum(x[1:n] / mu[1:n] - eta_r(c(x, muold, phiold))[1:n])
        }
        
        grad2 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum(x[2:n] * x[1:(n-1)] / mu[2:n] - x[1:(n-1)] * eta_r(c(x, muold, phiold))[2:n])
        }
        
        grad3 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum(x[3:n] * x[1:(n-2)] / mu[3:n] - x[1:(n-2)] * eta_r(c(x, muold, phiold))[3:n])
        }
        
        c(grad1(c(x, mu, muold)), grad2(c(x, mu, muold)), grad3(c(x, mu, muold)))
      }
      
      # Formulation of the constraints: ui %*% theta >= ci,
      # where theta = (alpha0, alpha1, alpha2)
      ui <- diag(1,3)
      ci <- c(0, 0, 0)
      fit <- try(constrOptim(theta = c(alpha0old, alpha1old, alpha2old), f = q.function,
                             grad = grad, ui = ui, ci = ci,
                             method = method,
                             control = list(fnscale = -1)))
      
      if(class(fit) != "try-error"){
        
        alpha0new <- fit$par[1]
        alpha1new <- fit$par[2]
        alpha2new <- fit$par[3]
        
        # Updating phi    
        phinew <- n / (2 * sum(eta_r(c(x, muold, phiold))/2 - 1 - kappa_r(c(x, muold, phiold))))
        
        # tolerance criterion
        tolerance <- max(abs(c(alpha0new, alpha1new, alpha2new) - c(alpha0old, alpha1old, alpha2old)))
        
        # updating parameters values
        alpha0old <- alpha0new
        alpha1old <- alpha1new
        alpha2old <- alpha2new
        phiold <- phinew
      }
      
      if((class(fit) == "try-error")){
        tolerance = tol
        stop("Error in the estimation process!")
      }
    }
    # end of the estimation process
    
    #Summary
    thetahat <- c(alpha0new, alpha1new, alpha2new, phinew)
    
    im <- function(z){
      
      x <- z[1:n]
      alpha0 <- z[n+1]
      alpha1 <- z[n+2]
      alpha2 <- z[n+3]
      phi <- z[n+4]
      
      n <- length(x)
      
      # Auxiliary functions
      qsi0 <- -1/2
      B <- function(x) {-sqrt((-2 * x))}
      d2B <- function(x) {(2 * sqrt(2) * (-x)^(3 / 2))^(-1)}
      T <- function(x) {x}
      D <- function(x) {log(x) / 2}
      dD <- function(x) {1 / (2 * x)}
      d2D <- function(x) {-1 / (2 * x^2)}
      H <- function(x) {-log(2 * pi * x^3) / 2}
      G <- function(x) {-1 / (2 * x)}
      
      # Conditional expectations: E(Z|F_(t-1))
      eta_t <- function(x){
        
        y <- x[1:n]
        mu <- x[(n+1):(2*n)]
        phi <- x[2*n+1]
        
        (y + 1) * dPIG(y + 1, mu, 1 / phi) / (mu * dPIG(y, mu, 1 / phi))
      }
      
      # E(Z^2|F_(t-1))
      gamma_t <- function(x){
        
        y <- x[1:n]
        mu <- x[(n+1):(2*n)]
        phi <- x[2*n+1]
        
        (y + 1) * (y + 2) * dPIG(y + 2, mu, 1 / phi) / (mu^2 * dPIG(y, mu, 1 / phi))
      }
      
      # E[g(Z)|F_(t-1)]
      kappa_t <- function(x){
        
        y <- x[1:n]
        mu <- x[(n+1):(2*n)]
        phi <- x[(2*n+1)]
        
        res <- (y == 0) * ((1 + sqrt(phi * (2 * mu + phi))) / phi) +
          (y > 0) * (mu * dPIG((y > 0) * (y - 1), mu, 1 / phi) /
                       (((y == 0) + y) * dPIG(y, mu, 1 / phi))) 
        - .5 * res
      }
      
      # Estimated mu
      mu <- NULL
      mu[1] <- alpha0 / (1 - alpha1 - alpha2)
      mu[2] <- alpha0 + alpha1 * x[1] 
      
      for (i in 2:(n-1)){
        mu[i+1] <- alpha0 + alpha1 * x[i] + alpha2 * x[i-1]
      }
      
      # First derivatives of mu_t
      # d_mu/d_ao
      dmut.a0 <- function(z){
        
        a1 <- z[1]
        a2 <- z[2]
        b1 <- z[3]
        
        dmut <- NULL
        dmut[1] <- 1 / (1 - a1 - a2 - b1)
        
        for(i in 2:n){
          dmut[i] <- 1 + b1 * dmut[i-1] 
        }
        dmut
      }
      # d_mu/d_a1
      dmut.a1 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        a2 <- z[3]
        b1 <- z[4]
        
        dmut <- NULL
        dmut[1] <- a0 / (1 - a1 - a2 - b1)^2
        
        for(i in 2:n){
          dmut[i] <- x[i-1] + b1 * dmut[i-1]
        }
        dmut
      }
      # d_mu/d_a2
      dmut.a2 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        a2 <- z[3]
        
        dmut <- NULL
        dmut[1] <- a0 / (1 - a1 - a2)^2
        dmut[2] <- a0 / (1 - a1 - a2)^2
        
        for(i in 3:n){
          dmut[i] <- x[i-2] 
        }
        dmut
      }
      
      E_d2l_a0 <- sum((x/mu^2) * dmut.a0(c(alpha1, alpha2, 0))^2)
      E_d2l_a1 <- sum((x/mu^2) * dmut.a1(c(alpha0, alpha1, alpha2, 0))^2)
      E_d2l_a2 <- sum((x/mu^2) * dmut.a2(c(alpha0, alpha1, alpha2))^2)
      E_d2l_phi <- - n * d2D(phi)
      E_d2l_a0_a1 <- sum((x/mu^2) * dmut.a0(c(alpha1, alpha2, 0)) * dmut.a1(c(alpha0, alpha1, alpha2, 0)))
      E_d2l_a0_a2 <- sum((x/mu^2) * dmut.a0(c(alpha1, alpha2, 0)) * dmut.a2(c(alpha0, alpha1, alpha2)))
      E_d2l_a1_a2 <- sum((x/mu^2) * dmut.a1(c(alpha0, alpha1, alpha2, 0)) * dmut.a2(c(alpha0, alpha1, alpha2)))
      E_d2l_a0_phi <- 0
      E_d2l_a1_phi <- 0
      E_d2l_a2_phi <- 0
      
      z1 <- dD(phi) + eta_t(c(x, mu, phi)) * qsi0 - B(qsi0) + kappa_t(c(x, mu, phi))
      z2_a0 <- (x/mu - eta_t(c(x, mu, phi))) * dmut.a0(c(alpha1, alpha2, 0))
      z2_a1 <- (x/mu - eta_t(c(x, mu, phi))) * dmut.a1(c(alpha0, alpha1, alpha2, 0))
      z2_a2 <- (x/mu - eta_t(c(x, mu, phi))) * dmut.a2(c(alpha0, alpha1, alpha2))
      z3 <- (x/mu)^2 - 2 * (x/mu) * eta_t(c(x, mu, phi)) + gamma_t(c(x, mu, phi))
      
      m_a0 <- outer(as.numeric(z2_a0), as.numeric(z2_a0))
      m_a1 <- outer(as.numeric(z2_a1), as.numeric(z2_a1))
      m_a2 <- outer(as.numeric(z2_a2), as.numeric(z2_a2))
      m_a0_a1 <- outer(as.numeric(z2_a0), as.numeric(z2_a1))
      m_a0_a2 <- outer(as.numeric(z2_a0), as.numeric(z2_a2))
      m_a1_a2 <- outer(as.numeric(z2_a1), as.numeric(z2_a2))
      m_a0_phi <- outer(as.numeric(z2_a0), as.numeric(z1))
      m_a1_phi <- outer(as.numeric(z2_a1), as.numeric(z1))
      m_a2_phi <- outer(as.numeric(z2_a2), as.numeric(z1))
      m_phi_phi <- outer(as.numeric(z1), as.numeric(z1))
      
      E_dl_a0_a0 <- sum(z3 * dmut.a0(c(alpha1, alpha2, 0))^2) + sum(m_a0) - sum(diag(m_a0))
      E_dl_a0_a1 <- sum(z3 * dmut.a0(c(alpha1, alpha2, 0)) * dmut.a1(c(alpha0, alpha1, alpha2, 0))) + sum(m_a0_a1) - sum(diag(m_a0_a1))
      E_dl_a0_a2 <- sum(z3 * dmut.a0(c(alpha1, alpha2, 0)) * dmut.a2(c(alpha0, alpha1, alpha2))) + sum(m_a0_a2) - sum(diag(m_a0_a2))
      E_dl_a1_a1 <- sum(z3 * dmut.a1(c(alpha0, alpha1, alpha2, 0))^2) + sum(m_a1) - sum(diag(m_a1))
      E_dl_a1_a2 <- sum(z3 * dmut.a1(c(alpha0, alpha1, alpha2, 0)) * dmut.a2(c(alpha0, alpha1, alpha2))) + sum(m_a1_a2) - sum(diag(m_a1_a2))
      E_dl_a2_a2 <- sum(z3 * dmut.a2(c(alpha0, alpha1, alpha2))^2) + sum(m_a2) - sum(diag(m_a2))
      E_dl_a0_phi <- sum(m_a0_phi)
      E_dl_a1_phi <- sum(m_a1_phi)
      E_dl_a2_phi <- sum(m_a2_phi)
      E_dl_phi_phi <- sum(m_phi_phi)
      
      I_a0_a0 <- E_d2l_a0 - E_dl_a0_a0
      I_a0_a1 <- E_d2l_a0_a1 - E_dl_a0_a1
      I_a0_a2 <- E_d2l_a0_a2 - E_dl_a0_a2
      I_a0_phi <- E_d2l_a0_phi - E_dl_a0_phi
      I_a1_a1 <- E_d2l_a1 - E_dl_a1_a1
      I_a1_a2 <- E_d2l_a1_a2 - E_dl_a1_a2
      I_a1_phi <- E_d2l_a1_phi - E_dl_a1_phi
      I_a2_a2 <- E_d2l_a2 - E_dl_a2_a2
      I_a2_phi <- E_d2l_a2_phi - E_dl_a2_phi
      I_phi_phi <- E_d2l_phi - E_dl_phi_phi
      
      I <- matrix(
        c(I_a0_a0, I_a0_a1, I_a0_a2, I_a0_phi, 
          I_a0_a1, I_a1_a1, I_a1_a2, I_a1_phi,
          I_a0_a2, I_a1_a2, I_a2_a2, I_a2_phi,
          I_a0_phi, I_a1_phi, I_a2_phi, I_phi_phi),
        nrow = 4, ncol = 4)
      
      se <- round(sqrt(diag(solve(I))), 3)
      return(list("se" = se))
    }
    
    # estimated mu
    muhat <- NULL
    muhat[1] <- alpha0new / (1 - alpha1new - alpha2new)
    muhat[2] <- alpha0new + alpha1new * x[1]
    
    for (i in 2:(n-1)){
      muhat[i+1] <- alpha0new + alpha1new * x[i] + alpha2new * x[i-1]
    }
    
    loglik <- sum(log(dPIG(x, muhat, 1 / phinew)))
    
    se <- im(c(x, alpha0new, alpha1new, alpha2new, phinew))$se
    
    # AIC: -2 * l(thetahat) + 2 * (p+q+1)
    aic <- -2 * loglik + 2 * (2 + 0 + 1)
    # BIC: -2 * l(thetahat) + (p+q+1) * log(n)
    bic <- -2 * loglik + (2 + 0 + 1) * log(n)
    
    # Output
    return(list("Parameters" = thetahat,
                "SEs" = se,
                "Likelihood" = loglik,
                "AIC" = aic,
                "BIC" = bic))
  }# end of the pigingarch2 function
  
  # PIGINGARCH(1,1) case
  
  pigingarch11 <- function(x, method, tol){
    
    n <- length(x)
    
    # Auxiliary functions
    qsi0 <- -1/2
    B <- function(x) {-sqrt((-2 * x))}
    d2B <- function(x) {(2 * sqrt(2) * (-x)^(3/2))^(-1)}
    T <- function(x) {x}
    D <- function(x) {log(x) / 2}
    dD <- function(x) {1/(2 * x)}
    H <- function(x) {-log(2 * pi * x^3) / 2}
    G <- function(x) {-1/(2 * x)}
    
    # Conditional expectations
    eta_r <- function(x){
      
      y <- x[1:n]
      mu <- x[(n+1):(2*n)]
      phi <- x[2*n+1]
      
      (y + 1) * dPIG(y + 1, mu, 1 / phi) / (mu * dPIG(y, mu, 1 / phi))
    }
    
    kappa_r <- function(x){
      
      y<-x[1:n]
      mu<-x[(n+1):(2*n)]
      phi<-x[2*n+1]
      
      res <- (y == 0) * ((1 + sqrt(phi * (2 * mu + phi))) / phi) +
        (y > 0) * (mu * dPIG((y > 0) * (y - 1), mu, 1 / phi) /
                     (((y == 0) + y) * dPIG(y, mu, 1 / phi))) 
      -.5 * res
    }
    
    # Conditional least square (CLS) estimation
    fitcls <- arma(x, order = c(1, 1), lag = NULL, coef = NULL,
                   include.intercept = TRUE, series = NULL)
    
    alpha0cls <- fitcls$coef[3]
    alpha1cls <- fitcls$coef[1]
    beta1cls <- fitcls$coef[2]
    
    mu <- NULL
    mu[1] <- alpha0cls / (1 - alpha1cls - beta1cls)
    
    for (i in 1:(n-1)){
      mu[i+1] <- alpha0cls + alpha1cls * x[i] + beta1cls * mu[i]
    }
    
    phicls <- (sum(mu^4)) / (sum(x^2 * mu^2 - mu^3 - mu^4))
    
    ig <- c(alpha0cls, alpha1cls, beta1cls, phicls)
    ig[which(ig < 0 | ig > 10)] <- 0.1
    
    alpha0old <- ig[1]
    alpha1old <- ig[2]
    beta1old <- ig[3]
    phiold <- ig[4]
    
    # Begin of the estimation process
    tolerance <- 1
    while(tolerance > tol){
      
      # previous mu
      muold <- NULL
      muold[1] <- alpha0old / (1 - alpha1old - beta1old)  
      
      for (i in 1:(n-1)){
        muold[i+1] <- alpha0old + alpha1old * x[i] + beta1old * muold[i]
      }  
      
      # expected likelihood function    
      q.function <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        beta1 <- z[3]
        
        # Current mu
        mu <- NULL
        mu[1] <- alpha0 / (1 - alpha1 - beta1)
        
        for (i in 1:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i] + beta1 * mu[i]
        } 
        
        # Expected likelihood function
        sum(x * log(mu) - mu * eta_r(c(x, muold, phiold)) + 
              D(phiold) + phiold * (eta_r(c(x, muold, phiold)) * qsi0 - 
                                      B(qsi0) + kappa_r(c(x, muold, phiold))))
      }
      
      # First derivative of mu_t with respect to alpha_0
      dmut.a0 <- function(z){
        
        a1 <- z[1]
        a2 <- z[2]
        b1 <- z[3]
        
        dmut <- array(0, c(n,1))
        dmut[1] <- 1 / (1 - a1 - a2 - b1)^2
        dmut[2:n] <- 1 + b1 * dmut[1:(n-1)]
        dmut
      }
      
      # First derivative of mu_t with respect to alpha_1
      dmut.a1 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        a2 <- z[3]
        b1 <- z[4]
        
        dmut <- array(0, c(n,1))
        dmut[1] <- a0 / (1 - a1 - a2 - b1)^2
        dmut[2:n] <- x[1:(n-1)] + b1 * dmut[1:(n-1)]
        dmut
      }
      
      # First derivative of mu_t with respect to beta1
      dmut.b1 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        b1 <- z[3]
        
        dmut <- array(0, c(n,1))
        dmut[1] <- a0 / (1 - a1 - b1)^2
        dmut[2:n] <- mu[1:(n-1)] + b1 * dmut[1:(n-1)]
        dmut
      }
      
      # Gradient vector
      grad <- function(z){
        
        alpha0 <- z[1]
        alpha1 <- z[2]
        beta1 <- z[3]
        
        # Current mu
        mu <- NULL
        mu[1] <- alpha0 / (1 - alpha1 - beta1)
        
        for (i in 1:(n-1)){
          mu[i+1] <- alpha0 + alpha1 * x[i] + beta1 * mu[i]
        }
        # grad functions (alpha0)
        grad1 <- function(par){ 
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum((x[1:n] / mu[1:n] - eta_r(c(x, muold, phiold))[1:n]) * dmut.a0(c(alpha1, 0, beta1))[1:n])
        }
        grad2 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum((x[1:n] / mu[1:n] - eta_r(c(x, muold, phiold))[1:n]) * dmut.a1(c(alpha0, alpha1, 0, beta1))[1:n])
        }
        grad3 <- function(par){
          
          x <- par[1:n]
          mu <- par[(n+1):(2*n)]
          muold <- par[(2*n+1):(3*n)]
          
          sum((x[1:n] / mu[1:n] - eta_r(c(x, muold, phiold))[1:n]) * dmut.b1(c(alpha0, alpha1, beta1))[1:n])
        }
        
        c(grad1(c(x, mu, muold)), grad2(c(x, mu, muold)), grad3(c(x, mu, muold)))
      }
      
      # Formulation of the constraints: ui %*% theta >= ci,
      # where theta = (alpha0, alpha1, beta1)
      ui <- diag(1,3)
      ci <- c(0, 0, 0)
      fit <- try(constrOptim(theta = c(alpha0old, alpha1old, beta1old), f = q.function,
                             grad = grad, ui = ui, ci = ci,
                             method = method,
                             control = list(fnscale = -1)))
      
      if(class(fit) != "try-error"){ 
        
        alpha0new <- fit$par[1]
        alpha1new <- fit$par[2]
        beta1new <- fit$par[3]
        
        # Updating phi      
        phinew <- n / (2*sum(eta_r(c(x, muold, phiold))/2 - 1 - kappa_r(c(x, muold, phiold))))
        
        # tolerance criterion
        tolerance <- max(abs(c(alpha0new, alpha1new, beta1new) - c(alpha0old, alpha1old, beta1old)))
        
        # updating parameters values
        alpha0old <- alpha0new
        alpha1old <- alpha1new
        beta1old <- beta1new
        phiold <- phinew
      }
      
      if((class(fit) == "try-error")){
        tolerance = tol
        stop("Error in the estimation process!")
      }
    }# end of the estimation process
    
    im <- function(z){
      
      x <- z[1:n]
      alpha0 <- z[n+1]
      alpha1 <- z[n+2]
      beta1 <- z[n+3]
      phi <- z[n+4]
      
      n <- length(x)
      
      # Auxiliary functions
      qsi0 <- -1/2
      B <- function(x) {-sqrt((-2 * x))}
      d2B <- function(x) {(2 * sqrt(2) * (-x)^(3 / 2))^(-1)}
      T <- function(x) {x}
      D <- function(x) {log(x) / 2}
      dD <- function(x) {1 / (2 * x)}
      d2D <- function(x) {-1 / (2 * x^2)}
      H <- function(x) {-log(2 * pi * x^3) / 2}
      G <- function(x) {-1 / (2 * x)}
      
      #Conditional expectations: E(Z|F_(t-1))
      eta_t <- function(x){
        
        y <- x[1:n]
        mu <- x[(n+1):(2*n)]
        phi <- x[2*n+1]
        
        (y + 1) * dPIG(y + 1, mu, 1 / phi) / (mu * dPIG(y, mu, 1 / phi))
      }
      # E(Z^2|F_(t-1))
      gamma_t <- function(x){
        
        y <- x[1:n]
        mu <- x[(n+1):(2*n)]
        phi <- x[2*n+1]
        
        (y + 1) * (y + 2) * dPIG(y + 2, mu, 1 / phi) / (mu^2 * dPIG(y, mu, 1 / phi))
      }
      # E[g(Z)|F_(t-1)]
      kappa_t <- function(x){
        
        y<-x[1:n]
        mu<-x[(n+1):(2*n)]
        phi<-x[(2*n+1)]
        
        res <- (y == 0) * ((1 + sqrt(phi * (2 * mu + phi))) / phi) +
          (y > 0) * (mu * dPIG((y > 0) * (y - 1), mu, 1 / phi) /
                       (((y == 0) + y) * dPIG(y, mu, 1 / phi))) 
        - .5 * res
      }
      
      n <- length(x)
      
      # Estimated mu
      mu <- NULL
      mu[1] <- alpha0 / (1 - alpha1 - beta1)
      
      for (i in 1:(n-1)){
        mu[i+1] <- alpha0 + alpha1 * x[i] + beta1 * mu[i]
      }
      
      # First derivatives
      # d_mu/d_a0
      dmut.a0 <- function(z){
        
        a1 <- z[1]
        b1 <- z[2]
        
        dmut <- NULL
        dmut[1] <- 1 / (1 - a1 - b1)
        
        for(i in 2:n){
          dmut[i] <- 1 + b1 * dmut[i-1]
        }
        dmut
      }
      # d_mu/d_a1
      dmut.a1 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        b1 <- z[3]
        
        dmut <- NULL
        dmut[1] <- a0 / (1 - a1 - b1)^2
        
        for(i in 2:n){
          dmut[i] <- x[i-1] + b1 * dmut[i-1]
        }
        dmut
      }
      # d_mu/d_b1
      dmut.b1 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        b1 <- z[3]
        
        dmut <- NULL
        dmut[1] <- a0 / (1 - a1 - b1)^2
        
        for(i in 2:n){
          dmut[i] <- mu[i-1] + b1 * dmut[i-1]
        }
        dmut
      }
      
      #Second derivatives: d2_mu/d_a0^2
      d2mut_a0_a0 <- 0
      
      # d2_mu/da0_da1
      d2mut_a0_a1 <- function(z){
        
        a1 <- z[1]
        b1 <- z[2]
        
        d2mut <- NULL
        d2mut[1] <- 1 / (1 - a1 - b1)^2
        
        for(i in 2:n){
          d2mut[i] <- b1 * d2mut[i-1]
        }
        d2mut
      }
      # d2_mu/da0_db1
      d2mut_a0_b1 <- function(z){
        
        a1 <- z[1]
        b1 <- z[2]
        
        d2mut <- NULL
        d2mut[1] <- 1 / (1 - a1 - b1)^2
        
        for(i in 2:n){
          d2mut[i] <- dmut.a0(c(a1, b1))[i-1] + b1 * d2mut[i-1]
        }
        d2mut
      }
      # d2_mu/da1^2
      d2mut_a1_a1 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        b1 <- z[3]
        
        d2mut <- NULL
        d2mut[1] <- (2*a0) / (1 - a1 - b1)^3
        
        for(i in 2:n){
          d2mut[i] <- b1 * d2mut[i-1] 
        }
        d2mut
      }
      # d2_mu/da1_db1
      d2mut_a1_b1 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        b1 <- z[3]
        
        d2mut <- NULL
        d2mut[1] <- (2*a0) / (1 - a1 - b1)^3
        
        for(i in 2:n){
          d2mut[i] <- dmut.a1(c(alpha0, alpha1, beta1))[i-1] + b1 * d2mut[i-1] 
        }
        d2mut
      }
      # d2_mu/db1^2
      d2mut_b1_b1 <- function(z){
        
        a0 <- z[1]
        a1 <- z[2]
        b1 <- z[3]
        
        d2mut <- NULL
        d2mut[1] <- (2*a0) / (1 - a1 - b1)^3
        
        for(i in 2:n){
          d2mut[i] <- 2 * dmut.b1(c(alpha0, alpha1, beta1))[i-1] + b1 * d2mut[i-1] 
        }
        d2mut
      }
      # d2_mu/da0_dphi
      d2mut_a0_phi <- 0
      # d2_mu/da1_dphi
      d2mut_a1_phi <- 0
      # d2_mu/db1_dphi
      d2mut_b1_phi <- 0
      # d2_mu/dphi_dphi
      d2mut_phi_phi <- 0
      
      E_d2l_a0 <- sum((x/mu^2) * dmut.a0(c(alpha1, beta1))^2)
      E_d2l_a1 <- sum((x/mu^2) * dmut.a1(c(alpha0, alpha1, beta1))^2 -
                        (x/mu - eta_t(c(x, mu, phi))) * d2mut_a1_a1(c(alpha0, alpha1, beta1)))
      E_d2l_b1 <- sum((x/mu^2) * dmut.b1(c(alpha0, alpha1, beta1))^2 -
                        (x/mu - eta_t(c(x, mu, phi))) * d2mut_b1_b1(c(alpha0, alpha1, beta1)))
      E_d2l_phi <- - n * d2D(phi)
      E_d2l_a0_a1 <- sum((x/mu^2) * dmut.a0(c(alpha1, beta1)) * dmut.a1(c(alpha0, alpha1, beta1)) -
                           (x/mu - eta_t(c(x, mu, phi))) * d2mut_a0_a1(c(alpha1, beta1)))
      E_d2l_a0_b1 <- sum((x/mu^2) * dmut.a0(c(alpha1, beta1)) * dmut.b1(c(alpha0, alpha1, beta1)) -
                           (x/mu - eta_t(c(x, mu, phi))) * d2mut_a0_b1(c(alpha1, beta1)))
      E_d2l_a1_b1 <- sum((x/mu^2) * dmut.a1(c(alpha0, alpha1, beta1)) * dmut.b1(c(alpha0, alpha1, beta1)) -
                           (x/mu - eta_t(c(x, mu, phi))) * d2mut_a1_b1(c(alpha0, alpha1, beta1)))
      E_d2l_a0_phi <- 0
      E_d2l_a1_phi <- 0
      E_d2l_b1_phi <- 0
      
      z1 <- dD(phi) + eta_t(c(x, mu, phi)) * qsi0 - B(qsi0) + kappa_t(c(x, mu, phi))
      z2_a0 <- (x/mu - eta_t(c(x, mu, phi))) * dmut.a0(c(alpha1, beta1))
      z2_a1 <- (x/mu - eta_t(c(x, mu, phi))) * dmut.a1(c(alpha0, alpha1, beta1))
      z2_b1 <- (x/mu - eta_t(c(x, mu, phi))) * dmut.b1(c(alpha0, alpha1, beta1))
      z3 <- (x/mu)^2 - 2 * (x/mu) * eta_t(c(x, mu, phi)) + gamma_t(c(x, mu, phi))
      
      m_a0 <- outer(as.numeric(z2_a0), as.numeric(z2_a0))
      m_a1 <- outer(as.numeric(z2_a1), as.numeric(z2_a1))
      m_b1 <- outer(as.numeric(z2_b1), as.numeric(z2_b1))
      m_a0_a1 <- outer(as.numeric(z2_a0), as.numeric(z2_a1))
      m_a0_b1 <- outer(as.numeric(z2_a0), as.numeric(z2_b1))
      m_a1_b1 <- outer(as.numeric(z2_a1), as.numeric(z2_b1))
      m_a0_phi <- outer(as.numeric(z2_a0), as.numeric(z1))
      m_a1_phi <- outer(as.numeric(z2_a1), as.numeric(z1))
      m_b1_phi <- outer(as.numeric(z2_b1), as.numeric(z1))
      m_phi_phi <- outer(as.numeric(z1), as.numeric(z1))
      
      E_dl_a0_a0 <- sum(z3 * dmut.a0(c(alpha1, beta1))^2) + sum(m_a0) - sum(diag(m_a0))
      E_dl_a0_a1 <- sum(z3 * dmut.a0(c(alpha1, beta1)) * dmut.a1(c(alpha0, alpha1, beta1))) + sum(m_a0_a1) - sum(diag(m_a0_a1))
      E_dl_a0_b1 <- sum(z3 * dmut.a0(c(alpha1, beta1)) * dmut.b1(c(alpha0, alpha1, beta1))) + sum(m_a0_b1) - sum(diag(m_a0_b1))
      E_dl_a1_a1 <- sum(z3 * dmut.a1(c(alpha0, alpha1, beta1))^2) + sum(m_a1) - sum(diag(m_a1))
      E_dl_a1_b1 <- sum(z3 * dmut.a1(c(alpha0, alpha1, beta1)) * dmut.b1(c(alpha0, alpha1, beta1))) + sum(m_a1_b1) - sum(diag(m_a1_b1))
      E_dl_b1_b1 <- sum(z3 * dmut.b1(c(alpha0, alpha1, beta1))^2) + sum(m_b1) - sum(diag(m_b1))
      E_dl_a0_phi <- sum(m_a0_phi)
      E_dl_a1_phi <- sum(m_a1_phi)
      E_dl_b1_phi <- sum(m_b1_phi)
      E_dl_phi_phi <- sum(m_phi_phi)
      
      I_a0_a0 <- E_d2l_a0 - E_dl_a0_a0
      I_a0_a1 <- E_d2l_a0_a1 - E_dl_a0_a1
      I_a0_b1 <- E_d2l_a0_b1 - E_dl_a0_b1
      I_a0_phi <- E_d2l_a0_phi - E_dl_a0_phi
      I_a1_a1 <- E_d2l_a1 - E_dl_a1_a1
      I_a1_b1 <- E_d2l_a1_b1 - E_dl_a1_b1
      I_a1_phi <- E_d2l_a1_phi - E_dl_a1_phi
      I_b1_b1 <- E_d2l_b1 - E_dl_b1_b1
      I_b1_phi <- E_d2l_b1_phi - E_dl_b1_phi
      I_phi_phi <- E_d2l_phi - E_dl_phi_phi
      
      I <- matrix(
        c(I_a0_a0, I_a0_a1, I_a0_b1, I_a0_phi, 
          I_a0_a1, I_a1_a1, I_a1_b1, I_a1_phi,
          I_a0_b1, I_a1_b1, I_b1_b1, I_b1_phi,
          I_a0_phi, I_a1_phi, I_b1_phi, I_phi_phi),
        nrow = 4, ncol = 4)
      
      se <- round(sqrt(diag(solve(I))), 3)
      return(list("se" = se))
    }
    
    #Summary
    thetahat <- c(alpha0new, alpha1new, beta1new, phinew)
    se <- im(c(x, alpha0new, alpha1new, beta1new, phinew))$se
    
    # estimated mu
    muhat <- NULL
    muhat[1] <- alpha0new / (1 - alpha1new - beta1new)
    
    for (i in 1:(n-1)){
      muhat[i+1] <- alpha0new + alpha1new * x[i] + beta1new * muhat[i]
    }
    
    loglik <- sum(log(dPIG(x, muhat, 1/phinew)))
    
    # AIC: -2 * l(thetahat) + (p+q+1) * 2
    aic <- -2 * loglik + (1 + 1 + 1) * 2
    # BIC: -2 * l(thetahat) + (p+q+1) * log(n)
    bic <- -2 * loglik + (1 + 1 + 1) * log(n)
    
    # Output
    return(list("Parameters" = thetahat,
                "SEs" = se,
                "Likelihood" = loglik,
                "AIC" = aic,
                "BIC" = bic))
  }# end of the pigingarch11 function
  
  order <- as.character(order)
  switch(type, 
         "PIG" = switch(order,
                        "1" = pigingarch1(x, method, tol = tol),
                        "2" = pigingarch2(x, method, tol = tol),
                        "11" = pigingarch11(x, method, tol = tol),
                        stop("Order not expected!")),
         "NB" = switch(order,
                       "1" = nbingarch1(x, method, tol = tol),
                       "2" = nbingarch2(x, method, tol = tol),
                       "11" = nbingarch11(x, method, tol = tol),
                       stop("Order not expected!")),
         stop("Type not expected!")
  )
  
}# end of the fit.mpingarch function


