#' Loglikelihood for diversity-independent diversification model
#'
#' This function computes loglikelihood of a diversity-independent
#' diversification model for a given set of branching times and parameter
#' values.
#'
#'
#' @param pars1 Vector of parameters:
#' \cr \cr \code{pars1[1]} corresponds to
#' lambda0 (speciation rate)
#' \cr \code{pars1[2]} corresponds to mu0 (extinction
#' rate)
#' \cr \code{pars1[3]} corresponds to lambda1 (decline parameter in
#' speciation rate) or K in diversity-dependence-like models
#' \cr \code{pars1[4]} corresponds to mu1 (decline parameter in extinction rate)
#' @param pars2 Vector of model settings:
#' \cr \cr \code{pars2[1]} sets the
#' model of time-dependence:
#' \cr - \code{pars2[1] == 0} no time dependence
#' \cr - \code{pars2[1] == 1} speciation and/or extinction rate is exponentially
#' declining with time
#' \cr - \code{pars2[1] == 2} stepwise decline in
#' speciation rate as in diversity-dependence without extinction
#' \cr - \code{pars2[1] == 3} decline in speciation rate following deterministic
#' logistic equation for ddmodel = 1
#' \cr - \code{pars2[1] == 4} decline in
#' speciation rate such that the expected number of species matches with that
#' of ddmodel = 1 with the same mu
#' \cr \cr \code{pars2[2]} sets the
#' conditioning:
#' \cr - \code{pars[2] == 0} conditioning on stem or crown age
#' \cr - \code{pars[2] == 1} conditioning on stem or crown age and
#' non-extinction of the phylogeny
#' \cr - \code{pars[2] == 2} conditioning on
#' stem or crown age and on the total number of extant taxa (including missing
#' species)
#' \cr - \code{pars[2] == 3} conditioning on the total number of
#' extant taxa (including missing species)
#' \cr \cr \code{pars2[3]} sets whether
#' the likelihood is for the branching times (0) or the phylogeny (1)
#' \cr \cr \code{pars2[4]} sets whether the parameters and likelihood should be 
#' shown
#' on screen (1) or not (0)
#' \cr \cr \code{pars2[5]} sets whether the first data
#' point is stem age (1) or crown age (2)
#' @param brts A set of branching times of a phylogeny, all positive
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @param methode The method used to solve the master equation, default is
#' 'odeint::runge_kutta_cash_karp54'.
#' @return The loglikelihood
#' @author Rampal S. Etienne, Bart Haegeman & Cesar Martinez
#' @seealso \code{\link{bd_ML}}
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' bd_loglik(pars1 = c(0.5,0.1), pars2 = c(0,1,1,0,2), brts = 1:10, 
#' missnumspec = 0)
#' @export bd_loglik
bd_loglik <- function(pars1, pars2, brts, missnumspec, methode = "odeint::runge_kutta_cash_karp54"){
  # pars1 contains model parameters
  # - pars1[1] = la0 = speciation rate
  # - pars1[2] = mu0 = extinction rate
  # - pars1[3] = la1 = parameter in exponential decay of speciation rate, OR K in 
  # - diversity-dependence-like models (default 0)
  # - pars1[4] = mu1 = parameter in exponential decay of extinction rate (
  #   default 0)
  # - pars1[5] = T0 = age at which lambda is lambda0 (default 
  #   T0 = age of phylogeny)
  # pars2 contains settings
  # - pars2[1] = tdmodel = model of time-dependence
  # . tdmodel == 0: no time-dependence
  # . tdmodel == 1: exponential decline in speciation or extinction rate
  # . tdmodel == 2: stepwise decline following diversity-dependence when 
  #   extinction = 0
  # . tdmodel == 3: decline in speciation rate following deterministic logistic 
  #   equation for ddmodel = 1
  # . tdmodel == 4...8: change in speciation/extinction rate following ddmodel = 
  #   1...5, with n = expected number of species
  # - pars2[2] = cond = conditioning on age (0), non-extinction of the clade and 
  #   age (1), number of taxa and age (2), number of taxa (3)
  # - pars2[3] = btorph = likelihood for branching times (0) or phylogeny (1)
  # - pars2[4] = printing of parameters and likelihood (1) or not (0)
  # - pars2[5] = likelihood is for a tree with crown age (2) or stem age (1) - 
  #   stem age not yet implemented for tdmodel == 2
  # - pars2[6] = lx = length of ODE equation (only for tdmodel = 4..8)
  
  rhotaut <- function(tau, t1, pars){
    la0 <- pars[1]
    mu0 <- pars[2]
    la1 <- pars[3]
    mu1 <- pars[4]
    if (la1 == 0 & mu1 == 0) {
      rtt <- mu0 * (tau - t1) - la0 * (tau - t1)
    }
    if (la1 > 0 & mu1 == 0) {
      rtt <- mu0 * (tau - t1) - la0 / la1 * (exp(-la1 * t1) - exp(-la1 * tau))
    }
    if (la1 == 0 & mu1 > 0) {
      rtt <- mu0 / mu1 * (exp(-mu1 * t1) - exp(-mu1 * tau)) - la0 * (tau - t1)
    }
    if (la1 > 0 & mu1 > 0) {
      rtt <- mu0 / mu1 * (exp(-mu1 * t1) - exp(-mu1 * tau)) - la0 / la1 * 
        (exp(-la1 * t1) - exp(-la1 * tau))
    }
    return(rtt)
  }
  
  PtTint <- function(x, t1, pars){
    PtTint <- pars[2] * exp(-pars[4] * x) * exp(rhotaut(x, t1, pars))
    return(PtTint)
  }
  
  ff <- function(t1, pars){
    ff <- 1 / pars[3] + (1 / 2 - 1 / pars[3]) * exp(-(pars[1] - pars[2]) * t1)
    return(ff)
  }
  
  exprhotaut2 <- function(tau, t1, pars){
    rtt <- ff(tau, pars) / ff(t1, pars)
    return(rtt)
  }
  
  intPtTint2 <- function(t1, T, pars){
    intPtTint2 <- pars[2] / ff(t1,pars) * 
      ( (T - t1) / pars[3] + (ff(t1, pars) - ff(T, pars)) / 
          (pars[1] - pars[2]) )
    return(intPtTint2)
  }
  
  tdmodel <- pars2[1]
  if (min(pars1) < 0 | (tdmodel == 4 & pars1[1] <= pars1[2])){
    cat("Negative parameter values are not allowed.\n")
    loglik <- -Inf
    return(as.numeric(loglik))
  } else {
    cond <- pars2[2]
    soc <- pars2[5]
    if (cond == 3){
      soc <- 2
    }
    btorph <- pars2[3]
    verbose <- pars2[4]
    m <- missnumspec
    if (is.na(pars2[6])){
      pars2[6] = Inf
    }
    
    if (soc == 1 & tdmodel == 2){
      pars1new <- c(pars1[1],1e-10,pars1[3])
      pars2new <- c(min(10 * (1 + missnumspec + length(brts)),1000),1,pars2[2:5])
      loglik <- dd_loglik(pars1 = pars1new, pars2 = pars2new, brts = brts, missnumspec = missnumspec, methode = methode)
      return(as.numeric(loglik))
    } else {        
      brts <- sort(abs(as.numeric(brts)), decreasing = TRUE)
      if (brts[length(brts)] == 0){
        brts <- brts[-length(brts)]
      }
      brts2 <- c(-sort(abs(as.numeric(brts)), decreasing = TRUE), 0)
      for(i in 2:length(brts2)) {
        if (brts2[i] == brts2[i - 1]){
          brts2[i] <- brts2[i] + 1E-10
        }
      }
      if (soc == 2){
        brts <- c(brts[1], brts)
      }
      TT <- brts[1]
      t <- TT - brts
      
      S <- length(brts)
      N <- S + m
      np <- length(brts2)
      PtT <- rep(0, S)
      ux <- rep(0, S)
      
      abstol <- 1e-16
      reltol <- 1e-10
      
      T0 <- TT
      la0 <- pars1[1]
      mu0 <- pars1[2]
      if (length(pars1) > 2){
        la1 <- pars1[3]
        K <- pars1[3]
        mu1 <- pars1[4]
        if (length(pars1) == 5){
          T0 <- pars1[5]
        } 
      } else {
        la1 <- 0
        K <- Inf
        mu1 <- 0
        pars1 <- c(pars1, c(0,0))
      }
      
      if(tdmodel == 1){
        la0 <- la0 * exp(-la1 * (T0 - TT))
        mu0 <- mu0 * exp(-mu1 * (T0 - TT))
      }
      
      loglik <- (btorph == 0) * lgamma(S)
      
      if(tdmodel == 0 | (tdmodel == 1 & la1 == 0 & mu1 == 0))
      {
        if(abs(la0 - mu0) < 1E-10){
          lamu <- max(la0, mu0)
          PtT <- 1 / (1 + lamu * (TT - t))
          ux <- lamu * (TT - t) / (1 + lamu * (TT - t))
        } else {
          PtT <- (la0 - mu0) / (la0 - mu0 * exp(-(la0 - mu0) * (TT - t)))
      ux <- la0 * (1 - exp(-(la0 - mu0) * (TT - t))) /
        (la0 - mu0 * exp(-(la0 - mu0) * (TT - t)))
        }
      }
      if (tdmodel == 1 & (la1 != 0 | mu1 != 0)){
        for(i in 1:S){
          PtT[i] <- (1 + stats::integrate(
            PtTint, lower = t[i], upper = TT, 
            t1 = t[i], pars = pars1, subdivisions = 10000L
          )$value)^(-1)
          ux[i] <- 1 - PtT[i] * exp(rhotaut(TT, t[i], pars1))
        }
      } else if (tdmodel == 2){
        if (N > ceiling(K)){
          loglik <- -Inf
          return(as.numeric(loglik))
        }
        la <- pmax(0,la0 * (1 - (2:S) / K))
        dx <- c(abs(diff(brts)),brts[S])[2:S]
        # mpfr(la, precBits = 500)
        # mpfr(dx, precBits = 500)
        ladx <- la * dx
        PtT <- rep(1,S)
        for(i in 2:S){
          ux[i] <- 1 - exp(-sum(ladx[(i - 1):(S - 1)]))
        }
        ux[1] <- ux[2]
      } else if (tdmodel == 3){
        PtT <- 1 / (1 + intPtTint2(t, TT, pars1))
        ux <- 1 - PtT * exprhotaut2(TT, t, pars1)
      } else if (tdmodel == 4){
        Kprime <- la0 / (la0 - mu0) * K
        lx <- min(max(1 + missnumspec, 1 + ceiling(Kprime)), round(pars2[6]))
        if (
          ceiling(Kprime) < N | 
          K >= 0.9 * round(pars2[6]) | 
          K < 1 | 
          la0 < mu0 | 
          mu0 < 0
        ){
          loglik <- -Inf
          return(as.numeric(loglik))
        }
        variables <- rep(0, lx + np - 1)
        variables[soc + 1] <- 1
        expn <- rep(0, np)
        expn[1] <- soc
        latd <- rep(0, np)
        latd[1] <- mu0 + 
          ((la0 - mu0) * expn[1] - expn[1]^2 * (la0 - mu0) / K) / expn[1]
        for(k in 2:np){
          t1 <- brts2[k - 1]
          t2 <- brts2[k]
          variables[lx + k - 1] <- 0
          parsvec = c(pars1[1:min(4, length(pars1))], tdmodel - 3, lx)

          if(startsWith(methode, "odeint::")) {
            variables = dd_integrate_td_odeint(variables, c(t1, t2), parsvec, 1e-10, 1e-16, methode)
          }
          else {
            stop("f95 is gone")
          }
          if (is.na(sum(variables)) | is.nan(sum(variables))) {
            if (verbose) {
              cat('NA or NaN issues encountered.\n')
            }
            loglik <- -Inf
            variables <- rep(-Inf, length(variables))
          } else if (sum(variables) <= 0) {
            if(verbose) cat('Probabilities smaller than 0 encountered\n')
            loglik <- -Inf
            variables <- rep(-Inf,length(variables))
          } 
          #variables2 = y2[2, 2:(lx + k)]
          #if(any(variables2 != variables))
          #{ 
          #   print(variables2 - variables)
          #   unequal <- which(variables != variables2)
          #  if(length(unequal) > 0) stop('Stop here')
          #}
          expn[k] <- sum((0:(lx - 1)) * variables[1:lx])
          expn2 <- sum((0:(lx - 1))^2 * variables[1:lx])
          dEN_dt <- (la0 - mu0) * expn[k] - expn2 * (la0 - mu0) / K
          latd[k] <- mu0 + dEN_dt / expn[k]
        }
        if (sum(variables[1:lx]) < 0.99 ) {
          #print(lx)
          #print(sum(variables[1:lx]))
          cat('Leaking probabilities detected.')
          loglik <- -Inf
        }
        sig <- expn[1:(np - 1)] * variables[(lx + 1):(lx + np - 1)]
        PtT <- 1 / (1 + sig)
        ux <- 1 - PtT * expn[1:(np - 1)] / expn[np]
        if (soc == 2) {
          PtT <- c(PtT[1], PtT)
          ux <- c(ux[1], ux)
        }
      }
      
      if (S > soc | cond == 3) {
        if(tdmodel == 0) {
          lavec <- rep(la0, S - 1)
        } else if (tdmodel == 1) { 
          lavec <- la0 * exp(-la1 * t[2:S])
        } else if (tdmodel == 2) {
          lavec <- la0 * (1 - (1:(S - 1)) / K)
        } else if(tdmodel == 3) {
          lavec <- la0 * (1 - (1 - mu0 / la0) / (K * ff(t[2:S], pars1)))
        } else if(tdmodel == 4) {
          lavec <- latd[1:(np - 1)]
        } else
        {
          stop('The tdmodel you selected does not exist. Please check pars2.')
        }
        if (S > soc){
          loglik <- loglik + sum(log(lavec[soc:length(lavec)])) 
        }
      }
      
      logp <- 0
      if(cond == 1) {
        logp = soc * log(PtT[1])
      }
      if(cond == 2 | cond == 3) {
        if(tdmodel == 2) {
          eps <- 1E-6
          K2 <- K
          if (floor(K) == ceiling(K)) { 
            K2 <- K + eps
          } else if (floor(K + eps) == ceiling(K)) { 
            K2 <- K - eps
          } else if (ceiling(K - eps) == floor(K)) {
            K2 <- K + eps
          }
          s <- (1:N) * pmax(0, la0 * (1 - (1:N) / K2))
          logsdiff <- rep(0, N)
          sgnsdiff <- rep(0, N)
          logsdiff2 <- rep(0, N)
          sgnsdiff2 <- rep(0, N)
          logsdiff3 <- rep(0, N)
          #if(1/la0 > 30)
          #{
          #pB = max(50,round(3/la0))
          #s = mpfr(s, precBits = pB)      
          #sgnsdiff = mpfr(sgnsdiff, precBits = pB)
          #logsdiff = mpfr(logsdiff, precBits = pB)
          #sgnsdiff2 = mpfr(sgnsdiff2, precBits = pB)
          #logsdiff2 = mpfr(logsdiff2, precBits = pB)
          #logsdiff3 = mpfr(logsdiff3, precBits = pB)
          #}
          for (i in 1:N) { 
            sgnsdiff[i] <- prod(sign(s[-c(1,i)] - s[i]))
            logsdiff[i] <- sum(log(abs(s[-c(1,i)] - s[i])))
            sgnsdiff2[i] <- prod(sign(s[-i] - s[i]))
            logsdiff2[i] <- sum(log(abs(s[-i] - s[i])))
            logsdiff3[i] <- -sum(log(abs(1 - s[i]/s[-c(1,i,N)]))) - 
              log(abs(s[N] / s[i] - 1))
          }
          logsdiff3[N] <- -sum(log(abs(1 - s[N] / s[-c(1, N)])))
          if (cond == 2) {
            #         logp <- sum(log(s[2:(N-1)])) + log(sum(sgnsdiff[2:N] * 
            #         exp(-s[2:N] * T - logsdiff[2:N])))
            #         logp <- log(sum(sgnsdiff[2:N] * exp(sum(log(s[2:(N-1)])) - 
            #         s[2:N] * T - logsdiff[2:N])))
            logp <- log(sum(sgnsdiff[2:N] * exp(-s[2:N] * TT + logsdiff3[2:N])))
          } else if (cond == 3) {
            logp <- sum(log(s[1:(N-1)])) + 
              log(sum(sgnsdiff2[1:N] * 1 / s[1:N] * exp(-logsdiff2[1:N])))
          }
        } else {
          # if tdmodel is not 2
          if(cond == 2) {
            logp <- (soc == 2) * log(S + m - 1) + soc * log(PtT[1]) + 
              soc * log(1 - ux[1]) + (S + m - soc) * log(ux[1])
          } else if (cond == 3) {
            logp <- log(PtT[1]) - log(S + missnumspec) - 
              (soc == 2) * log(lavec[1])
          }
        }
      }
      
      loglik <- loglik + sum(log(PtT)) + sum(log(1 - ux)) - logp
      
      if (m > 0) {
        if(tdmodel == 2) {
          cat(
            "Missing species in diversity-dependence models 
              cannot be dealt with by bd_loglik.R\n"
          )
          return(-Inf)
        }
        if(cond == 3) {
          x <- (1 - ux[1]^((0:m)+1)) / (1 - ux[1])
        } else {
          x <- (1:(m + 1)) * ux[1]^(0:m)
        }
        for(j in 2:S) {
          #x = convolve(x,rev((1:(m + 1)) * ux[j]^(0:m)),type = 'open')[1:(m + 1)]
          x <- conv(x, (1:(m + 1)) * ux[j]^(0:m))[1:(m+1)]
        }
        loglik <- loglik + lgamma(S + 1) + lgamma(m + 1) - 
          lgamma(S + m + 1) + log(x[m + 1])
        #loglik = loglik - log(S + m + 1) + log(S + 1) + log(S + m - 1) - log(S - 1)
      }
      if (verbose) {
        if(tdmodel == 0) {
          s1 <- sprintf('Parameters: %f %f', pars1[1], pars1[2])
        } else if (tdmodel == 1) {
          s1 <- sprintf(
            'Parameters: %f %f %f %f', pars1[1], pars1[2], pars1[3], pars1[4]
          )
        } else if (tdmodel >= 2) {
          s1 <- sprintf('Parameters: %f %f %f', pars1[1], pars1[2], pars1[3])
        }
        s2 <- sprintf(', Loglikelihood: %f', loglik)
        cat(s1, s2, "\n", sep = "")
        utils::flush.console()
      }
    }
    loglik <- as.numeric(loglik)
    if (is.nan(loglik) | is.na(loglik)) {
      loglik <- -Inf
    }
    return(loglik)
  }
}
