edd_update_lamu <- function(ED, ED_max, params, model) {
  N <- params[1]
  la0 <- params[2]
  mu0 <- params[3]
  beta_N <- params[4]
  beta_phi <- params[5]
  if (model == "dsce2") {
    if (length(params) != 5) {
      stop("incorrect parameter(s)")
    }
    #dependent speciation, constant extinction
    if (beta_phi < 0) {
      newlas <- pmax(0,la0 + beta_N * N + beta_phi * ED)
    } else {
      newlas <- pmax(0,la0 + beta_N * N + beta_phi * ED_max)
    }
    newmus <- rep(mu0,length(newlas))
  } else if (model == "dsde2") {
    if (length(params) != 7) {
      stop("incorrect parameter(s)")
    }
    #dependent speciation, dependent extinction
    gamma_N <- params[6]
    gamma_phi <- params[7]
    if (beta_phi < 0) {
      newlas <- pmax(0,la0 + beta_N * N + beta_phi * ED)
    } else {
      newlas <- pmax(0,la0 + beta_N * N + beta_phi * ED_max)
    }
    if (gamma_phi < 0) {
      newmus <- pmax(0,mu0 + gamma_N * N + gamma_phi * ED)
    } else {
      newmus <- pmax(0,mu0 + gamma_N * N + gamma_phi * ED_max)
    }
  }
  return(list(newlas = newlas, newmus = newmus))
}

edd_sum_rates <- function(las, mus) {
  return(sum(las) + sum(mus))
}

edd_sample_event <- function(las, mus, linlist) {
  events <- 1:(2 * length(linlist))
  return(DDD::sample2(events, 1, prob = c(las,mus)))
}

edd_sim <- function (pars,
                     age,
                     model = "dsce2",
                     metric = "pd",
                     offset = "none") {
  if (pars[1] <= 0 | pars[2] <= 0) {
    stop('per species rates should be positive')
  }
  if (model == "dsce2" && length(pars) != 4) {
    stop('incorrect parameters')
  }
  if (model == "dsde2" && length(pars) != 6) {
    stop('incorrect parameters')
  }
  if (pars[2] <= 0) {
    stop('coefficient for extinction should be positive')
  }
  if (metric != "pd" && offset != "none"){
    stop('only pd metric has offset methods')
  }

    done <- 0
    while (done == 0) {
      # initialization
      t <- rep(0, 1)
      L <- matrix(0, 2, 4)
      i <- 1
      t[1] <- 0
      N <- 2
      L[1, 1:4] <- c(0, 0, -1, -1)
      L[2, 1:4] <- c(0, -1, 2, -1)
      Phi <- rep(0, 1) 
      
      if (metric == "ed") {
        params <- c(N, pars)
      } else{
        params <- c(N, pars[-c(1, 2)])
      }

      linlist <- c(-1, 2)
      newL <- 2
      
      # controlling significant digits in tibble objects
      options(pillar.sigfig = 10)
      
      if (metric == "ed") {
        ED <- c(0, 0)
        ED_max <- as.vector(L2ED(L,age))
        lamu <- edd_update_lamu(ED, ED_max, params, model)
      } else {
        lamu <- matrix(c(pars[1], pars[2]), ncol = 2)
        Phi[i] <- 0
      }
      # ED values determine the total rates
      if (metric == "ed") {
        t[i + 1] <-
          t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))
      } else {
        t[i + 1] <-
          t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))
      }
      
      # main simulation circle
      while (t[i + 1] <= age) {
        i <- i + 1
        #print(i)
        #print(t[i])
        if (metric == "pd"){
          if (offset == "none") {
            Phi[i] <- L2Phi(L, t[i], metric)
          } else if (offset == "simtime") {
            Phi[i] <- L2Phi(L, t[i], metric) - t[i]
          } else if (offset == "nspecies") {
            Phi[i] <- L2Phi(L, t[i], metric) / N[i - 1]
          } else{
            stop("no such offset method")
          }
        }
        
        if (metric == "ed") {
          ED <- as.vector(L2ED(L, t[i]))
          lamu_real <- edd_update_lamu(ED, ED, params, model)
          # sample whether event is fake or real
          #print(lamu$newlas)
          event <- sample(c('real','fake'),
                          1,
                          prob = c(sum(lamu_real$newlas + lamu_real$newmus),sum(lamu$newlas - lamu_real$newlas + lamu$newmus - lamu_real$newmus))
                          )
          #print(event)
          if (event == 'real') {
            # sample an event containing info of focal lineage and event type
            event <- edd_sample_event(lamu_real$newlas, lamu_real$newmus, linlist)
            #print(event)
            ranL <- c(linlist,linlist)[event]
            
            if (event <= length(linlist)) {
              N[i] <- N[i - 1] + 1
              newL <- newL + 1
              L <- rbind(L, c(t[i], ranL, sign(ranL) * newL, -1))
              linlist <- c(linlist, sign(ranL) * newL)
            } else {
              N[i] <- N[i - 1] - 1
              L[abs(ranL), 4] <- t[i]
              w <- which(linlist == ranL)
              linlist <- linlist[-w]
              linlist <- sort(linlist)
            }
          } else {
            N[i] <- N[i - 1]
          }  
          #print(linlist)
          if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
            t[i + 1] <- Inf
          } else {
            ED <- as.vector(L2ED(L, t[i]))
            ED_max <- as.vector(L2ED(L, age))
            params[1] <- N[i]
            lamu <- edd_update_lamu(ED, ED_max, params, model)
            #print(lamu)
            t[i + 1] <- t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))
          }
        } else {
          lamu <-
            rbind(lamu, pdd_update_lamu(lamu, Phi[i], params, model))
          event <- pdd_sample_event(lamu, N, i)
          
          if (event == "spec") {
            N[i] <- N[i - 1] + 1
            newL <- newL + 1
            L <- rbind(L, c(t[i], ranL, sign(ranL) * newL,
                            -1))
            linlist <- c(linlist, sign(ranL) * newL)
          } else if (event == "ext") {
            N[i] <- N[i - 1] - 1
            L[abs(ranL), 4] <- t[i]
            w <- which(linlist == ranL)
            linlist <- linlist[-w]
            linlist <- sort(linlist)
            
            if (sum(linlist < 0) == 0 |
                sum(linlist > 0) == 0) {
              # when one whole crown branch is extinct, do nothing
            } else {
              Phi[i] <- L2Phi(L, t[i], metric)
              lamu[i,] <-
                pdd_update_lamu(lamu, Phi[i], params, model)
            }
          } else if (event == "fake_spec" |
                     event == "fake_ext") {
            N[i] <- N[i - 1]
          }
          if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
            t[i + 1] <- Inf
          } else {
            t[i + 1] <- t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))
          }
        }
      }
      
      if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
        done <- 0
      } else {
        done <- 1
      }
    }
    
    L[, 1] <- age - c(L[, 1])
    notmin1 <- which(L[, 4] != -1)
    L[notmin1, 4] <- age - c(L[notmin1, 4])
    L[which(L[, 4] == age + 1), 4] <- -1
    tes <- L2phylo(L, dropextinct = T)
    tas <- L2phylo(L, dropextinct = F)
    brts <- L2brts(L, dropextinct = T)
    
    if (metric == "ed") {
      LTT <-
        data.frame("time" = t[-i],
                   "N" = N)
      out <-
        list(
          tes = tes,
          tas = tas,
          L = L,
          brts = brts,
          LTT = LTT
        )
    } else {
      lamuphis <-
        data.frame(
          "time" = t[-i],
          "lambda" = lamu[, 1],
          "mu" = lamu[, 2],
          "Phi" = Phi,
          "N" = N
        )
      out <-
        list(
          tes = tes,
          tas = tas,
          L = L,
          brts = brts,
          lamuphis = lamuphis
        )
    }

    return(out)
} 