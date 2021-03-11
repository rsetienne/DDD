edd_update_lamu <- function(lamu, Phi, pars, model) {
  if (model == "dsce1") {
    if (length(pars) != 1) {
      stop("incorrect parameter(s)")
    }
    #dependent speciation, constant extinction
    K <- pars
    newla <- max(0, lamu[1, 1] * (1 - Phi / K))
    newmu <- lamu[1, 2]
  } else if (model == "dsce2") {
    if (length(pars) != 3) {
      stop("incorrect parameter(s)")
    }
    #dependent speciation, constant extinction
    N <- pars[1]
    beta_N <- pars[2]
    beta_phi <- pars[3]
    newla <- max(0, lamu[1, 1] + beta_N * N + beta_phi * Phi)
    newmu <- lamu[1, 2]
  } else if (model == "dsde1") {
    if (length(pars) != 2) {
      stop("incorrect parameter(s)")
    }
    mu0 <- pars[1]
    K <- pars[2]
    newla <- max(0, lamu[1, 1] * (1 - Phi / K))
    newmu <- max(0, mu0 * (Phi / K))
  } else if (model == "dsde2") {
    if (length(pars) != 5) {
      stop("incorrect parameter(s)")
    }
    #dependent speciation, constant extinction
    N <- pars[1]
    beta_N <- pars[2]
    beta_phi <- pars[3]
    gamma_N <- pars[4]
    gamma_phi <- pars[5]
    newla <- max(0, lamu[1, 1] + beta_N * N + beta_phi * Phi)
    newmu <- max(0, lamu[1, 2] + gamma_N * N + gamma_phi * Phi)
  }
  
  return(c(newla, newmu))
}

edd_sum_rates <- function(lamu, ED, i) {
  # calculate the rates of each lineages and sum them up
  return(sum(ED) * sum(lamu[i, ]))
}

edd_sample_event <- function(lamu, ED, i, EDPars, EDModel, EDMethod) {
  if (EDModel == "dsce") {
    if (EDMethod == "scale") {
      if (is.null(EDPars) == FALSE) {
        stop("incorrect parameter(s)")
      }
      
      ED <- range01(ED)
      
      if ((lamu[i - 1, 1] - lamu[i, 1]) >= 0) {
        rspec <- lamu[i, 1] * ED
        rfake_spec <- (lamu[i - 1, 1] - lamu[i, 1]) * ED
      } else {
        rspec <- lamu[i , 1] * ED
        rfake_spec <- 0
      }
      
      if ((lamu[i - 1, 2] - lamu[i, 2]) >= 0) {
        rext <- lamu[i, 2] + ED - ED
        rfake_ext <- (lamu[i - 1, 2] - lamu[i, 2]) + ED - ED
      } else {
        rext <- rlamu[i, 2] + ED - ED
        rfake_ext <- 0
      }
    }
    
    if (EDMethod == "add") {
      if (length(EDPars) != 1) {
        stop("incorrect parameter(s)")
      }
      
      rho_la <- EDPars[1]
      
      if ((lamu[i - 1, 1] - lamu[i, 1]) >= 0) {
        rspec <- lamu[i, 1] + (rho_la * ED)
        rspec[rspec < 0] <- 0
        rfake_spec <-
          (lamu[i - 1, 1] - lamu[i, 1]) + (rho_la * ED)
        rfake_spec[rfake_spec < 0] <- 0
      } else {
        rspec <- lamu[i, 1] + (rho_la * ED)
        rspec[rspec < 0] <- 0
        rfake_spec <- 0
      }
      
      if ((lamu[i - 1, 2] - lamu[i, 2]) >= 0) {
        rext <- lamu[i, 2] + ED - ED
        rfake_ext <- (lamu[i - 1, 2] - lamu[i, 2]) + ED - ED
      } else {
        rext <- rlamu[i, 2] + ED - ED
        rfake_ext <- 0
      }
    }
  }
  
  if (EDModel == "dsde") {
    if (EDMethod == "scale") {
      if (is.null(EDPars) != FALSE) {
        stop("incorrect parameter(s)")
      }
      
      ED <- range01(ED)
      
      if ((lamu[i - 1, 1] - lamu[i, 1]) >= 0) {
        rspec <- lamu[i, 1] * ED
        rfake_spec <- (lamu[i - 1, 1] - lamu[i, 1]) * ED
      } else {
        rspec <- lamu[i , 1] * ED
        rfake_spec <- 0
      }
      
      if ((lamu[i - 1, 2] - lamu[i, 2]) >= 0) {
        rext <- lamu[i, 2] / ED
        rfake_ext <- (lamu[i - 1, 2] - lamu[i, 2]) / ED
      } else {
        rext <- lamu[i, 2] / ED
        rfake_ext <- 0
      }
    }
    
    if (EDMethod == "add") {
      if (length(EDPars) != 2) {
        stop("incorrect parameter(s)")
      }
      
      rho_la <- EDPars[1]
      rho_mu <- EDPars[2]
      
      if ((lamu[i - 1, 1] - lamu[i, 1]) >= 0) {
        rspec <- lamu[i, 1] + (rho_la * ED)
        rspec[rspec < 0] <- 0
        rfake_spec <-
          (lamu[i - 1, 1] - lamu[i, 1]) + (rho_la * ED)
        rfake_spec[rfake_spec < 0] <- 0
      } else {
        rspec <- lamu[i, 1] + (rho_la * ED)
        rspec[rspec < 0] <- 0
        rfake_spec <- 0
      }
      
      if ((lamu[i - 1, 2] - lamu[i, 2]) >= 0) {
        rext <-  lamu[i, 2] + (rho_mu * ED)
        rext[rext < 0] <- 0
        rfake_ext <-
          (lamu[i - 1, 2] - lamu[i, 2]) + (rho_mu * ED)
        rfake_ext[rfake_ext < 0] <- 0
      } else {
        rext <- lamu[i, 2] + (rho_mu * ED)
        rext[rext < 0] <- 0
        rfake_ext <- 0
      }
    }
  }
  
  # add prefix to names of the rates to distinguish between different event types
  names(rspec) <- paste("rs", names(rspec), sep = "_")
  names(rext) <- paste("re", names(rext), sep = "_")
  names(rfake_spec) <- paste("rfs", names(rfake_spec), sep = "_")
  names(rfake_ext) <- paste("rfe", names(rfake_ext), sep = "_")
  
  rates <- c(rspec, rext, rfake_spec, rfake_ext)
  events <- names(rates)
  
  return(DDD::sample2(events, 1, prob = rates))
}

edd_sim <- function (pars,
                     age,
                     model = "dsce1",
                     metric = "pd",
                     offset = "none",
                     EDPars = NULL,
                     EDModel = "dsce",
                     EDMethod = "add") {
  
  if (pars[1] <= 0 | pars[2] <= 0) {
    stop('per species rates should be positive')
  }
  
  if (length(pars) == 3) {
    if (pars[3] <= 0) {
      stop('clade level carrying capacity should be positive')
    }
  }
  
  if (length(pars) == 4) {
    if (pars[2] <= 0) {
      stop('coefficient for extinction should be positive')
    }
  }
  
  if (length(EDPars) == 1 && EDModel == "dsde") {
    stop('incorrect parameter(s)')
  }
  
  if (length(EDPars) == 2 && EDModel == "dsce") {
    stop('incorrect parameter(s)')
  }
  
  if (length(pars) == 3 && model == "dsce1") {
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
      Phi <- rep(0, 1) # PD
      K <- pars[3]
      linlist <- c(-1, 2)
      newL <- 2
      lamu <- matrix(c(pars[1], pars[2]), ncol = 2)
      Phi[i] <- 0
      # add new intial parameter ED
      ED <- c(1, 1)
      names(ED) <- c("t1", "t2")
      
      # ED values determine the total rates
      t[i + 1] <-
        t[i] + stats::rexp(1, edd_sum_rates(lamu, ED, i))
      
      # main simulation circle
      while (t[i + 1] <= age) {
        i <- i + 1
        
        # new algorithm to deal with non-constant rates
        t_new <- t[i - 1] + stats::rexp(1, edd_sum_rates(lamu, ED, 1))
        
        if (t_new <= t[i])
          t[i] <- t_new
        
        # new offset methods
        if (offset == "none") {
          Phi[i] <- L2Phi(L, t[i], metric)
        } else if (offset == "simtime") {
          Phi[i] <- L2Phi(L, t[i], metric) - t[i]
        } else if (offset == "nspecies") {
          Phi[i] <- L2Phi(L, t[i], metric) / N[i - 1]
        } else{
          stop("no such offset method")
        }
        
        # also update ED at the start of each iteration
        ED <- L2ED(L, t[i])
        
        # currently the same as pdd_sim
        lamu <-
          rbind(lamu, edd_update_lamu(lamu, Phi[i], K, model))
        # sample an event containing info of focal lineage and event type
        event <- edd_sample_event(lamu, ED, i, EDPars, EDModel, EDMethod)
        # extract lineage number from event string
        numL <- sub(".*_t", "", event)
        # find out the lineage to match the sign in linlist
        ranL <- linlist[match(numL,abs(linlist))]
        
        # detect event type
        if (grepl("rs", event, fixed = TRUE)) {
          N[i] <- N[i - 1] + 1
          newL <- newL + 1
          L <- rbind(L, c(t[i], ranL, sign(ranL) * newL,
                          -1))
          linlist <- c(linlist, sign(ranL) * newL)
        } else if (grepl("re", event, fixed = TRUE)) {
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
              edd_update_lamu(lamu, Phi[i], K, model)
          }
        } else if (grepl("rfs", event, fixed = TRUE) |
                   grepl("rfe", event, fixed = TRUE)) {
          N[i] <- N[i - 1]
        }
        if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
          t[i + 1] <- Inf
        } else {
          t[i + 1] <- t[i] + stats::rexp(1, edd_sum_rates(lamu, ED, i))
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
    return(out)
  }
  
  if (length(pars) == 4 && model == "dsce2") {
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
      Phi <- rep(0, 1) # Phylogenetic metrices
      Nbetas <- c(N, pars[3], pars[4])
      linlist <- c(-1, 2)
      newL <- 2
      lamu <- matrix(c(pars[1], pars[2]), ncol = 2)
      Phi[i] <- 0
      # add new intial parameter ED
      ED <- c(1, 1)
      names(ED) <- c("t1", "t2")
      
      # ED values determine the total rates
      t[i + 1] <-
        t[i] + stats::rexp(1, edd_sum_rates(lamu, ED, i))
      
      # main simulation circle
      while (t[i + 1] <= age) {
        i <- i + 1
        
        # new algorithm to deal with non-constant rates
        t_new <- t[i - 1] + stats::rexp(1, edd_sum_rates(lamu, ED, 1))
        
        if (t_new <= t[i])
          t[i] <- t_new
        
        if (offset == "none") {
          Phi[i] <- L2Phi(L, t[i], metric)
        } else if (offset == "simtime") {
          Phi[i] <- L2Phi(L, t[i], metric) - t[i]
        } else if (offset == "nspecies") {
          Phi[i] <- L2Phi(L, t[i], metric) / N[i - 1]
        } else{
          stop("no such offset method")
        }
        
        # also update ED at the start of each iteration
        ED <- L2ED(L, t[i])
        
        # currently the same as pdd_sim
        lamu <-
          rbind(lamu, edd_update_lamu(lamu, Phi[i], Nbetas, model))
        # sample an event containing info of focal lineage and event type
        event <- edd_sample_event(lamu, ED, i, EDPars, EDModel, EDMethod)
        # extract lineage number from event string
        numL <- sub(".*_t", "", event)
        # find out the lineage to match the sign in linlist
        ranL <- linlist[match(numL,abs(linlist))]
        
        if (grepl("rs", event, fixed = TRUE)) {
          N[i] <- N[i - 1] + 1
          newL <- newL + 1
          L <- rbind(L, c(t[i], ranL, sign(ranL) * newL,
                          -1))
          linlist <- c(linlist, sign(ranL) * newL)
        } else if (grepl("re", event, fixed = TRUE)) {
          N[i] <- N[i - 1] - 1
          L[abs(ranL), 4] <- t[i]
          w <- which(linlist == ranL)
          linlist <- linlist[-w]
          linlist <- sort(linlist)
          
          if (sum(linlist < 0) == 0 |
              sum(linlist > 0) == 0) {
            
          } else {
            Phi[i] <- L2Phi(L, t[i], metric)
            lamu[i,] <-
              edd_update_lamu(lamu, Phi[i], Nbetas, model)
          }
        } else if (grepl("rfs", event, fixed = TRUE) |
                   grepl("rfe", event, fixed = TRUE)) {
          N[i] <- N[i - 1]
        }
        if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
          t[i + 1] <- Inf
        } else {
          t[i + 1] <- t[i] + stats::rexp(1, edd_sum_rates(lamu, ED, i))
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
    return(out)
  }
  
  if (length(pars) == 3 && model == "dsde1") {
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
      Phi <- rep(0, 1) # PD
      K <- pars[3]
      mu0 <- pars[2]
      linlist <- c(-1, 2)
      newL <- 2
      lamu <- matrix(c(pars[1], 0), ncol = 2)
      Phi[i] <- 0
      # add new intial parameter ED
      ED <- c(1, 1)
      names(ED) <- c("t1", "t2")
      
      t[i + 1] <-
        t[i] + stats::rexp(1, edd_sum_rates(lamu, ED, i))
      
      # main simulation circle
      while (t[i + 1] <= age) {
        i <- i + 1
        
        # new algorithm to deal with non-constant rates
        t_new <- t[i - 1] + stats::rexp(1, edd_sum_rates(lamu, ED, 1))
        
        if (t_new <= t[i])
          t[i] <- t_new
        
        if (offset == "none") {
          Phi[i] <- L2Phi(L, t[i], metric)
        } else if (offset == "simtime") {
          Phi[i] <- L2Phi(L, t[i], metric) - t[i]
        } else if (offset == "nspecies") {
          Phi[i] <- L2Phi(L, t[i], metric) / N[i - 1]
        } else{
          stop("no such offset method")
        }
        
        # also update ED at the start of each iteration
        ED <- L2ED(L, t[i])
        
        # currently the same as pdd_sim
        lamu <-
          rbind(lamu, edd_update_lamu(lamu, Phi[i], c(mu0, K), model))
        # sample an event containing info of focal lineage and event type
        event <- edd_sample_event(lamu, ED, i, EDPars, EDModel, EDMethod)
        # extract lineage number from event string
        numL <- sub(".*_t", "", event)
        # find out the lineage to match the sign in linlist
        ranL <- linlist[match(numL,abs(linlist))]
        
        if (grepl("rs", event, fixed = TRUE)) {
          N[i] <- N[i - 1] + 1
          newL <- newL + 1
          L <- rbind(L, c(t[i], ranL, sign(ranL) * newL,
                          -1))
          linlist <- c(linlist, sign(ranL) * newL)
        } else if (grepl("re", event, fixed = TRUE)) {
          N[i] <- N[i - 1] - 1
          L[abs(ranL), 4] <- t[i]
          w <- which(linlist == ranL)
          linlist <- linlist[-w]
          linlist <- sort(linlist)
          
          if (sum(linlist < 0) == 0 |
              sum(linlist > 0) == 0) {
            
          } else {
            Phi[i] <- L2Phi(L, t[i], metric)
            lamu[i,] <-
              edd_update_lamu(lamu, Phi[i], c(mu0, K), model)
          }
        } else if (grepl("rfs", event, fixed = TRUE) |
                   grepl("rfe", event, fixed = TRUE)) {
          N[i] <- N[i - 1]
        }
        if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
          t[i + 1] <- Inf
        } else {
          t[i + 1] <- t[i] + stats::rexp(1, edd_sum_rates(lamu, ED, i))
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
    return(out)
  }
  
  if (length(pars) == 6 && model == "dsde2") {
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
      Phi <- rep(0, 1) # Phylogenetic metrices
      Nbg <- c(N, pars[3], pars[4], pars[5], pars[6])
      linlist <- c(-1, 2)
      newL <- 2
      lamu <- matrix(c(pars[1], pars[2]), ncol = 2)
      Phi[i] <- 0
      ED <- c(1, 1)
      names(ED) <- c("t1", "t2")
      
      t[i + 1] <-
        t[i] + stats::rexp(1, edd_sum_rates(lamu, ED, i))
      
      # main simulation circle
      while (t[i + 1] <= age) {
        i <- i + 1
        
        # new algorithm to deal with non-constant rates
        t_new <- t[i - 1] + stats::rexp(1, edd_sum_rates(lamu, ED, 1))
        
        if (t_new <= t[i])
          t[i] <- t_new
        
        if (offset == "none") {
          Phi[i] <- L2Phi(L, t[i], metric)
        } else if (offset == "simtime") {
          Phi[i] <- L2Phi(L, t[i], metric) - t[i]
        } else if (offset == "nspecies") {
          Phi[i] <- L2Phi(L, t[i], metric) / N[i - 1]
        } else{
          stop("no such offset method")
        }
        
        # also update ED at the start of each iteration
        ED <- L2ED(L, t[i])
        
        # currently the same as pdd_sim
        lamu <-
          rbind(lamu, edd_update_lamu(lamu, Phi[i], Nbg, model))
        # sample an event containing info of focal lineage and event type
        event <- edd_sample_event(lamu, ED, i, EDPars, EDModel, EDMethod)
        # extract lineage number from event string
        numL <- sub(".*_t", "", event)
        # find out the lineage to match the sign in linlist
        ranL <- linlist[match(numL,abs(linlist))]
        
        if (grepl("rs", event, fixed = TRUE)) {
          N[i] <- N[i - 1] + 1
          newL <- newL + 1
          L <- rbind(L, c(t[i], ranL, sign(ranL) * newL,
                          -1))
          linlist <- c(linlist, sign(ranL) * newL)
        } else if (grepl("re", event, fixed = TRUE)) {
          N[i] <- N[i - 1] - 1
          L[abs(ranL), 4] <- t[i]
          w <- which(linlist == ranL)
          linlist <- linlist[-w]
          linlist <- sort(linlist)
          
          if (sum(linlist < 0) == 0 |
              sum(linlist > 0) == 0) {
            
          } else {
            Phi[i] <- L2Phi(L, t[i], metric)
            lamu[i,] <-
              edd_update_lamu(lamu, Phi[i], Nbg, model)
          }
        } else if (grepl("rfs", event, fixed = TRUE) |
                   grepl("rfe", event, fixed = TRUE)) {
          N[i] <- N[i - 1]
        }
        if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
          t[i + 1] <- Inf
        } else {
          t[i + 1] <- t[i] + stats::rexp(1, edd_sum_rates(lamu, ED, i))
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
    return(out)
  }
} 