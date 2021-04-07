edd_update_lamu <- function(ED, params, model, i) {
  if (model == "dsce2") {
    if (length(params) != 5) {
      stop("incorrect parameter(s)")
    }
    #dependent speciation, constant extinction
    N <- params[1]
    la0 <- params[2]
    mu0 <- params[3]
    beta_N <- params[4]
    beta_phi <- params[5]
    newlas <- la0 + beta_N * N + beta_phi * dplyr::select_if(ED[i, ], !is.na(ED[i, ]))
    newlas[newlas < 0] <- 0
    newmus <- NA
  } else if (model == "dsde2") {
    if (length(params) != 7) {
      stop("incorrect parameter(s)")
    }
    #dependent speciation, constant extinction
    N <- params[1]
    la0 <- params[2]
    mu0 <- params[3]
    beta_N <- params[4]
    beta_phi <- params[5]
    gamma_N <- params[6]
    gamma_phi <- params[7]
    newlas <- la0 + beta_N * N + beta_phi * dplyr::select_if(ED[i, ], !is.na(ED[i, ]))
    newlas[newlas < 0] <- 0
    newmus <- mu0 + gamma_N * N + gamma_phi * dplyr::select_if(ED[i, ], !is.na(ED[i, ]))
    newmus[newmus < 0] <- 0
  }
  
  return(list(newlas = newlas, newmus = newmus))
}

edd_sum_rates <- function(las, mus, i) {
  return(sum(dplyr::select_if(las[i, ], !is.na(las[i, ]))) + sum(dplyr::select_if(mus[i, ], !is.na(mus[i, ]))))
}

edd_sample_event <- function(las, mus, i) {
  # add prefix to names of the rates to distinguish between different event types
  rspec <- dplyr::select_if(las[i, ], !is.na(las[i, ]))
  rext <- dplyr::select_if(mus[i, ], !is.na(mus[i, ]))
  names(rspec) <- paste("rs", names(rspec), sep = "_")
  names(rext) <- paste("re", names(rext), sep = "_")
  
  rates <- c(rspec, rext)
  events <- names(rates)
  
  return(DDD::sample2(events, 1, prob = rates))
}

edd_sim <- function (pars,
                     age,
                     model = "dsce1",
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
      
      if (metric == "ed"){
        ED <- c(1, 1)
        names(ED) <- c("t1", "t2")
        las <- ED * pars[1]
        las <- tibble::as_tibble_row(las)
        mus <- ED * pars[2]
        mus <- tibble::as_tibble_row(mus)
      }else{
        lamu <- matrix(c(pars[1], pars[2]), ncol = 2)
        Phi[i] <- 0
      }

      # ED values determine the total rates
      if (metric == "ed") {
        t[i + 1] <-
          t[i] + stats::rexp(1, edd_sum_rates(las, mus, i))
      } else{
        t[i + 1] <-
          t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))
      }
      
      # main simulation circle
      while (t[i + 1] <= age) {
        i <- i + 1
        
        # new algorithm to deal with non-constant rates
        if (metric == "ed") {
          t_new <-
            t[i - 1] + stats::rexp(1, edd_sum_rates(las, mus, i - 1))
        } else{
          t_new <- t[i - 1] + stats::rexp(1, pdd_sum_rates(lamu, N, i - 1))
        }
        
        if (t_new <= t[i])
          t[i] <- t_new
        
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
          ED <- dplyr::bind_rows(ED, L2ED(L, t[i]))
          newlamus <- edd_update_lamu(ED, params, model, i)
          las <- dplyr::bind_rows(las, newlamus$newlas)
          
          if(model == "dsce2"){
            mus <- dplyr::bind_rows(mus, mus)
          }else{
            mus <- dplyr::bind_rows(mus, newlamus$newmus)
          }
          
          # sample an event containing info of focal lineage and event type
          event <-
            edd_sample_event(las, mus, i)
          # extract lineage number from event string
          numL <- sub(".*_t", "", event)
          # find out the lineage to match the sign in linlist
          ranL <- linlist[match(numL, abs(linlist))]
          
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
          } else if (grepl("rfs", event, fixed = TRUE) |
                     grepl("rfe", event, fixed = TRUE)) {
            N[i] <- N[i - 1]
          }
          
          if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
            t[i + 1] <- Inf
          } else {
            t[i + 1] <- t[i] + stats::rexp(1, edd_sum_rates(las, mus, i))
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
    
    if (metric == "ed") {
      ED <- dplyr::bind_rows(ED, L2ED(L, t[i + 1]))
      newlamus <- edd_update_lamu(ED, params, model, i + 1)
      las <- dplyr::bind_rows(las, newlamus$newlas)
      las <- las[-1, ]
      if(model == "dsce2"){
        mus <- dplyr::bind_rows(mus, mus)
        mus <- mus[-1, ]
      }else{
        mus <- dplyr::bind_rows(mus, newlamus$newmus)
        mus <- mus[-1, ]
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
          ED = ED,
          las = las,
          mus = mus,
          LTT = LTT
        )
    }else{
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