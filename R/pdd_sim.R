pdd_update_lamu <- function(lamu, Phi, pars, model) {
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

pdd_sum_rates <- function(lamu, N, i) {
    return(lamu[i, 1] * N + lamu[i, 2] * N)
}

pdd_sample_event <- function(lamu, N, i) {
    events = c("spec", "ext", "fake_spec", "fake_ext")
    
    if ((lamu[i - 1, 1] - lamu[i, 1]) >= 0) {
        rspec <- lamu[i, 1]
        rfake_spec <- lamu[i - 1, 1] - lamu[i, 1]
    } else {
        rspec <- lamu[i , 1]
        rfake_spec <- 0
    }
    if ((lamu[i - 1, 2] - lamu[i, 2]) >= 0) {
        rext <- lamu[i, 2]
        rfake_ext <- lamu[i - 1, 2] - lamu[i, 2]
    } else {
        rext <- lamu[i, 2]
        rfake_ext <- 0
    }
    
    rates <- c(rspec, rext, rfake_spec, rfake_ext)
    
    return(DDD::sample2(events, 1, prob = rates))
}

pdd_sim <- function (pars,
                     age,
                     model = "dsce1",
                     metric = "pd",
                     offset = "none") {
    if (pars[1] < 0 | pars[2] < 0) {
        stop('per species rates should be positive')
    }
    
    if (length(pars) == 3) {
        if (pars[3] < 0) {
            stop('clade level carrying capacity should be positive')
        }
    }
    
    if (length(pars) == 4) {
        if (pars[2] < 0) {
            stop('coefficient for extinction should be positive')
        }
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
            
            t[i + 1] <-
                t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))
            
            # main simulation circle
            while (t[i + 1] <= age) {
                i <- i + 1
                ranL <- sample2(linlist, 1)
                
                # new algorithm to deal with non-constant rates
                t_new <- t[i - 1] + stats::rexp(1, pdd_sum_rates(lamu, N, 1))
                
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
                
                lamu <-
                    rbind(lamu, pdd_update_lamu(lamu, Phi[i], K, model))
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
                            pdd_update_lamu(lamu, Phi[i], K, model)
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
            
            t[i + 1] <-
                t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))
            
            # main simulation circle
            while (t[i + 1] <= age) {
                i <- i + 1
                ranL <- sample2(linlist, 1)
                
                # new algorithm to deal with non-constant rates
                t_new <- t[i - 1] + stats::rexp(1, pdd_sum_rates(lamu, N, 1))
                
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
                
                lamu <-
                    rbind(lamu,
                          pdd_update_lamu(lamu, Phi[i], Nbetas, model))
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
                        
                    } else {
                        Phi[i] <- L2Phi(L, t[i], metric)
                        lamu[i,] <-
                            pdd_update_lamu(lamu, Phi[i], Nbetas, model)
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
            
            t[i + 1] <-
                t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))
            
            # main simulation circle
            while (t[i + 1] <= age) {
                i <- i + 1
                ranL <- sample2(linlist, 1)
                
                # new algorithm to deal with non-constant rates
                t_new <- t[i - 1] + stats::rexp(1, pdd_sum_rates(lamu, N, 1))
                
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
                
                lamu <-
                    rbind(lamu,
                          pdd_update_lamu(lamu, Phi[i], c(mu0, K), model))
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
                        
                    } else {
                        Phi[i] <- L2Phi(L, t[i], metric)
                        lamu[i,] <-
                            pdd_update_lamu(lamu, Phi[i], c(mu0, K), model)
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
            
            t[i + 1] <-
                t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))
            
            # main simulation circle
            while (t[i + 1] <= age) {
                i <- i + 1
                ranL <- sample2(linlist, 1)
                
                # new algorithm to deal with non-constant rates
                t_new <- t[i - 1] + stats::rexp(1, pdd_sum_rates(lamu, N, 1))
                
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
                
                lamu <-
                    rbind(lamu,
                          pdd_update_lamu(lamu, Phi[i], Nbg, model))
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
                        
                    } else {
                        Phi[i] <- L2Phi(L, t[i], metric)
                        lamu[i,] <-
                            pdd_update_lamu(lamu, Phi[i], Nbg, model)
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