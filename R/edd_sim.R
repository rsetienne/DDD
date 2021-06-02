edd_pars_check <- function(pars, age, model, metric, offset) {
  #pars range check
  if (pars[1] <= 0 | pars[2] <= 0) {
    stop('per species rates should be positive')
  }
  if (pars[2] <= 0) {
    stop('coefficient for extinction should be positive')
  }
  #pars and model match check
  if (model == "dsce2" && length(pars) != 4) {
    stop('this model requires four parameters')
  }
  if (model == "dsde2" && length(pars) != 6) {
    stop('this model requires six parameters')
  }
  #metric and offset match check
  if (metric != "pd" && offset != "none") {
    stop('only pd metric has offset methods')
  }
}

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
      newlas <- pmax(0, la0 + beta_N * N + beta_phi * ED)
    } else {
      newlas <- pmax(0, la0 + beta_N * N + beta_phi * ED_max)
    }
    newmus <- rep(mu0, length(newlas))
  } else if (model == "dsde2") {
    if (length(params) != 7) {
      stop("incorrect parameter(s)")
    }
    #dependent speciation, dependent extinction
    gamma_N <- params[6]
    gamma_phi <- params[7]
    if (beta_phi < 0) {
      newlas <- pmax(0, la0 + beta_N * N + beta_phi * ED)
    } else {
      newlas <- pmax(0, la0 + beta_N * N + beta_phi * ED_max)
    }
    if (gamma_phi < 0) {
      newmus <- pmax(0, mu0 + gamma_N * N + gamma_phi * ED)
    } else {
      newmus <- pmax(0, mu0 + gamma_N * N + gamma_phi * ED_max)
    }
  }
  return(list(newlas = newlas, newmus = newmus))
}

edd_get_edmax <- function(N, L, age, metric, offset) {
  if (metric == "ed") {
    ED_max <- as.vector(L2ED(L, age))
  } else if (metric == "pd") {
    if (offset == "none") {
      ED_max <- rep(as.vector(L2Phi(L, age, metric)), N)
    } else if (offset == "simtime") {
      ED_max <- rep(as.vector(L2Phi(L, age, metric) - age), N)
    } else {
      stop("no such offset method")
    }
  } else {
    stop("no such metric")
  }
  return(ED_max)
}

edd_get_ed <- function(N, L, t, metric, offset){
  if (metric == "ed") {
    ED <- as.vector(L2ED(L, t))
  } else if (metric == "pd") {
    if (offset == "none") {
      ED <- rep(as.vector(L2Phi(L, t, metric)), N)
    } else if (offset == "simtime") {
      ED <- rep(as.vector(L2Phi(L, t, metric) - t), N)
    } else {
      stop("no such offset method")
    }
  } else {
    stop("no such metric")
  }
  return(ED)
}

edd_sum_rates <- function(las, mus) {
  return(sum(las) + sum(mus))
}

edd_sample_event <- function(las, mus, linlist) {
  events <- 1:(2 * length(linlist))
  return(DDD::sample2(events, 1, prob = c(las, mus)))
}

#' Function to simulate the evolutionary distinctiveness dependent 
#' diversification process
#' 
#' Simulating the evolutionary distinctiveness dependent 
#' diversification process
#' 
#' 
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to
#' lambda (speciation rate) \cr \code{pars[2]} corresponds to mu (extinction
#' rate) \cr \code{pars[3]} corresponds to beta-N (coefficient for species 
#' number effect on speciation) \cr \code{pars[4]} corresponds to beta-Phi 
#' (coefficient for evolutionary distinctiveness effect on speciation) 
#' \cr \code{pars[5]} corresponds to gamma-N (coefficient for species number 
#' effect on speciation) \cr \code{pars[6]} corresponds to gamma-Phi 
#' (coefficient for evolutionary distinctiveness effect on extinction) 
#' @param age Sets the crown age for the simulation
#' @param model Sets the model of diversity-dependence: \cr \code{model ==
#' dsce2} : linear dependence in speciation rate with parameters beta-N and
#' beta-Phi\cr \code{model == dsde2} : linear dependence
#' in both speciation rate and extinction rate with parameters beta-N, beta-Phi,
#' gamma-N and gamma-Phi
#' @return \item{ out }{ A list with the following four elements: The first
#' element is the tree of extant species in phylo format \cr The second element
#' is the tree of all species, including extinct species, in phylo format \cr
#' The third element is a matrix of all species where \cr - the first column is
#' the time at which a species is born \cr - the second column is the label of
#' the parent of the species; positive and negative values only indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values only indicate whether the species belongs to the left or
#' right crown lineage \cr - the fourth column is the time of extinction of the
#' species. If this equals -1, then the species is still extant.\cr The fourth
#' element is the set of branching times of the tree of extant species.\cr The
#' fifth element is the lineage-through-time plot. \cr The sixth element is a
#' matrix of all evolutionary distinctiveness values of all lineages. \cr The
#' seventh element is a matrix of all the speciation and extinction rates of all
#' lineages. \cr The eighth element is a matrix of}
#' @author Tianjian Qin, Rampal S. Etienne
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#' edd_sim(pars = c(0.5, 0.1,-0.001,-0.001,0.001,0.001), age = 6, 
#' model = "dsde2", metric = "ed", offset = "none")
#' @export edd_sim
edd_sim <- function (pars,
                     age,
                     model = "dsce2",
                     metric = "ed",
                     offset = "none") {
  
  edd_pars_check(pars,age,model,metric,offset)
  
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
    linlist <- c(-1, 2)
    newL <- 2
    params <- c(N, pars)
    ED <- c(0, 0)
    ED_max <- edd_get_edmax(N, L, age, metric, offset)
    lamu <- edd_update_lamu(ED, ED_max, params, model)
    
    # controlling significant digits in tibble objects
    options(pillar.sigfig = 10)
    
    # store EDs and lamus and associated lineages
    EDs <- list(ED)
    lamus <- list(lamu)
    linlists <- list(linlist)
    
    # get time interval
    t[i + 1] <-
      t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))
    
    # main simulation circle
    while (t[i + 1] <= age) {
      # time step index
      i <- i + 1
      ED <- edd_get_ed(N[i - 1], L, t[i], metric, offset)
      lamu_real <- edd_update_lamu(ED, ED, params, model)
      event_type <- sample(c('real', 'fake'),
                      1,
                      prob = c(
                        sum(lamu_real$newlas + lamu_real$newmus),
                        sum(
                          lamu$newlas - lamu_real$newlas + lamu$newmus - lamu_real$newmus
                        )
                      ))
      if (event_type == 'real') {
        event <-
          edd_sample_event(lamu_real$newlas, lamu_real$newmus, linlist)
        ranL <- c(linlist, linlist)[event]
        
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
      
      if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
        t[i + 1] <- Inf
      } else {
        ED <- edd_get_ed(N[i], L, t[i], metric, offset)
        ED_max <- edd_get_edmax(N[i], L, age, metric, offset)
        params[1] <- N[i]
        lamu <- edd_update_lamu(ED, ED_max, params, model)
        
        # append EDs and lamus
        EDs <- c(EDs, list(ED))
        lamus <- c(lamus, list(lamu))
        linlists <-c(linlists, list(linlist))
        
        t[i + 1] <-
          t[i] + stats::rexp(1, edd_sum_rates(lamu$newlas, lamu$newmus))
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
  
  LTT <-
    data.frame("time" = t[-i],
               "N" = N)
  out <-
    list(
      tes = tes,
      tas = tas,
      L = L,
      brts = brts,
      LTT = LTT,
      ED = EDs,
      lamu = lamus,
      linlist = linlists
    )
  
  return(out)
} 