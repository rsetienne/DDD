#' Simulation of a diversity-dependent-like time-dependent process
#' 
#' Simulates a phylogenetic tree branching according to a 
#' time-dependent process calibrated on the expected number 
#' of species under a diversity-dependent process over time.
#' 
#' @param pars Vector of parameters:
#' \cr \cr \code{pars[1]} corresponds to
#' lambda0 (speciation rate)
#' \cr \code{pars[2]} corresponds to mu0 (extinction
#' rate)
#' \cr \code{pars[3]} corresponds to lambda1 (decline parameter in
#' speciation rate) or K in diversity-dependence-like models
#' \cr \code{pars[4]} corresponds to mu1 (decline parameter in extinction rate)
#' @param age crown age of the tree to simulate, i.e. the simulation time.
#' @param ddmodel the diversity-dependent model used as reference for the 
#' time-dependent model.
#' @param methode The method used to solve the master equation. 
#' See \code{deSolve::ode()} documentation for possible inputs
#' 
#' @return A list with the following four elements: The first element is the 
#' tree of extant species in phylo format \cr
#' The second element is the tree of all species, including extinct species, 
#' in phylo format \cr
#' The third element is a matrix of all species where 
#' - the first column is the time at which a species is born \cr
#' - the second column is the label of the parent of the species; 
#' positive and negative values only indicate whether the species belongs to 
#' the left or right crown lineage \cr
#' - the third column is the label of the daughter species itself; 
#' positive and negative values only indicate whether the species belongs to 
#' the left or right crown lineage \cr
#' - the fourth column is the time of extinction of the species. If this 
#' equals -1, then the species is still extant.
#'
#' @author César Martinez, Rampal S. Etienne
#' 
#' @export

td_sim = function(pars, age, ddmodel = 1, methode = 'ode45')
  # Other methodes: 'ode45', 'lsodes'
{
  
  la0 = pars[1]
  mu0 = pars[2]
  K = pars[3]
  
  if(la0 < mu0) {stop('the function is designed for lambda_0 > mu_0')}
  if(mu0 < 0) {stop('per species rates should be positive')}
  if(K < 1) {stop('clade level carrying capacity should be positive')}
  
  Kprime = la0 / (la0 - mu0) * K
  
  soc = 2
  nbeq = 1000
  tdmodel = 4
  
  done = 0
  while(done == 0)
  {
    # number of species N at time t
    # i = index running through t and N
    t = rep(0, 1)  
    L = matrix(0, 2, 4)
    
    i = 1 # counter variable for the potential events    
    ev = 1 #counter variable for the events events (conditioning is at crown age)
    t[1] = 0    #initialisation of a time vector for all potential events
    N = 2
    
    lx = min(1 + ceiling(Kprime), nbeq)
    variables = rep(0, lx)
    variables[soc + 1] = 1
    
    # lambda(t) at t=0, time dependent model
    latd = mu0 + ((la0 - mu0) * soc - soc ^ 2 * (la0-mu0) / K) / soc 
    mutd = mu0
    
    #firt upper values for the rates
    lamax = latd 
    mumax = mu0             
    # L = data structure for lineages,
    # . L[,1] = branching times
    # . L[,2] = index of parent species
    # . L[,3] = index of daughter species
    # . L[,4] = time of extinction
    # j = index running through L
    L[1, 1:4] = c(0, 0, -1, -1)
    L[2, 1:4] = c(0, -1, 2, -1)
    linlist = c(-1, 2)
    newL = 2
    
    # parameter next potential event
    denom = (lamax + mumax) * N[i]
    
    #update the time with the first potential event
    t[i + 1] = t[i] + stats::rexp(1, denom)
    
    while(t[i + 1] <= age){
      i = i + 1     
      
      # sample an individual 
      ranL = sample(linlist, 1)
      
      # calculation of the rate latd at the new time step
      y = deSolve::ode(
        variables,
        c(t[i-1], t[i]),
        td_loglik_rhs_sim,
        c(pars[1:min(4, length(pars))], tdmodel - 3, lx),
        rtol = 1e-10,
        atol = 1e-16,
        method = methode
      )
      
      variables = y[2, 2:(lx + 1)]
      expn = sum((0:(lx - 1)) * variables[1:lx])
      expn2 = sum((0:(lx - 1)) ^ 2 * variables[1:lx])
      dEN_dt = (la0 - mu0) * expn - expn2 * (la0 - mu0) / K
      
      latd = mu0 + dEN_dt / expn
      
      if(lamax - latd < (-1e-10))
      {
        stop('latd should be a decreasing function of time')
      }      
      if( ((latd + mutd) / (lamax + mumax)) >= stats::runif(1) )  #does the next potential event occur?                    
      {
        ev = ev + 1
        lamax = latd
        
        if( latd / (latd + mutd) >= stats::runif(1) )  # Is it a speciation event?
        {           
          N[ev] = N[ev - 1] + 1
          newL = newL + 1
          L = rbind(L,c(t[i],ranL,sign(ranL) * newL,-1))  
          linlist = c(linlist,sign(ranL) * newL)
        } else {         
          N[ev] = N[ev - 1] - 1
          L[abs(ranL),4] = t[i]
          w = which(linlist == ranL)
          linlist = linlist[-w]
          linlist = sort(linlist)          
        }
      } 
      
      # Is one of the two crown lineages extinct ?
      if(sum(linlist < 0) == 0 | sum(linlist > 0) == 0)
      {
        t[i + 1] = Inf
      } 
      else 
      {
        #parameter for the next potential event is updated with new lambda_max, mu_max, and number of species      
        denom = (lamax + mumax) * N[ev]
        #time is updated by adding a new potential step
        t[i + 1] = t[i] + stats::rexp(1,denom)
      } 
    }
    
    if(sum(linlist < 0) == 0 | sum(linlist > 0) == 0)
    {
      done = 0
    } else {
      done = 1
    }
  }
  
  L[,1] = age - c(L[,1])
  notmin1 = which(L[,4] != -1)
  L[notmin1,4] = age - c(L[notmin1,4])
  L[which(L[,4] == age + 1),4] = -1
  tes = L2phylo(L,dropextinct = T)
  tas = L2phylo(L,dropextinct = F)
  out = list(tes,tas,L)
  return(out)  
}
