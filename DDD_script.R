methode_vec <- c(
  'analytical',
  'lsoda',
  'ode45',
  'lsodes',
  'odeint::runge_kutta_cash_karp54',
  'odeint::runge_kutta_fehlberg78',    # default
  'odeint::runge_kutta_dopri5',
  'odeint::bulirsch_stoer'
)



DDD_script <- function(i,j)
{
  library("ape")
  library("geiger")
  library("DDD")
  library("picante")
  library("pspline")
  library("TreePar")
  
  #--------------------------------------------------------------------------------------------------------------------------------------##
  ############################################## Initialization ############################33
  #--------------------------------------------------------------------------------------------------------------------------------------##
  
  tree <- read.nexus("5000.BicyclusTrees.tre")
  brts <- sort(branching.times(tree))
  tree$node.label <- ((length(tree$tip)+1):((length(tree$tip)*2)-1))
  sub <- extract.clade(tree, node=99)
  main <- read.nexus("5000.BicyclusTrees.tre")
  total_richness <- Ntip(main) + 8
  main <- drop.tip(main,tip = sub$tip.label)
  brtsM <- as.vector(sort(branching.times(main)))
  brtsS <- as.vector(sort(branching.times(sub)))
  brtsMS <- sort(c(brtsM,brtsS))
  tsplit <- NULL
  for(k in 1:length(brts))
  {
    if(!any(abs(brtsMS - brts[k]) < 1E-8)) tsplit <- c(tsplit,brts[k])
  }
  if(length(tsplit) != 1) stop('Something\'s wrong with the branching times.')
  brtsM <- as.vector(sort(c(tsplit,brtsM)))
  
  pars_M <- bd_ML(brtsM)
  pars_S <- bd_ML(brtsS)
  lambda_M <- pars_M[[1]]
  lambda_S <- pars_S[[1]]
  mu_M <- pars_M[[2]]
  mu_S <- pars_S[[2]]
  K_M <- 500
  t_d <- 19 
  pars_all <- bd_ML(brts)
  lambda <- pars_all[[1]]
  mu <- pars_all[[2]]
    
  tsplit <- 20.1718
  missnumspec <- 7
  res <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
  soc <- 2
  cond <- 1
  correction <- TRUE
  optimmethod <- "simplex"
  num_cycles <- 3
  methode <- methode_vec[j]

  btorph <- 1

  print(methode)
  
  #--------------------------------------------------------------------------------------------------------------------------------------##
  ######################################################### Key Innovation Models ###############################################################
  #--------------------------------------------------------------------------------------------------------------------------------------##
  
  # KI0 = diversity indepedent (DI) null model
  # KI1 = diversity dependent shift in K
  # SR2 = diversity dependent shift in K and  mu
  # SR3 = diversity dependent shift in K and lambda
  # SR4 = diversity dependent shift in K, lambda, and mu
  
  if(i == 1) {
  # KI0 = diversity indepedent (DI) null model
  print("KI_0: null model with dummy shift")
  KI_0a <- dd_KI_ML(brtsM = brtsM,
                    brtsS = brtsS,
                    tsplit = tsplit,
                    initparsopt = c(lambda_M, mu_M,  19.244),
                    parsfix = c(Inf, Inf),
                    idparsopt = c(1,2,7),
                    idparsfix = c(3,6),
                    idparsnoshift = c(4,5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI_0a)
  }
  
  if(i == 2) {
  print("KI_0: null model with dummy shift")
  KI_0b <- dd_KI_ML(brtsM = brtsM,
                    brtsS = brtsS,
                    tsplit = tsplit,
                    initparsopt = c(0.121662, 2.160097e-06, 19.3),
                    parsfix = c(Inf,Inf),
                    idparsopt = c(1,2,7),
                    idparsfix = c(3,6),
                    idparsnoshift = c(4,5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI_0b)
  }
  
  ##### KI1 = diversity dependent shift in K
  if(i == 3)
  {
  print("KI1_a: diversity dependent Key innovation diff in K")
  KI1_a <- dd_KI_ML(brtsM = brtsM,
                    brtsS = brtsS,
                    tsplit = tsplit,
                    initparsopt = c(lambda_M, mu_M, length(brtsM)+8, length(brtsS)+8, 18),
                    idparsopt = c(1:3,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4,5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI1_a)
  }
  
  if(i == 4)
  {
  print("KI1_b: diversity dependent Key innovation diff in K")
  KI1_b <- dd_KI_ML(brtsM = brtsM,
                    brtsS = brtsS,
                    tsplit = tsplit,
                    initparsopt = c(0.288, 1.4e-06, 53.6, 80.5, 19),
                    idparsopt = c(1:3,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4,5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI1_b)
  }
  
  if(i == 5)
  {
  print("KI1_c: diversity dependent Key innovation diff in K")
  KI1_c <- dd_KI_ML(brtsM = brtsM,
                    brtsS=brtsS,
                    tsplit = tsplit,
                    initparsopt = c(0.2786422, 1.52E-06, 54.51665, 75.15848, 17.6),
                    idparsopt = c(1:3,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4,5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI1_c)
  }

  ##### KI2 = diversity dependent shift in K and mu

  if(i == 6)
  {
  print("KI2_a: diversity dependent Key innovation diff in K, mu")
  KI2_a <- dd_KI_ML(brtsM = brtsM,
                    brtsS = brtsS,
                    tsplit = tsplit,
                    initparsopt = c(lambda_M, mu_M, length(brtsM)+8, mu_S, length(brtsS)+8, 17.7),
                    idparsopt = c(1:3,5:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI2_a)
  }
  
  if(i == 7)
  {
  print("KI2_b: diversity dependent Key innovation diff in K, mu")
  KI2_b <- dd_KI_ML(brtsM = brtsM,
                    brtsS = brtsS,
                    tsplit = tsplit,
                    initparsopt = c(0.312, 0.037, 50.6, 5.0E-10, 76.3, 17.7),
                    idparsopt = c(1:3,5:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI2_b)
  }
  
  if(i == 8)
  {
  print("KI2_c: diversity dependent Key innovation diff in K, mu")
  KI2_c <- dd_KI_ML(brtsM = brtsM,
                    brtsS = brtsS,
                    tsplit = tsplit,
                    initparsopt = c(0.2936094, 0.03910276, 52.18853, 5.35E-10, 72.95699, 19.24469),
                    idparsopt = c(1:3,5:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI2_c)
  }

  ##### KI3 = diversity dependent shift in K and lambda

  if(i == 9)
  {
  print("KI3_a: diversity dependent Key innovation diff in K, la")
  KI3_a <- dd_KI_ML(brtsM = brtsM,
                    brtsS = brtsS,
                    tsplit = tsplit,
                    initparsopt = c(lambda_M, mu_M, length(brtsM)+8, lambda_S, length(brtsS)+8, 17.7),
                    idparsopt = c(1:4,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI3_a)
  }
  
  if(i == 10)
  {
  print("KI3_b: diversity dependent Key innovation diff in K, la")
  KI3_b <- dd_KI_ML(brtsM = brtsM,
                    brtsS=brtsS,
                    tsplit = tsplit,
                    initparsopt = c(0.242, 2.0e-12, 59.5, 0.337, 73.1, 17.7),
                    idparsopt = c(1:4,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI3_b)
  }
  
  if(i == 11)
  {
  print("KI3_c: diversity dependent Key innovation diff in K, la")
  KI3_c <- dd_KI_ML(brtsM = brtsM,
                    brtsS=brtsS,
                    tsplit = tsplit,
                    initparsopt = c(0.2277065, 0.000000000002071503, 62.16587, 0.3066715, 71.44503, 19.24469),
                    idparsopt = c(1:4,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI3_c)
  }
  
  ##### KI4 = diversity dependent shift in K, mu and lambda

  if(i == 12)
  {
  print("KI4_a: diversity dependent Key innovation diff in K, la, mu")
  KI4_a <- dd_KI_ML(brtsM = brtsM,
                    brtsS=brtsS,
                    tsplit = tsplit,
                    initparsopt = c(lambda_M, mu_M, length(brtsM)+8, lambda_S, mu_S, length(brtsS)+8, 17.7),
                    idparsopt = c(1:7),
                    idparsfix = NULL,
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI4_a)
  }
  
  if(i == 13)
  {
  print("KI4_b: diversity dependent Key innovation diff in K, la, mu")
  KI4_b <- dd_KI_ML(brtsM = brtsM,
                    brtsS=brtsS,
                    tsplit = tsplit,
                    initparsopt = c(0.241, 2.6e-13, 59.5, 0.410, 0.038, 66, 17.7),
                    idparsopt = c(1:7),
                    idparsfix = NULL,
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI4_b)
  }
  
  if(i == 14)
  {
  print("KI4_c: diversity dependent Key innovation diff in K, la, mu")
  KI4_c <- dd_KI_ML(brtsM = brtsM,
                    brtsS = brtsS,
                    tsplit = tsplit,
                    initparsopt = c(0.2276901, 1.840177e-07, 62.17377, 0.3066758, 2.755141e-18, 71.44633, 19.24469),
                    idparsopt = c(1:7),
                    idparsfix = NULL,
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode,
                    correction = correction)
  print(KI4_c)
  }
  
  #--------------------------------------------------------------------------------------------------------------------------------------##
  ######################################################### shifting Rates Models ###############################################################
  #--------------------------------------------------------------------------------------------------------------------------------------##
  
  # SR1 = diversity dependent shift in K
  # SR2 = diversity dependent shift in K and  mu
  # SR3 = diversity dependent shift in K and lambda
  # SR4 = diversity dependent shift in K, lambda, and mu

  res <- 10 * (1 + length(brts) + sum(missnumspec))
    
  if(i == 15)
  {
  print("SR_1: diversity dependent Key shift in K")
  SR_1a <- dd_SR_ML(brts = brts,
                    initparsopt = c(lambda, mu, length(brts)+8, 140, 15), 
                    idparsopt = c(1:3,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4,5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_1a)
  }
  
  if(i == 16)
  {
  print("SR_1: diversity dependent Key shift in K")
  SR_1b <- dd_SR_ML(brts = brts,
                    initparsopt = c(0.8, 0.1, 200, 140, 18), 
                    idparsopt = c(1:3,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4,5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_1b)
  }
  
  if(i == 17)
  {
  print("SR_1: diversity dependent Key shift in K")
  SR_1c <- dd_SR_ML(brts = brts,
                    initparsopt = c(0.18877, 8.557095e-03, 50, 120.2008, 14), 
                    idparsopt = c(1:3,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4,5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_1c)
  }

  if(i == 17.5)
  {
    print("SR_1: diversity dependent Key shift in K")
    SR_1d <- dd_SR_ML(brts = brts,
                      initparsopt = c(0.3485, 2.6E-06, 29.7, 120.2, 9.6), 
                      idparsopt = c(1:3,6:7),
                      idparsfix = NULL,
                      idparsnoshift = c(4,5),
                      res = res,
                      soc = soc,
                      cond = cond,
                      missnumspec = missnumspec,
                      optimmethod = optimmethod,
                      num_cycles = num_cycles,
                      methode = methode)
    print(SR_1d)
  }

  if(i == 17.6)
  {
    print("SR_1: diversity dependent Key shift in K")
    SR_1e <- dd_SR_ML(brts = brts,
                      initparsopt = c(0.3919175, 5.718238e-09, 52.36168, 116.024, 9.933965),
                      idparsopt = c(1:3,6:7),
                      idparsfix = NULL,
                      idparsnoshift = c(4,5),
                      res = res,
                      soc = soc,
                      cond = cond,
                      missnumspec = missnumspec,
                      optimmethod = optimmethod,
                      num_cycles = num_cycles,
                      methode = 'ode45')
    print(SR_1e)
  }

  if(i == 17.6)
  {
    print("SR_1: diversity dependent Key shift in K")
    SR_1e <- dd_SR_ML(brts = brts,
                      initparsopt = c(0.3919175, 5.718238e-09, 52.36168, 116.024, 9.933965),
                      idparsopt = c(1:3,6:7),
                      idparsfix = NULL,
                      idparsnoshift = c(4,5),
                      res = res,
                      soc = soc,
                      cond = cond,
                      missnumspec = missnumspec,
                      optimmethod = optimmethod,
                      num_cycles = num_cycles,
                      methode = 'lsoda')
    print(SR_1e)
  }

  if(i == 17.6)
  {
    print("SR_1: diversity dependent Key shift in K")
    SR_1f <- dd_SR_ML(brts = brts,
                      initparsopt = c(0.3470108, 6.147171e-10, 38.2134, 120.4255, 9.585393),
                      idparsopt = c(1:3,6:7),
                      idparsfix = NULL,
                      idparsnoshift = c(4,5),
                      res = res,
                      soc = soc,
                      cond = cond,
                      missnumspec = missnumspec,
                      optimmethod = optimmethod,
                      num_cycles = num_cycles,
                      methode = 'ode45')
    print(SR_1f)
  }

  if(i == 17.6)
  {
    print("SR_1: diversity dependent Key shift in K")
    SR_1f <- dd_SR_ML(brts = brts,
                      initparsopt = c(0.3470108, 6.147171e-10, 38.2134, 120.4255, 9.585393),
                      idparsopt = c(1:3,6:7),
                      idparsfix = NULL,
                      idparsnoshift = c(4,5),
                      res = res,
                      soc = soc,
                      cond = cond,
                      missnumspec = missnumspec,
                      optimmethod = optimmethod,
                      num_cycles = num_cycles,
                      methode = 'lsoda')
    print(SR_1f)
  }
  
  if(i == 18)
  {
  print("SR_2: diversity dependent Key shift in K, mu")
  SR_2a <- dd_SR_ML(brts = brts,
                    initparsopt = c(lambda, mu, length(brts)+8, 0.065, 140, 15), 
                    idparsopt = c(1:3,5:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_2a)
  }

  if(i == 19)
  {
  print("SR_2: diversity dependent Key shift in K, mu")
  SR_2b <- dd_SR_ML(brts = brts,
                    initparsopt=c(0.44931, 4.98E-08, 29.44462, 9.03E-09, 118.6032, 10.2358), 
                    idparsopt = c(1:3,5:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_2b)
  }
  
  if(i == 20)
  {
  print("SR_2: diversity dependent Key shift in K, mu")
  SR_2c <- dd_SR_ML(brts = brts,
                    initparsopt = c(0.3484879, 5.37E-08, 29.68422, 1.05E-08, 120.2019, 9.585393), 
                    idparsopt = c(1:3,5:7),
                    idparsfix = NULL,
                    idparsnoshift = c(4),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_2c)
  }
  
  if(i == 21)
  {
  print("SR_3: diversity dependent Key shift in K, lambda")
  SR_3a <- dd_SR_ML(brts = brts,
                    initparsopt = c(lambda, mu, length(brts)+8, 0.49, 140, 15), 
                    idparsopt = c(1:4,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_3a)
  }

  if(i == 22)
  {
  print("SR_3: diversity dependent Key shift in K, lambda")
  SR_3b <- dd_SR_ML(brts = brts,
                    initparsopt = c(0.2484879, 5.37E-08, 29.68422, 0.49, 140, 17), 
                    idparsopt = c(1:4,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_3b)
  }
  
  if(i == 23)
  {
  print("SR_3: diversity dependent Key shift in K, lambda")
  SR_3c <- dd_SR_ML(brts = brts,
                    initparsopt = c(0.3959745, 3.45E-07, 28.08213, 0.3349147, 121.9991, 9.585393), 
                    idparsopt = c(1:4,6:7),
                    idparsfix = NULL,
                    idparsnoshift = c(5),
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_3c)
  }
  
  if(i == 24)
  {
  print("SR_4: diversity dependent Key shift in K, lambda, mu")
  SR_4a <- dd_SR_ML(brts = brts, 
                    initparsopt = c(lambda, mu, length(brts)+8, 0.49, 0.065, 140, 15), 
                    idparsopt = c(1:7),
                    idparsfix = NULL,
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_4a)
  }  
  
  if(i == 25)
  {
  print("SR_4: diversity dependent Key shift in K, lambda, mu")
  SR_4b <- dd_SR_ML(brts = brts, 
                    initparsopt = c(0.3528819,  1.396088e-06, 26.60786, 0.3235725, 1.594448e-05, 123.5547, 16), 
                    idparsopt = c(1:7),
                    idparsfix = NULL,
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_4b)
  }
  
  if(i == 26)
  {
  print("SR_4: diversity dependent Key shift in K, lambda, mu")
  SR_4c <- dd_SR_ML(brts = brts, 
                    initparsopt = c(0.3,  0.01, 40, 0.13, 0.01, 130, 9.5), 
                    idparsopt = c(1:7),
                    idparsfix = NULL, 
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_4c)
  }
  
  if(i == 27)
  {
  print("SR_4: diversity dependent Key shift in K, lambda, mu")
  SR_4d <- dd_SR_ML(brts = brts, 
                    initparsopt=c(1.187448,  0.148156, 21.69782, 0.316318, 5.098066e-17, 124.6615, 10.2358), 
                    idparsopt = c(1:7),
                    idparsfix = NULL,
                    res = res,
                    soc = soc,
                    cond = cond,
                    missnumspec = missnumspec,
                    optimmethod = optimmethod,
                    num_cycles = num_cycles,
                    methode = methode)
  print(SR_4d)
  }
  
  
  #--------------------------------------------------------------------------------------------------------------------------------------##
  ######################################################### MS Models ###############################################################
  #--------------------------------------------------------------------------------------------------------------------------------------##

  res <- 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec))
    
  if(i == 28)
  {
    print("MS_0B: DD null model with dummy shift")
    MS_0B1 <- dd_MS_ML(brtsM = brtsM,
                       brtsS = brtsS,
                       tsplit = tsplit,
                       initparsopt = c(lambda_M, mu_M, K_M, t_d),
                       idparsopt = c(1,2,3,6),
                       idparsfix = NULL,
                       idparsnoshift = c(4,5),
                       res = res,
                       soc = soc,
                       cond = 0,
                       missnumspec = missnumspec,
                       correction = correction,
                       methode = methode,
                       optimmethod = optimmethod,
                       verbose = FALSE,
                       num_cycles = num_cycles)
    print(MS_0B1)
  }
  
  #--------------------------------------------------------------------------------------------------------------------------------------##
  ######################################################### Standard Density Diversity  Models ###############################################################
  #--------------------------------------------------------------------------------------------------------------------------------------##
  
  res <- 10 * (1 + length(brts) + sum(missnumspec))
  
  # CRO = DDD0 = DI;    ddmodel == 1 diversity independent model
  # CR1 = DDD1 = DDL;   ddmodel == 1 speciation declines linearly with diversity and no extinction
  # CR2 = DDD2 = DDL+E; ddmodel == 1 speciation declines linearly with diversity and non-zero extinction
  # CR3 = DDD3 = DDX+E; ddmodel == 2 speciation declines exponentially with diversity and non-zero extinction
  # CR4 = DDD3 = DD+EL; ddmodel == 3 extinction increases linearly with diversity
  # CR5 = DDD3 = DD+EX; ddmodel == 4 extinction increases exponentially with diversity
  
  if(i == 29)
  {
  print("Linear independence of speciation rate with extinction")
  DDD_0a <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.3, 0.1),
                  idparsopt = c(1,2),
                  idparsfix = 3,
                  parsfix = Inf,
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_0a)
  
  DDD_0b <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.1, 0.01),
                  idparsopt = c(1,2),
                  idparsfix = 3,
                  parsfix = Inf,
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_0b)
  
  DDD_0c <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.5, 0.1),
                  idparsopt=c(1,2),
                  idparsfix = 3,
                  parsfix = Inf,
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_0c)
  
  DDD_0d <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.1393492, 1.441418E-17),
                  idparsopt = c(1,2),
                  idparsfix = 3,
                  parsfix = Inf,
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_0d)
  
  ###
  print("Linear dependence of speciation rate without extinction")
  DDD_1a <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.3, total_richness),
                  idparsopt = c(1,3),
                  idparsfix = c(2),
                  parsfix = c(0),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_1a)
  
  DDD_1b <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.2, 500),
                  idparsopt = c(1,3),
                  idparsfix = c(2),
                  parsfix = c(0),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_1b)
  
  DDD_1c <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.1, 500),
                  idparsopt = c(1,3),
                  idparsfix = c(2),
                  parsfix = c(0),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_1c)
  
  DDD_1d <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.1756549, 540.6477),
                  idparsopt = c(1,3),
                  idparsfix = c(2),
                  parsfix = c(0),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_1d)
  
  
  ###
  print("Linear dependence of speciation rate with extinction") #tweak rates from 0.3 and 0.1 to 0.2 and 0.01 respectively
  DDD_2a <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.3, 0.01, total_richness),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_2a)
  
  DDD_2b <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.158, 5.753e-05, 600),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_2b)
  
  DDD_2c <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.18, 0.0015, 500),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_2c)
  
  DDD_2d <- dd_ML(brts = brts,
                  ddmodel = 1,
                  initparsopt = c(0.257, 0.001602569, 1384),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_2d)
  
  ###
  print("Exponential dependence of speciation rate with extinction")
  DDD_3a <- dd_ML(brts = brts,
                  ddmodel = 2,
                  initparsopt = c(0.8, 0.1, 500),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_3a)
  
  DDD_3b <- dd_ML(brts = brts,
                  ddmodel = 2, 
                  initparsopt = c(0.5, 0.1, 500),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_3b)
  
  DDD_3c <- dd_ML(brts = brts,
                  ddmodel = 2,
                  initparsopt = c(0.5, 0.1, total_richness),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_3c)
  
  DDD_3d <- dd_ML(brts = brts,
                  ddmodel = 2,
                  initparsopt = c(0.257, 0.01, 140),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_3d)
  
  ###
  print("Linear dependence of extinction rate") #tweak rates from 0.3 and 0.1 to 0.2 and 0.01 respectively
  DDD_4a <- dd_ML(brts = brts,
                  ddmodel = 3,
                  initparsopt = c(0.1941203, 0.0001050502, 400), 
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_4a)
  
  DDD_4b <- dd_ML(brts = brts,
                  ddmodel = 3,
                  initparsopt = c(0.1551, 9.72E-05, 391),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_4b)
  
  DDD_4c <- dd_ML(brts = brts,
                  ddmodel = 3,
                  initparsopt = c(0.3, 0.0004, 300),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_4c)
  
  ###
  print("Exponential dependence of extinction rate")
  DDD_5a <- dd_ML(brts = brts,
                  ddmodel = 4,
                  initparsopt=c(0.2, 0.01, 900),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_5a)
  
  DDD_5b <- dd_ML(brts = brts,
                  ddmodel = 4,
                  initparsopt = c(0.1, 0.001, total_richness),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_5b)
  
  DDD_5c <- dd_ML(brts = brts,
                  ddmodel = 4,
                  initparsopt=c(0.3, 0.00001, 300),
                  btorph = btorph,
                  res = res,
                  soc = soc,
                  cond = cond,
                  missnumspec = missnumspec,
                  optimmethod = optimmethod,
                  num_cycles = num_cycles,
                  methode = methode)
  print(DDD_5c)
  }
}  


DDD_test_all = function() {
	for (i in 1:1) {
	  print(paste('### i = ',i))
		for (j in 1:length(methode_vec)) {
		  print(paste('### ', methode_vec[j]))
		  DDD_script(i, j)
		}
	  print('### end')
	}
}



#DDD_test_all()
