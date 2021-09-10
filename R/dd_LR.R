#' Bootstrap likelihood ratio test of diversity-dependent diversification model
#' 
#' This function computes the maximum likelihood and the associated estimates
#' of the parameters of a diversity-dependent diversification model for a given
#' set of phylogenetic branching times. It then performs a bootstrap likelihood
#' ratio test of the diversity-dependent (DD) model against the constant-rates
#' (CR) birth-death model. Finally, it computes the power of this test.
#' 
#' The output is a list with 3 elements:
#' 
#' @param brts A set of branching times of a phylogeny, all positive
#' @param initparsoptDD The initial values of the parameters that must be
#' optimized for the diversity-dependent (DD) model: lambda_0, mu and K
#' @param initparsoptCR The initial values of the parameters that must be
#' optimized for the constant-rates (CR) model: lambda and mu
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @param outputfilename The name (and location) of the file where the output
#' will be saved. Default is no save.
#' @param seed The seed for the pseudo random number generator for simulating
#' the bootstrap data
#' @param endmc The number of bootstraps
#' @param alpha The significance level of the test
#' @param plotit Boolean to plot results or not
#' @param res Sets the maximum number of species for which a probability must
#' be computed, must be larger than 1 + length(brts)
#' @param ddmodel Sets the model of diversity-dependence: \cr \code{ddmodel ==
#' 1} : linear dependence in speciation rate with parameter K (= diversity
#' where speciation = extinction)\cr \code{ddmodel == 1.3} : linear dependence
#' in speciation rate with parameter K' (= diversity where speciation = 0)\cr
#' \code{ddmodel == 2} : exponential dependence in speciation rate with
#' parameter K (= diversity where speciation = extinction)\cr \code{ddmodel ==
#' 2.1} : variant of exponential dependence in speciation rate with offset at
#' infinity\cr \code{ddmodel == 2.2} : 1/n dependence in speciation rate\cr
#' \code{ddmodel == 2.3} : exponential dependence in speciation rate with
#' parameter x (= exponent)\cr \code{ddmodel == 3} : linear dependence in
#' extinction rate \cr \code{ddmodel == 4} : exponential dependence in
#' extinction rate \cr \code{ddmodel == 4.1} : variant of exponential
#' dependence in extinction rate with offset at infinity \cr \code{ddmodel ==
#' 4.2} : 1/n dependence in extinction rate with offset at infinity \cr
#' \code{ddmodel == 5} : linear dependence in speciation and extinction rate
#' \cr
#' @param cond Conditioning: \cr cond == 0 : conditioning on stem or crown age
#' \cr cond == 1 : conditioning on stem or crown age and non-extinction of the
#' phylogeny \cr cond == 2 : conditioning on stem or crown age and on the total
#' number of extant taxa (including missing species) \cr cond == 3 :
#' conditioning on the total number of extant taxa (including missing species)
#' \cr Note: cond == 3 assumes a uniform prior on stem age, as is the standard
#' in constant-rate birth-death models, see e.g. D. Aldous & L. Popovic 2004.
#' Adv. Appl. Prob. 37: 1094-1115 and T. Stadler 2009. J. Theor. Biol. 261:
#' 58-66.
#' @param btorph Sets whether the likelihood is for the branching times (0) or
#' the phylogeny (1)
#' @param soc Sets whether stem or crown age should be used (1 or 2)
#' @param tol Sets the tolerances in the optimization. Consists of: \cr reltolx
#' = relative tolerance of parameter values in optimization \cr reltolf =
#' relative tolerance of function value in optimization \cr abstolx = absolute
#' tolerance of parameter values in optimization
#' @param maxiter Sets the maximum number of iterations in the optimization
#' @param changeloglikifnoconv if TRUE the loglik will be set to -Inf if ML
#' does not converge
#' @param optimmethod Method used in optimization of the likelihood. Current
#' default is 'subplex'. Alternative is 'simplex' (default of previous
#' versions)
#' @param methode The method used to solve the master equation, default is
#' 'analytical' which uses matrix exponentiation; alternatively numerical ODE
#' solvers can be used, such as 'odeint::runge_kutta_cash_karp54'. These were used in the
#' package before version 3.1.
#' @return \item{treeCR}{a list of trees generated under the constant-rates
#' model using the ML parameters under the CR model} \item{treeDD}{a list of
#' trees generated under the diversity-dependent model using the ML parameters
#' under the diversity-dependent model} \item{out}{a dataframe with the
#' parameter estimates and maximum likelihoods for diversity-dependent and
#' constant-rates models \code{$model} - the model used to generate the data. 0
#' = unknown (for real data), 1 = CR, 2 = DD \cr \code{$mc} - the simulation
#' number for each model \cr \code{$lambda_CR} - speciation rate estimated
#' under CR \cr \code{$mu_CR} - extinction rate estimated under CR \cr
#' \code{$LL_CR} - maximum likelihood estimated under CR \cr \code{$conv_CR} -
#' convergence code for likelihood optimization; conv = 0 means convergence \cr
#' \code{$lambda_DD1} - initial speciation rate estimated under DD for first
#' set of initial values\cr \code{$mu_DD1} - extinction rate estimated under DD
#' for first set of initial values \cr \code{$K_DD1} - clade-wide
#' carrying-capacity estimated under DD for first set of initial values \cr
#' \code{$LL_DD1} - maximum likelihood estimated under DD for first set of
#' initial values \cr \code{$conv_DD1} - convergence code for likelihood
#' optimization for first set of initial values; conv = 0 means convergence \cr
#' \code{$lambda_DD2} - initial speciation rate estimated under DD for second
#' set of initial values\cr \code{$mu_DD2} - extinction rate estimated under DD
#' for second set of initial values \cr \code{$K_DD2} - clade-wide
#' carrying-capacity estimated under DD for second set of initial values \cr
#' \code{$LL_DD2} - maximum likelihood estimated under DD for second set of
#' initial values \cr \code{$conv_DD2} - convergence code for likelihood
#' optimization for second set of initial values; conv = 0 means convergence
#' \cr \code{$LR} - likelihood ratio between DD and CR } \item{pvalue}{p-value
#' of the test} \item{LRalpha}{Likelihood ratio at the signifiance level alpha}
#' \item{poweroftest}{power of the test for significance level alpha}
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{dd_loglik}}, \code{\link{dd_ML}}
#' @references - Etienne, R.S. et al. 2016. Meth. Ecol. Evol. 7: 1092-1099,
#' doi: 10.1111/2041-210X.12565 \cr - Etienne, R.S. et al. 2012, Proc. Roy.
#' Soc. B 279: 1300-1309, doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B.
#' Haegeman 2012. Am. Nat. 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @export dd_LR
dd_LR = function(             
   brts,
   initparsoptDD,
   initparsoptCR,
   missnumspec,   
   outputfilename = NULL,
   seed = 42,
   endmc = 1000,
   alpha = 0.05,
   plotit = TRUE,
   res = 10 * (1 + length(brts) + missnumspec),
   ddmodel = 1,
   cond = 1,
   btorph = 1,
   soc = 2,
   tol = c(1E-3,1E-4,1E-6),
   maxiter = 2000,
   changeloglikifnoconv = FALSE,
   optimmethod = 'subplex',
   methode = 'analytical'   
   )
{
  if(!is.null(seed))
  {
     set.seed(roundn(seed))
  }
  if(cond > 1)
  {
    cat("Conditioning on number of tips is not possible.\n")
    return(NULL)
  }
  age = max(brts)
  cat("\nEstimating parameters under the constant-rate model ...\n")  
  outCRO = dd_ML(brts = brts,initparsopt = initparsoptCR,idparsopt = 1:2,idparsfix = 3,parsfix = Inf,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod, methode = methode)
  cat("\nEstimating parameters under the diversity-dependent model ...\n")  
  outDDO = dd_ML(brts = brts,initparsopt = initparsoptDD,idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod, methode = methode)
  LRO = outDDO$loglik - outCRO$loglik
  out = cbind(NA,NA,outCRO,outDDO,NA,NA,NA,NA,NA,LRO)
  out = out[,-c(5,7,13)]
  newnames = c("model","mc","lambda_CR","mu_CR","LL_CR","conv_CR","lambda_DD1","mu_DD1","K_DD1","LL_DD1","conv_DD1","lambda_DD2","mu_DD2","K_DD2","LL_DD2","conv_DD2","LR")
  names(out) = newnames
  if(!is.null(outputfilename))
  {
      save(seed,brts,out,file = outputfilename)
  }
  parsCR = as.numeric(outCRO[1:2])
  parsDD = as.numeric(outDDO[1:3])
  treeCR = list()
  treeDD = list()
  cat('\nSimulating trees under CR and DD models ...\n')
  for(mc in 1:endmc)
  {
     treeCR[[mc]] = dd_sim(pars = c(parsCR,Inf),age = age,ddmodel = ddmodel)
     treeDD[[mc]] = dd_sim(pars = parsDD,age = age,ddmodel = ddmodel)
  }
  if(!is.null(outputfilename))
  {
      save(seed,brts,out,treeCR,treeDD,file = outputfilename)
  }
  cat('\nPerforming bootstrap to determine critical LR ...\n')  
  for(mc in 1:endmc)
  {
     cat('\nAnalyzing simulation:',mc,'\n')
     brtsCR = ape::branching.times(treeCR[[mc]][[1]])
     outCR = dd_ML(brtsCR,initparsopt = parsCR,idparsopt = 1:2,idparsfix = 3,parsfix = Inf,res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod, methode = methode)
     outDD1 = dd_ML(brtsCR,initparsopt = parsDD,idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod, methode = methode)
     outDD2 = dd_ML(brtsCR,initparsopt = c(parsCR + 0.05,length(brts) + 1000),idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod, methode = methode)
     if(outDD1$conv == -1 & outDD2$conv == -1)
     {
        maxLLDD = outCR$loglik
     } else if(outDD1$conv == -1 & outDD2$conv != -1)
     {
        maxLLDD = outDD2$loglik
     } else if(outDD1$conv != -1 & outDD2$conv == -1)
     {
        maxLLDD = outDD1$loglik
     } else {
        maxLLDD = max(outDD1$loglik,outDD2$loglik)
     }   
     LR = pmax(0,maxLLDD - outCR$loglik)
     outff = cbind(1,mc,outCR,outDD1,outDD2,LR)
     outff = outff[,-c(5,7,13,19)]     
     names(outff) = newnames
     out = rbind(out,outff)
     if(!is.null(outputfilename))
     {
        save(seed,brts,out,treeCR,treeDD,file = outputfilename)
     }
  }
  opt = rep(0,endmc)
  cat('\nPerforming bootstrap to determine power ...\n')
  for(mc in 1:endmc)
  {
     cat('\nAnalyzing simulation:',mc,'\n')
     brtsDD = ape::branching.times(treeDD[[mc]][[1]])
     outCR = dd_ML(brtsDD,initparsopt = parsCR,idparsopt = 1:2, idparsfix = 3,parsfix = Inf,res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod, methode = methode)
     outDD1 = dd_ML(brtsDD,initparsopt = parsDD,idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod, methode = methode)
     outDD2 = dd_ML(brtsDD,initparsopt = c(parsCR + 0.05,length(brts) + 1000), idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod, methode = methode)
     if(outDD1$conv == -1 & outDD2$conv == -1)
     {
        maxLLDD = outCR$loglik
        opt[mc] = 1
     } else if(outDD1$conv != -1 & outDD2$conv == -1)
     {
        maxLLDD = outDD1$loglik
        opt[mc] = 2
     } else if(outDD1$conv == -1 & outDD2$conv != -1)
     {
        maxLLDD = outDD2$loglik
        opt[mc] = 3
     } else {
        maxLLDD = max(outDD1$loglik,outDD2$loglik)
        opt[mc] = 1 + min(which(c(outDD1$loglik,outDD2$loglik) == maxLLDD))
     }
     LR = pmax(0,maxLLDD - outCR$loglik)
     outff = cbind(1,mc,outCR,outDD1,outDD2,LR)
     outff = outff[,-c(5,7,13,19)]
     names(outff) = newnames
     out = rbind(out,outff)
     if(!is.null(outputfilename))
     {
        save(seed,brts,out,opt,treeCR,treeDD,file = outputfilename)
     }
  }
  inverse_quantile = function(samples,x)
  {
     samplessort = sort(samples)
     pup = which(samplessort > x)
     if(length(pup) > 0)
     {
        if(length(pup) < length(samplessort))
        {
           pup = min(pup)
           invquant = (pup + (x - samplessort[pup])/(samplessort[pup - 1] - samplessort[pup]))/length(samples)
        } else {
           invquant = 0
        }
     } else {
        invquant = 1
     }
     return(invquant)
  }
  funpvalue = function(samples,x)
  {
     samplessort = sort(samples)
     pup = which(samplessort > x)
     pvalue = (length(pup) + 1)/ (length(samples) + 1)
     return(pvalue)
  }
  funpoweroftest = function(samples,x)
  {
     samplessort = sort(samples)
     pup = which(samplessort > x)
     poweroftest = length(pup)/(length(samples) + 1)
     return(poweroftest)     
  }
  #pvalue = 1 - inverse_quantile(out$LR[2:(endmc + 1)],out$LR[1])
  pvalue = funpvalue(out$LR[2:(endmc + 1)],out$LR[1])
  LRalpha = as.numeric(stats::quantile(out$LR[2:(endmc + 1)],1 - alpha,type = 4))
  #poweroftest = 1 - inverse_quantile(out$LR[(endmc + 2):(2 * endmc + 1)],LRalpha)
  poweroftest = funpoweroftest(out$LR[(endmc + 2):(2 * endmc + 1)],LRalpha)
  if(plotit == TRUE)
  {
      try(grDevices::dev.off())
      try(grDevices::dev.off())
      pdffilename = paste(getwd(),'/LR.pdf',sep = '')
      grDevices::pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
      al = 0.03
      alw = 2
      alw2 = 1.7
      aa = 45
      graphics::par(mfrow = c(2,2),cex = 1, mar = c(5, 4, 3, 1) + 0.1)
      graphics::hist(out$LR[2:(1 + endmc)],main = 'Distribution of LLR under CR',xlab = 'LLR', ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(out$LR[1:(endmc + 1)])))
      graphics::arrows(out$LR[1],-1E+120, x1 = out$LR[1],y1 = 0, length = al, angle = aa,lwd = alw, col = 'black')
      graphics::arrows(LRalpha,-1E+120, x1 = LRalpha,y1 = 0, length = al, angle = aa,lwd = alw, col = 'blue')
      graphics::box()
      graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", col = "white", col.lab = 'white', col.axis = 'white')
      graphics::hist(out$LR[(endmc + 2):(1 + 2 * endmc)], main = 'Distribution of LLR under DD',xlab = 'LLR', ylab = 'Frequency', col = 'red',probability = T,nclass = 30)
      graphics::box()
      graphics::arrows(out$LR[1],-1E+120, x1 = out$LR[1],y1 = 0, length = al, angle = aa,lwd = alw, col = 'black')
      graphics::arrows(LRalpha,-1E+120, x1 = LRalpha,y1 = 0, length = al, angle = aa,lwd = alw, col = 'blue')

      graphics::par(mfrow = c(2,3),cex = 1, mar = c(5, 4, 3, 1) + 0.1)
      lambda = out$lambda_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$lambda_DD1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$lambda_DD2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
      mu = out$mu_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$mu_DD1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$mu_DD2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
      K = 1E+120 * (opt == 1) + pmin(1E+120,out$K_DD1[(endmc + 2):(2 * endmc + 1)]) * (opt == 2) + pmin(1E+120,out$K_DD2[(endmc + 2):(2 * endmc + 1)]) * (opt == 3)
      graphics::hist(lambda,main = NULL, xlab = expression(lambda), ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(lambda)))
      graphics::arrows(out$lambda_DD1[1],-1E+120, x1 = out$lambda_DD1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      graphics::box()
      graphics::hist(mu,main = NULL, xlab = expression(mu), ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(mu)))
      graphics::arrows(out$mu_DD1[1],-1E+120, x1 = out$mu_DD1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      graphics::box()
      graphics::hist(K,main = NULL, xlab = 'K', ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(min(K),max(K)))
      graphics::arrows(out$K_DD1[1],-1E+120, x1 = out$K_DD1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      graphics::box()
      try(grDevices::dev.off())
      try(grDevices::dev.off())
      os = .Platform$OS.type
      if(os == "windows")
      {
          shell.exec(pdffilename)
      }
      if(os == "unix")
      {
          system(paste("open",pdffilename,sep = " "))
      }
  }
  if(!is.null(outputfilename))
  {
      save(seed,brts,out,opt,treeCR,treeDD,pvalue,LRalpha,poweroftest,file = outputfilename)
  }
  return(list(treeCR = treeCR,treeDD = treeDD,out = out,pvalue = pvalue,LRalpha = LRalpha,poweroftest = poweroftest))
}
