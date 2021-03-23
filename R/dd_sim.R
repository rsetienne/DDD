dd_lamuN = function(ddmodel,pars,N)
{
    la = pars[1]
    mu = pars[2]
    K = pars[3]
    n0 = (ddmodel == 2 | ddmodel == 4)
    if(ddmodel == 1)
    {
        # linear dependence in speciation rate
        laN = max(0,la - (la - mu) * N/K)
        muN = mu
    }
    if(ddmodel == 1.3)
    {
        # linear dependence in speciation rate
        laN = max(0,la * (1 - N/K))
        muN = mu
    }
    if(ddmodel == 2 | ddmodel == 2.1 | ddmodel == 2.2)
    {
        # exponential dependence in speciation rate
        al = (log(la/mu)/log(K+n0))^(ddmodel != 2.2)
        laN = la * (N + n0)^(-al)
        muN = mu
    }
    if(ddmodel == 2.3)
    {
        # exponential dependence in speciation rate
        al = K
        laN = la * (N + n0)^(-al)
        muN = mu
    }
    if(ddmodel == 3)
    {
        # linear dependence in extinction rate
        laN = la
        muN = mu + (la - mu) * N/K
    }
    if(ddmodel == 4 | ddmodel == 4.1 | ddmodel == 4.2)
    {
        # exponential dependence in extinction rate
        al = (log(la/mu)/log(K+n0))^(ddmodel != 4.2)
        laN = la
        muN = mu * (N + n0)^al
    }
    if(ddmodel == 5)
    {
        r = pars[4]
        # linear dependence in speciation rate and extinction rate
        laN = max(0,la - 1/(r+1)*(la-mu) * N/K)
        muN = mu + r/(r+1)*(la-mu)/K * N
    }
    if (ddmodel == 6) {
        alpha <- pars[4]
        # linear dependence in speciation rate, exponential dependence in extinction rate
        al = (log(1 + alpha * (la  - mu) / mu) / log(K))
        laN = max(0, la - (1 - alpha) * (la - mu) * N / K)
        muN = mu * N ^ al
    }
    if (ddmodel == 7) {
        alpha <- pars[4]
        # exponential dependence in speciation rate and extinction rate
        al1 = log(la / alpha * (la - mu) + mu) / log(K)
        al2 = (log(1 + alpha * (la  - mu) / mu) / log(K))
        laN = la * N ^ (-al1)
        muN = mu * N ^ al2
    }
    if (ddmodel == 8) {
        alpha <- pars[4]
        # exponential dependence in speciation rate, linear dependence in extinction rate
        al = log(la / alpha * (la - mu) + mu) / log(K)
        laN = la * N ^ (-al)
        muN =  mu + alpha * (la - mu) * N / K
    }  
    if (ddmodel == 9) {
        # alternative exponential dependence in speciation rate
        laN = la * (mu / la) ^ (N / K)
        muN = mu
    } 
    if (ddmodel == 10) {
        # alternative exponential dependence in extinction rate
        laN = la
        muN = mu * (la / mu) ^ (N / K)
    } 
    if (ddmodel == 11) {
        alpha <- pars[4]
        # linear dependence in speciation rate, alternative exponential dependence in extinction rate
        laN = max(0, la - (1 - alpha) * (la - mu) * N / K)
        muN =  mu + alpha * (la - mu) * N / K
    }
    if (ddmodel == 12) {
        alpha <- pars[4]
        # alternative exponential dependence in speciation rate and extinction rate
        laN =  la * (alpha + (1 - alpha) * mu / la) ^ (N / K)
        muN = mu * ((1 - alpha) + alpha * la / mu) ^ (N / K)
    }
    if (ddmodel == 13) {
        alpha <- pars[4]
        # alternative exponential dependence in speciation rate, linear dependence in extinction rate
        laN =  la * (alpha + (1 - alpha) * mu / la) ^ (N / K)
        muN =  mu + alpha * (la - mu) * N / K
    }
    return(c(laN, muN))
}

#' Function to simulate the diversity-dependent diversification process
#' 
#' Simulating the diversity-dependent diversification process
#' 
#' 
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to
#' lambda0 (speciation rate) \cr \code{pars[2]} corresponds to mu0 (extinction
#' rate) \cr \code{pars[3]} corresponds to K (clade-level carrying capacity) \cr
#' \code{pars[4]} is only relevant for DD models 5 through 8 and 11 through 13 
#' (models where both speciation and extinction are diversity-dependent, see \code{ddmodel} below), 
#' and controls where the speciation and extinction rate intersect at N = K.
#' For \code{ddmodel = 5} (linear dependence in both rates), \code{pars[4]} is r, the
#' ratio of the slopes of the extinction and speciation rate. In all other models,
#' both rates intersect at a weighted mean of lambda0 and mu0, so that 
#' `lambda(K) = mu(K) = alpha * lambda0 + (1 - alpha) * mu0`. In that case, 
#' \code{pars[4]} is alpha. Note that r = alpha / (1 - alpha).
#' @param age Sets the crown age for the simulation
#' @param ddmodel Sets the model of diversity-dependence: \code{ddmodel ==
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
#' \code{ddmodel == 5} : linear dependence in speciation and extinction rate \cr 
#' \code{ddmodel == 6} : linear dependence in speciation rate, exponential dependence in extinction rate \cr
#' \code{ddmodel == 7} : exponential dependence in speciation and extinction rate \cr
#' \code{ddmodel == 8} : exponential dependence in speciation rate, linear dependence in extinction rate \cr
#' DD models 9, 10, 11, 12, 13 are equivalent to models 2, 3, 6, 7, 8 respectively,
#' but with an alternative formulation of the exponential dependence. 
#' In DD models 2, 3, 6, 7, and 8 the exponential model is a power function of N,
#' while in DD models 9, 10, 11, 12, and 13 the exponential model is based on an
#' exponential function of N.
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
#' element is the set of branching times of the tree of extant species.\cr }
#' @author Rampal S. Etienne
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#'  dd_sim(c(0.2,0.1,20),10) 
#' @export dd_sim
dd_sim = function(pars,age,ddmodel = 1)
{
    # Simulation of diversity-dependent process
    #  . start from crown age
    #  . no additional species at crown node
    #  . no missing species in present
    # pars = [la mu K]
    #  . la = speciation rate per species
    #  . mu = extinction rate per species
    #  . K = diversification carrying capacity
    #  . r = ratio of diversity-dependence in extinction rate over speciation rate
    # age = crown age
    # ddmodel = mode of diversity-dependence
    #  . ddmodel == 1 : linear dependence in speciation rate with parameter K
    #  . ddmodel == 1.3: linear dependence in speciation rate with parameter K'
    #  . ddmodel == 2 : exponential dependence in speciation rate
    #  . ddmodel == 2.1: variant with offset at infinity
    #  . ddmodel == 2.2: 1/n dependence in speciation rate
    #  . ddmodel == 2.3: exponential dependence in speciation rate with parameter x
    #  . ddmodel == 3 : linear dependence in extinction rate
    #  . ddmodel == 4 : exponential dependence in extinction rate
    #  . ddmodel == 4.1: variant with offset at infinity
    #  . ddmodel == 4.2: 1/n dependence in speciation rate
    #  . ddmodel == 5 : linear dependence in speciation rate and in extinction rate
    
    done = 0
    while(done == 0)
    {
        # number of species N at time t
        # i = index running through t and N
        t = rep(0,1)
        L = matrix(0,2,4)
        i = 1
        t[1] = 0
        N = 2
        # L = data structure for lineages,
        # . L[,1] = branching times
        # . L[,2] = index of parent species
        # . L[,3] = index of daughter species
        # . L[,4] = time of extinction
        # j = index running through L
        L[1,1:4] = c(0,0,-1,-1)
        L[2,1:4] = c(0,-1,2,-1)
        linlist = c(-1,2)
        newL = 2
        ff = dd_lamuN(ddmodel,pars,N[i])
        laN = ff[1]
        muN = ff[2]
        denom = (laN + muN) * N[i]
        t[i + 1] = t[i] + stats::rexp(1,denom)
        while(t[i + 1] <= age)
        {
            i = i + 1
            ranL = sample2(linlist,1)
            if((laN * N[i - 1] / denom) >= stats::runif(1))
            {
                # speciation event
                N[i] = N[i - 1] + 1
                newL = newL + 1
                L = rbind(L,c(t[i],ranL,sign(ranL) * newL,-1))
                linlist = c(linlist,sign(ranL) * newL)
            } else {
                # extinction event
                N[i] = N[i - 1] - 1
                L[abs(ranL),4] = t[i]
                w = which(linlist == ranL)
                linlist = linlist[-w]
                linlist = sort(linlist)
            }
            if(sum(linlist < 0) == 0 | sum(linlist > 0) == 0)
            {
                t[i + 1] = Inf
            } else {
                ff = dd_lamuN(ddmodel,pars,N[i])
                laN = ff[1]
                muN = ff[2]
                denom = (laN + muN) * N[i]
                if (denom == 0) {
                    t[i + 1] = Inf # no more event will happen -> conclude sim
                } else {
                    t[i + 1] = t[i] + stats::rexp(1,denom)
                }
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
    brts = L2brts(L,dropextinct = T)
    out = list(tes = tes,tas = tas,L = L,brts = brts)
    return(out)
    
}
