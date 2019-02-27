dd_SR_sim = function(pars,age,ddmodel = 1)
{
# Simulation of diversity-dependent process with a rate shift
#  . start from crown age
#  . no additional species at crown node
#  . no missing species in present
# pars = [la1 mu1 K1 la2 mu2 K2 tshift]
#  . la1 = speciation rate per species before the shift
#  . mu1 = extinction rate per species before the shift
#  . K1 = diversification carrying capacity before the shift
#  . la2 = speciation rate per species after the shift
#  . mu2 = extinction rate per species after the shift
#  . K2 = diversification carrying capacity after the shift
#  . tshift = time of shift in Mya
# age = crown age in Mya
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
    ff = dd_lamuN(ddmodel,pars[1:3],N[i])
    laN = ff[1]
    muN = ff[2]
    denom = (laN + muN) * N[i]
    t[i + 1] = t[i] - log(stats::runif(1))/denom
    shifted = 0
    tshift = pars[7]
    while(t[i + 1] <= age | shifted == 0)
    {
        if(t[i + 1] >= age - tshift & shifted == 0)
        {
            ff = dd_lamuN(ddmodel,pars[4:6],N[i])
            laN = ff[1]
            muN = ff[2]
            denom = (laN + muN) * N[i]
            t[i + 1] = age - tshift - log(stats::runif(1))/denom
            shifted = 1
        }
        i = i + 1
        ranL = sample2(linlist,1)
        if((laN * N[i - 1] / denom) >= runif(1))
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
            ff = dd_lamuN(ddmodel,pars[c(1:3) + 3 * (shifted == 1)],N[i])
            laN = ff[1]
            muN = ff[2]
            denom = (laN + muN) * N[i]
            t[i + 1] = t[i] - log(stats::runif(1))/denom
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