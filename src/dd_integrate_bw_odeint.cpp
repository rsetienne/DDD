// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
//' @export dd_integrate_odeint


#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <vector>
#include "odeint_helper.hpp"


dd_loglik_rhs = function(t,x,parsvec)
{
  lv = (length(parsvec) - 1)/3
  lavec = parsvec[1:lv]
  muvec = parsvec[(lv + 1):(2 * lv)]
  nn = parsvec[(2 * lv + 1):(3 * lv)]
  kk = parsvec[length(parsvec)]
  lx = length(x)
  xx = c(0,x,0)
  dx = lavec[(2:(lx+1))+kk-1] * nn[(2:(lx+1))+2*kk-1] * xx[(2:(lx+1))-1] + muvec[(2:(lx+1))+kk+1] * nn[(2:(lx+1))+1] * xx[(2:(lx+1))+1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
  return(list(dx))
}



dd_loglik_bw_rhs = function(t,x,parsvec)
{
#parsvec = c(dd_loglik_rhs_precomp(pars,x),pars[length(pars) - 1])
  lv = (length(parsvec) - 1)/3
  lavec = parsvec[1:lv]
  muvec = parsvec[(lv + 1):(2 * lv)]
  nn = parsvec[(2 * lv + 1):(3 * lv)]
  kk = parsvec[length(parsvec)]
  lx = length(x) - 1
  
  xx = c(0,x[1:lx],0)
  
  dx = lavec[(2:(lx+1))+kk] * nn[(2:(lx+1))+2*kk] * xx[(2:(lx+1))+1] + muvec[(2:(lx+1))+kk] * nn[(2:(lx+1))] * xx[(2:(lx+1))-1] - (c(lavec[(2:(lx))+kk],0) + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]

  dG = x[1 + (kk == 0)]
  
  return(list(c(dx,dG)))
}


class ode_rhs
{
public:
  ode_rhs(std::vector<double> parsvec)
  {
    const size_t lv = (parsvec.size() - 1) / 3;
    lavec.resize(lv, 0);
    muvec.resize(lv, 0);
    nn.resize(lv, 0);
    for (size_t i = 0; i < lv; ++i) {
      lavec[i] = parsvec[i];            // parsvec[1:lv]
      muvec[i] = parsvec[lv + i];       // parsvec[(lv + 1):(2 * lv)]
      nn[i] = parsvec[2 * lv + i];      // parsvec[(2 * lv + 1):(3 * lv)]
    }
    kk = static_cast<size_t>(parsvec.back());
  }
  
  void operator()(const std::vector<double>& x, std::vector<double>& dxdt, double /* t */)
  {
    // R code:
    // lx = length(x) - 1
    // xx = c(0, x[1:lx], 0)
    // dx = lavec[(2:(lx+1))+kk] * nn[(2:(lx+1))+2*kk] * xx[(2:(lx+1))+1] + muvec[(2:(lx+1))+kk] * nn[(2:(lx+1))] * xx[(2:(lx+1))-1] - (c(lavec[(2:(lx))+kk],0) + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
    //
    // dx = lavec[(2:(lx+1))+kk-1] * nn[(2:(lx+1))+2*kk-1] * xx[(2:(lx+1))-1] + muvec[(2:(lx+1))+kk+1] * nn[(2:(lx+1))+1] * xx[(2:(lx+1))+1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
    //
    // dG = x[1 + (kk == 0)]
    // return(list(c(dx,dG)))
    
    // unroll, avoiding xx
    const size_t lx = x.size();
    double x0 = 0.0;
    for (size_t i = 1; i < lx - 1; ++i) {
      const size_t i0 = i - 1;
      const size_t i1 = i + 1;
      dxdt[i0] = lavec[i0 + kk] * nn[i0 + 2*kk] * x0
        + muvec[i1 + kk] * nn[i1] * x[i]
      - (lavec[i + kk] + muvec[i + kk]) * nn[i + kk] * x[i0];
      x0 = x[i0];
    }
    dxdt[lx - 1] = lavec[lx - 1 + kk] * nn[lx - 1 + 2*kk] * x0
      + muvec[lx - 1 + kk] * nn[lx - 1] * x[lx - 1];
    // - 0.0    
  }
  
private:
  size_t kk;
  std::vector<double> lavec;
  std::vector<double> muvec;
  std::vector<double> nn;
};



using namespace Rcpp;



//' Driver for the boost::odeint solver
//'
//' @name dd_integrate_bw_odeint
RcppExport SEXP dd_integrate_bw_odeint(SEXP ry, SEXP rtimes, SEXP rpars, SEXP ratol, SEXP rrtol, SEXP rstepper) {
  BEGIN_RCPP
    auto y = as<std::vector<double>>(ry);
    auto times = as<std::vector<double>>(rtimes);
    auto pars = as<std::vector<double>>(rpars);
    auto atol = as<double>(ratol);
    auto rtol = as<double>(rrtol);
    auto stepper = as<std::string>(rstepper);
    
    auto rhs_obj = ode_rhs(pars);
    odeint_helper::integrate(stepper, std::ref(rhs_obj), y, times[0], times[1], 0.1 * (times[1] - times[0]), atol, rtol);
  
    return Rcpp::NumericVector(y.cbegin(), y.cend());
  END_RCPP
}
