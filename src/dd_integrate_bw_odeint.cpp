//' @export dd_integrate_bw_odeint


#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <vector>
#include "odeint_helper.hpp"


using namespace Rcpp;


class ode_bw_rhs
{
public:
  ode_bw_rhs(NumericVector parsvec)
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
    kk = static_cast<size_t>(parsvec[parsvec.size() - 1]);
  }
  
  void operator()(const std::vector<double>& xx, std::vector<double>& dx, double /* t */)
  {
    // R code:
    // lx = length(x) - 1
    // xx = c(0, x[1:lx], 0)
    // dx = lavec[(2:(lx+1))+kk] * nn[(2:(lx+1))+2*kk] * xx[(2:(lx+1))+1] 
    //    + muvec[(2:(lx+1))+kk] * nn[(2:(lx+1))] * xx[(2:(lx+1))-1] 
    //    - (c(lavec[(2:(lx))+kk],0) + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
    //
    // dG = x[1 + (kk == 0)]
    // return(list(c(dx,dG)))

    dx.front() = dx.back() = 0.0;
    const size_t lx = xx.size() - 1;
    for (size_t i = 1; i < lx - 2; ++i) {
      dx[i] = lavec[i + kk] * nn[i + 2 * kk] * xx[i + 1]
            + muvec[i + kk] * nn[i] * xx[i - 1]
            - (lavec[i + kk] + muvec[i + kk]) * nn[i + kk] * xx[i];
    }
    const size_t i = lx - 2;
    dx[i] = lavec[i + kk] * nn[i + 2 * kk] * xx[i + 1]
          + muvec[i + kk] * nn[i] * xx[i - 1]
          - (0.0 + muvec[i + kk]) * nn[i + kk] * xx[i];
    dx[lx-1] = xx[1 + (kk == 0)];
  }
  
private:
  size_t kk;
  std::vector<double> lavec;
  std::vector<double> muvec;
  std::vector<double> nn;
};


//' Driver for the boost::odeint solver
//'
//' @name dd_integrate_bw_odeint
RcppExport SEXP dd_integrate_bw_odeint(SEXP ry, SEXP rtimes, SEXP rpars, SEXP ratol, SEXP rrtol, SEXP rstepper) {
  BEGIN_RCPP
    auto y = as<NumericVector>(ry);
    std::vector<double> yy(y.size() + 2, 0.0);    // [0,y,0]
    std::copy(y.cbegin(), y.cend(), yy.begin() + 1); 
    auto times = as<std::vector<double>>(rtimes);
    auto pars = as<NumericVector>(rpars);
    auto atol = as<double>(ratol);
    auto rtol = as<double>(rrtol);
    auto stepper = as<std::string>(rstepper);
    
    auto rhs_obj = ode_bw_rhs(pars);
    auto steps = std::min(10000.0, (times[1] - times[0]) * 10.0);
    odeint_helper::integrate(stepper, std::ref(rhs_obj), yy, times[0], times[1], (times[1] - times[0]) / steps, atol, rtol);
    return Rcpp::NumericVector(yy.cbegin() + 1, yy.cend() -1);
  END_RCPP
}
