//' @useDynLib DDD


#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "odeint_helper.hpp"


using namespace Rcpp;


class ode_rhs
{
public:
  ode_rhs(NumericVector parsvec)
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
    // lx = length(x)
    // xx = c(0, x, 0)
    // dx = lavec[(2:(lx+1))+kk-1] * nn[(2:(lx+1))+2*kk-1] * xx[(2:(lx+1))-1] 
    //    - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
    //    + muvec[(2:(lx+1))+kk+1] * nn[(2:(lx+1))+1] * xx[(2:(lx+1))+1] 
    // return list(dx)
    
    dx.front() = dx.back() = 0.0;
    const size_t lx = xx.size() - 1;
	  for (size_t i = 1; i < lx; ++i) {
      const size_t i0 = i - 1;
      const size_t i1 = i + 1;
      dx[i] = lavec[i0 + kk] * nn[i0 + 2*kk] * xx[i0]
            + muvec[i1 + kk] * nn[i1] * xx[i1]
            - (lavec[i + kk] + muvec[i + kk]) * nn[i + kk] * xx[i];
    }
  }
  
private:
  size_t kk;
  std::vector<double> lavec;
  std::vector<double> muvec;
  std::vector<double> nn;
};


//' Driver for odeint
//' 
// [[Rcpp::export]]
NumericVector dd_integrate_odeint(NumericVector ry, 
                                  NumericVector times, 
                                  NumericVector pars, 
                                  double atol, 
                                  double rtol,
                                  std::string stepper) 
{
  std::vector<double> y(ry.size() + 2, 0.0);    // [0,y,0]
  std::copy(ry.begin(), ry.end(), y.begin() + 1); 

  auto rhs_obj = ode_rhs(pars);
  odeint_helper::integrate(stepper, std::ref(rhs_obj), y, times[0], times[1], 0.1 * (times[1] - times[0]), atol, rtol);
  return NumericVector(y.cbegin() + 1, y.cend() - 1);
}
