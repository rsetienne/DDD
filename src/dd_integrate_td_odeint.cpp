#define STRICT_R_HEADERS
#include "config.h"
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "odeint_helper.h"


using namespace Rcpp;


class ode_td_rhs
{
  using quad_t = boost::multiprecision::cpp_bin_float_quad;
  
public:
  ode_td_rhs(NumericVector parsvec)
  {
    lp = static_cast<size_t>(parsvec[parsvec.size() - 1]);
    double la = parsvec[0];
    mu = parsvec[1];
    double K = parsvec[2];
    lavec.resize(lp + 2, 0.0);
    for (size_t i = 0; i < lp + 2; ++i) {
      lavec[i] = std::max(0.0, la - (la - mu) / K * (static_cast<double>(i) - 1.0));
    }
  }
  
  void operator()(const std::vector<double>& xx, std::vector<double>& dx, double /* t */)
  {
    // R code
    // nn = -1:lp
    // muvec = mu
    // p = c(0, x, 0)
    // dp = lavec[(2:(lp+1))-1] * nn[(2:(lp + 1))-1] * p[(2:(lp + 1))-1] + muvec[(2:(lp+1))+1] * nn[(2:(lp+1))+1] * p[(2:(lp+1))+1] - (lavec[(2:(lp+1))] + muvec[(2:(lp+1))]) * nn[(2:(lp+1))] * p[(2:(lp+1))]

    p.clear();
    p.push_back(0.0);  
    for (size_t i = 0; i < lp; ++i) p.push_back(xx[i]);
    p.push_back(0.0);

    for (size_t i = 1; i <= lp; ++i) {
      const auto nn = static_cast<double>(i);
      dx[i - 1] = lavec[i - 1] * (nn - 2.0) * p[i - 1]
                + mu * nn * p[i + 1]
                - (lavec[i] + mu) * (nn - 1.0) * p[i];
    }

    // mutd = rep(mu,lrs)
    // En = sum((0:(lp - 1)) * x[1:lp] )
    // dsigdiv = mutd / En
    // return(list(c(dp,dsigdiv)))  
    quad_t En{0};
    for (size_t i = 0; i < lp; ++i) {
      En += quad_t{static_cast<double>(i) * xx[i]};
    }    
    for (size_t i = lp; i < xx.size(); ++i) {
      dx[i] = (quad_t{mu} / En).convert_to<double>();
    }
  }
  
private:
  size_t lp;    // number of probabilities
  double mu;
  std::vector<double> lavec;
  std::vector<double> p;
};



// [[Rcpp::export]]
std::vector<double> dd_integrate_td_odeint(std::vector<double> y, 
                                           NumericVector times, 
                                           NumericVector pars, 
                                           double atol, 
                                           double rtol,
                                           std::string stepper) 
{
  auto rhs_obj = ode_td_rhs(pars);
  odeint_helper::integrate(stepper, std::ref(rhs_obj), y, times[0], times[1], 0.1 * (times[1] - times[0]), atol, rtol);
  return y;
}
