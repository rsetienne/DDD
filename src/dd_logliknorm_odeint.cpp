#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <thread>
#include <chrono>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "odeint_helper.hpp"


using namespace Rcpp;
namespace ublas = boost::numeric::ublas;
using matrix_t = ublas::matrix<double>;
using mproxy_t = ublas::matrix_range<matrix_t>;
using cproxy_t = ublas::matrix_range<const matrix_t>;


class logliknorm1_rhs
{
public:
  logliknorm1_rhs(NumericVector parsvec) :
  m(parsvec.cbegin(), parsvec.cend())
  {
  }
  
  void operator()(const std::vector<double>& xx, std::vector<double>& dx, double /* t */)
  {
    // R code
    // nx = length(x)
    // m1 = m[1:nx]
    // m2 = m[(nx+1):(2*nx)]
    // m3 = m[(2*nx+1):(3*nx)]
    // xx = c(0,x,0)
    // dx = m1 * xx[1:nx] + m2 * xx[3:(nx+2)] - m3 * xx[2:(nx+1)]
    dx.front() = dx.back() = 0.0;
    const auto nx = xx.size() - 2;
    for (size_t i = 1; i <= nx; ++i) {
      const auto m1 = m[i - 1];
      const auto m2 = m[i - 1 + nx];
      const auto m3 = m[i - 1 + 2 * nx];
      dx[i] = m1 * xx[i - 1] + m2 * xx[i + 1] - m3 * xx[i];
    }
  }
  
private:
  std::vector<double> m;
};



// matrix copy
template <typename SRC, typename DST>
void mcopy(const SRC& src, DST& dst, size_t dim)
{
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      dst(i, j) = src(i, j);
    }
  }
}


matrix_t as(const NumericMatrix& nm)
{
  const auto dim = nm.ncol();
  matrix_t m(dim, dim);
  mcopy(nm, m, dim);
  return m;
}


mproxy_t mshift(matrix_t& m, size_t dim, size_t i, size_t j)
{
  return mproxy_t(m, ublas::range(i, dim + i), ublas::range(j, dim + j));
}


cproxy_t cshift(const matrix_t& m, size_t dim, size_t i, size_t j)
{
  return cproxy_t(m, ublas::range(i, dim + i), ublas::range(j, dim + j));
}


class logliknorm2_rhs
{
public:
  logliknorm2_rhs(size_t Dim, List pars) :
    dim(Dim),
    m1(as(pars[0])),
    m2(as(pars[1])),
    m3(as(pars[2])),
    m4(as(pars[3])),
    m5(as(pars[4])),
    m6(as(pars[5]))
  {
  }
  
  // xx: [dim+2, dim+2]
  // y:  [dim+2, dim+2]
  void operator()(const matrix_t& xx, matrix_t& dx, double t)
  {
    // nx = sqrt(length(x))
    // dim(x) = c(nx,nx)
    // xx = matrix(0,nx+2,nx+2)
    // xx[2:(nx+1),2:(nx+1)] = x
    // dx = m[[1]] * xx[1:nx,2:(nx+1)] + m[[2]] * xx[3:(nx+2),2:(nx+1)] + m[[4]] * xx[2:(nx+1),1:nx] + m[[5]] * xx[2:(nx+1),3:(nx+2)] - (m[[3]] + m[[6]]) * xx[2:(nx+1),2:(nx+1)]
    // dim(dx) = c(nx^2,1)
    // return(list(dx))
    dx = ublas::zero_matrix<double>(dim + 2, dim + 2);
    for (size_t i = 0; i < dim; i++) {
      for (size_t j = 0; j < dim; j++) {
        dx(i + 1, j + 1) = m1(i, j) * xx(i, j + 1)
                         + m2(i, j) * xx(i + 2, j + 1)
                         + m4(i, j) * xx(i + 1, j)
                         + m5(i, j) * xx(i + 1, j + 2)
                         - (m3(i, j) + m6(i, j)) * xx(i + 1, j + 1);
      }
    }
  }
  
private:
  const size_t dim;
  const matrix_t m1, m2, m3, m4, m5, m6;
};



// [[Rcpp::export]]
NumericVector dd_logliknorm1_odeint(NumericVector ry, 
                                    NumericVector times,
                                    NumericVector pars,
                                    double atol,
                                    double rtol,
                                    std::string stepper) 
{
  std::vector<double> y(ry.length() + 2, 0.0);    // [0,y,0]
  std::copy(ry.cbegin(), ry.cend(), y.begin() + 1); 

  auto rhs_obj = logliknorm1_rhs(pars);
  odeint_helper::integrate(stepper, std::ref(rhs_obj), y, times[0], times[1], 0.1 * (times[1] - times[0]), atol, rtol);
  return Rcpp::NumericVector(y.cbegin() + 1, y.cend() - 1);
}



// [[Rcpp::export]]
NumericMatrix dd_logliknorm2_odeint(NumericMatrix ry, 
                                    NumericVector times,
                                    List pars, 
                                    double atol, 
                                    double rtol, 
                                    std::string stepper) {
  BEGIN_RCPP
    size_t dim = std::sqrt(ry.length());
    matrix_t y(dim + 2, dim + 2, 0.0);
    auto py = mshift(y, dim, 1, 1);
    mcopy(ry, py, dim);

    auto rhs_obj = logliknorm2_rhs(dim, pars);
    odeint_helper::integrate(stepper, std::ref(rhs_obj), y, times[0], times[1], 0.1 * (times[1] - times[0]), atol, rtol);

    NumericMatrix ret(dim, dim);
    mcopy(py, ret, dim);
    return ret;
  END_RCPP
}


