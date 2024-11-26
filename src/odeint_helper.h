// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
#include "config.h"
#ifndef ODEINT_HELPER_HPP_INCLUDED
#define ODEINT_HELPER_HPP_INCLUDED


/*

  class MyOde
  {
  private:

  public:
    // x: current state
    // dxdt: result
    // t: current time
    //
    // state_type must match the state_type that was deduced for
    // 'y' in the call ODEINT::integrate (s. below).
    //
    void operator(const state_type& x, state_type& dxdt, double t);
  };


  and export a driver function e.g:


  RcppExport
  SEXP myode(SEXP rpars, SEXP rx, SEXP rt0, SEXP rt1, SEXP rstepper_name)
  {
    BEGIN_RCPP
      auto pars = as<List>(rpars);

      auto x = as<std::vector<double>>(rx);
      auto t0 as<double>(rt0);
      auto t1 as<double>(rt1);
      auto stepper_name = as<std::string>(Stepper);

      auto ode = std::make_unique<MyOde>(pars);
      ODEINT::integrate(stepper_name, std::move(ode), y, t0, t1);

      Rcpp::RObject rcpp_result_gen = y;
      return rcpp_result_gen;
    END_RCPP
  }


*/

// [[Rcpp::depends(BH)]]

#include <memory>
#include <exception>
#include <boost/numeric/odeint.hpp>


namespace odeint_helper {


#ifndef ODEINT_DEFAULT_STEPPER_NAME
#define ODEINT_DEFAULT_STEPPER_NAME "odeint::runge_kutta_fehlberg78"
#endif


#ifndef ODEINT_DEFAULT_ATOL
#define ODEINT_DEFAULT_ATOL 1E-6
#endif


#ifndef ODEINT_DEFAULT_RTOL
#define ODEINT_DEFAULT_RTOL 1E-6
#endif


  constexpr char default_stepper_name[] = ODEINT_DEFAULT_STEPPER_NAME;
  constexpr double default_atol = ODEINT_DEFAULT_ATOL;
  constexpr double default_rtol = ODEINT_DEFAULT_RTOL;


  template <
    typename ODE,
    typename STATE
  >
  inline void integrate(const std::string& stepper_name, ODE ode, STATE& start_state, double t0, double t1, double dt, double atol, double rtol)
  {
    namespace bno = boost::numeric::odeint;

    if ("odeint::runge_kutta_cash_karp54" == stepper_name) {
      bno::integrate_adaptive(bno::make_controlled<bno::runge_kutta_cash_karp54<STATE>>(atol, rtol), ode, start_state, t0, t1, dt);
    }
    else if ("odeint::runge_kutta_fehlberg78" == stepper_name) {
      bno::integrate_adaptive(bno::make_controlled<bno::runge_kutta_fehlberg78<STATE>>(atol, rtol), ode, start_state, t0, t1, dt);
    }
    else if ("odeint::runge_kutta_dopri5" == stepper_name) {
      bno::integrate_adaptive(bno::make_controlled<bno::runge_kutta_dopri5<STATE>>(atol, rtol), ode, start_state, t0, t1, dt);
    }
    else if ("odeint::bulirsch_stoer" == stepper_name) {
      bno::integrate_adaptive(bno::bulirsch_stoer<STATE>(atol, rtol), ode, start_state, t0, t1, dt);
    }
    else {
      throw std::runtime_error("ODEINT::integrate: unknown stepper");
    }
  }


  template <
    typename ODE,
    typename STATE
  >
  inline void integrate(const std::string& stepper_name, ODE ode, STATE& start_state, double t0, double t1, double dt)
  {
    integrate(stepper_name, ode, start_state, t0, t1, dt, default_atol, default_rtol);
  }


  template <
    typename ODE,
    typename STATE
  >
  inline void integrate(ODE ode, STATE& start_state, double t0, double t1, double dt, double atol, double rtol)
  {
    integrate(default_stepper_name, ode, start_state, t0, t1, dt, atol, rtol);
  }


  template <
    typename ODE,
    typename STATE
  >
  inline void integrate(ODE ode, STATE& start_state, double t0, double t1, double dt)
  {
    integrate(default_stepper_name, ode, start_state, t0, t1, dt);
  }

  
  // for debuging
  template <typename IT>
  inline void printf(IT first, IT last)
  {
    for (; first != last; ++first) {
      Rprintf("%f ", *first);
    }
    Rprintf("\n");
  }
  
}

#endif
