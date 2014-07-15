// Arguments: Ints, Doubles, Doubles, Doubles
#include <stan/prob/distributions/univariate/discrete/zero_inflated.hpp>
#include <stan/prob/distributions/univariate/discrete/binomial.hpp>

#include <stan/math/functions/multiply_log.hpp>
#include <stan/math/functions/log1m.hpp>
#include <stan/math/functions/binomial_coefficient_log.hpp>

using std::vector;
using std::numeric_limits;
using stan::agrad::var;

class AgradDistributionsZINegBinomialLog : public AgradDistributionTest {
public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(4);

    /* Density in R:
    require(VGAM); eta <- 2.0; phi <- 1; zeta <- -1.0; dzinegbin(0,
    munb=exp(eta), size=phi, pstr0=1/(1+exp(-zeta)), log=TRUE)
    */
    param[0] = 0;           // n
    param[1] = 2.0;          // eta
    param[2] = 1.0;          // phi
    param[3] = -1.0;          // zeta
    parameters.push_back(param);
    log_prob.push_back(-1.032584); // expected log_prob

    param[0] = 1;           // n
    param[1] = 2.0;          // eta
    param[2] = 1.0;          // phi
    param[3] = -1.0;          // zeta
    parameters.push_back(param);
    log_prob.push_back(-2.567118); // expected log_prob
    
  }
 
  void invalid_values(vector<size_t>& index, 
                      vector<double>& value) {
    // n
    index.push_back(0U);
    value.push_back(-1);
    
    // eta
    //index.push_back(1U);
    //value.push_back(0);
    
    // zeta
    //index.push_back(2U);
    //value.push_back(0);
  }

  template <class T_n, class T_log_location, class T_inv_scale,
            class T_zi, typename T4, typename T5, 
            typename T6, typename T7, typename T8, 
            typename T9>
  typename stan::return_type<T_log_location, T_inv_scale, T_zi>::type 
  log_prob(const T_n& n, const T_log_location& eta,
     const T_inv_scale& phi, const T_zi& zeta, const T4&, const T5&, 
     const T6&, const T7&, const T8&, 
     const T9&) {
    return stan::prob::zi_neg_binomial_log_log(n, eta, phi, zeta);
  }

  template <bool propto, 
            class T_n, class T_log_location, class T_inv_scale,
            class T_zi, typename T4, typename T5, 
            typename T6, typename T7, typename T8, 
            typename T9>
  typename stan::return_type<T_log_location, T_inv_scale, T_zi>::type 
  log_prob(const T_n& n, const T_log_location& eta,
     const T_inv_scale& phi, const T_zi& zeta, const T4&, const T5&, 
     const T6&, const T7&, const T8&, 
     const T9&) {
    return stan::prob::zi_neg_binomial_log_log<propto>(n, eta, phi, zeta);
  }

  template <class T_n, class T_log_location, class T_inv_scale,
            class T_zi, typename T4, typename T5, 
            typename T6, typename T7, typename T8, 
            typename T9>
  var log_prob_function(const T_n& n, const T_log_location& eta,
      const T_inv_scale& phi, const T_zi& zeta, const T4&, const T5&, 
      const T6&, const T7&, const T8&, 
      const T9&) {
    using stan::math::binomial_coefficient_log;
    using stan::math::log_sum_exp;
    using stan::math::log1p_exp;
    using boost::math::lgamma;
    using stan::prob::include_summand;
    using stan::math::multiply_log;
    using std::exp;
    using std::log;

    var logp(0);
    if (n==0) {
      if (include_summand<true, T_log_location,
                          T_inv_scale, T_zi>::value)
        logp += log(exp(zeta) + pow((phi/(exp(eta) + phi)), phi));
      if (include_summand<true, T_zi>::value)
        logp -= log1p_exp(zeta);
    } else {
      if (include_summand<true, T_inv_scale>::value)
        logp += multiply_log(phi, phi) - lgamma(phi);
      if (include_summand<true, T_log_location,T_inv_scale>::value)
        logp -= (n+phi)*log_sum_exp(eta, log(phi));
      if (include_summand<true, T_log_location>::value)
        logp += n*eta;
      if (include_summand<true, T_inv_scale>::value)
        logp += lgamma(n+phi);   
      if (include_summand<true, T_zi>::value)
        logp -= log1p_exp(zeta);
    }
    return logp;
  }
};
