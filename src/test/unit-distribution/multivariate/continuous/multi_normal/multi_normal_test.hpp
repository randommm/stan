// Arguments: vvectors, vvectors, matrices
#include <stan/prob/distributions/multivariate/continuous/multi_normal.hpp>
#include <stan/math/constants.hpp>
#include <stan/math/functions/square.hpp>

using std::vector;
using std::numeric_limits;
using stan::agrad::var;

class AgradDistributionMultiNormal : public AgradDistributionTest {
public:
  void valid_values(vector<vector<vector<double> > >& parameters,
                    vector<double>& log_prob) {
    {
      vector<vector<double> > param;
      vector<double> param_y(3);
      vector<double> param_mu(3);
      vector<double> param_sigma(9);

      param_y[0] = 2.0;
      param_y[1] = -2.;
      param_y[2] = 11.0;
      param_mu[0] = 1.0;
      param_mu[1] = -1.0;
      param_mu[2] = 3.0;
      param_sigma[0] = 9.0;
      param_sigma[1] = -3.0;
      param_sigma[2] = 0.0;
      param_sigma[3] = -3.0;
      param_sigma[4] = 4.0;
      param_sigma[5] = 0.0;
      param_sigma[6] = 0.0;
      param_sigma[7] = 0.0;
      param_sigma[8] = 5.0;
      
      param.push_back(param_y);
      param.push_back(param_mu);
      param.push_back(param_sigma);
      
      parameters.push_back(param);
      log_prob.push_back(-11.73908); // expected log_prob                 
    }
                  
  }
 
  void invalid_values(vector<size_t>& index, 
          vector<double>& value) {
  }

  template <typename T_y, typename T_loc, typename T_scale,
      typename T3, typename T4, typename T5, 
      typename T6, typename T7, typename T8, 
      typename T9>
  typename stan::return_type<T_y, T_loc, T_scale>::type 
  log_prob(const T_y& y, const T_loc& mu, const T_scale& sigma,
     const T3&, const T4&, const T5&, const T6&, const T7&, const T8&, const T9&) {
    return stan::prob::multi_normal_log(y, mu, sigma);
  }

  template <bool propto, 
      typename T_y, typename T_loc, typename T_scale,
      typename T3, typename T4, typename T5, 
      typename T6, typename T7, typename T8, 
      typename T9>
  typename stan::return_type<T_y, T_loc, T_scale>::type 
  log_prob(const T_y& y, const T_loc& mu, const T_scale& sigma,
     const T3&, const T4&, const T5&, const T6&, const T7&, const T8&, const T9&) {
    return stan::prob::multi_normal_log<propto>(y, mu, sigma);
  }
  
  
  template <typename T_y, typename T_loc, typename T_scale,
      typename T3, typename T4, typename T5, 
      typename T6, typename T7, typename T8, 
      typename T9>
  var log_prob_function(const T_y& y, const T_loc& mu, const T_scale& sigma,
      const T3&, const T4&, const T5&, const T6&, const T7&, const T8&, const T9&) {
    using stan::prob::include_summand;
    using stan::math::pi;
    using stan::math::square;
    var lp(3.3);
    /*if (include_summand<true,T_y,T_loc,T_scale>::value)
      lp -= 0.5 * (y - mu) * (y - mu) / (sigma * sigma);
    if (include_summand<true,T_scale>::value)
      lp -= log(sigma);
    if (include_summand<true>::value)
      lp -= log(sqrt(2.0 * pi()));*/
    return lp;
  }
};

