#ifndef GAUSSIANFIT_H_
#define GAUSSIANFIT_H_

#include "alglib/cpp/src/interpolation.h"
#include "Eigen/Dense"
#include <vector>

class gaussian_fitting
{
private:
public:
    gaussian_fitting();
    ~gaussian_fitting();
    
    /**
     * @brief Should only fit to data near gaussian peak, 
        as the result becomes unstable when fitting to data at gaussian tail even with low noise, 
        see article for more information, rather useless
     * @note   
     * @param  x: x axis values
     * @param  y: y axis values
     * @retval [ A, mu, sigma ]
     */
    std::vector<double> fit_caruana(std::vector<double> x,std::vector<double> y);
    
    /**
     * @brief  Standard non linear Levenberg-Marquardt least squares fit
     * @note   Intial guess u0 matters, so if your data is not centered at 0 calculate mu and use it as initial guess, 
        as finite floating point arithmetic can cause issues with gradient being 0, A and sigma can USUALLY be initialised as 1  
     * @param  x: x axis values
     * @param  y: y axis values
     * @param  u0: initial guess [ A, mu, sigma ]
     * @param  params: vector of estimated parameters [ A, mu, sigma ]
     * @param  rep: report on info like RMS error etc.
     * @param  tol: tolerance
     * @retval error code, should be err > 0
     */
    int fit_lev_marq(std::vector<double> x,std::vector<double> y,std::vector<double> u0,std::vector<double>& params,alglib::lsfitreport& rep,double tol=1e-6);
};

#endif