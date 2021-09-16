#include "gaussian_fitting.hpp"
#include "alglib/cpp/src/optimization.h"

gaussian_fitting::gaussian_fitting()
{
}

std::vector<double> gaussian_fitting::fit_caruana(std::vector<double> x,std::vector<double> y) {
        Eigen::Matrix<double,3,3> A;
        Eigen::Vector3d b;
        A(0,0)=x.size();

        for (size_t i = 0; i < x.size(); i++)
        {
            A(0,1)+=x[i];
            A(0,2)+=pow(x[i],2);
            A(1,2)+=pow(x[i],3);
            A(2,2)+=pow(x[i],4);
            b(0)+=log(y[i]);
            b(1)+=log(y[i])*x[i];
            b(2)+=log(y[i])*pow(x[i],2);
        }
        
        A(1,0)=A(0,1);
        A(1,1)=A(0,2);
        A(2,0)=A(0,2);
        A(2,1)=A(1,2);
        Eigen::Vector3d x_=A.ldlt().solve(b);

        double mu=-0.5*x_(1)/x_(2);
        double sigma=sqrt(-0.5/x_(2));
        double a=exp(x_(0)-0.25*pow(x_(1),2)/x_(2));
        std::vector<double> params;
        params.push_back(mu);
        params.push_back(sigma);
        params.push_back(a);

        return params;
    }

void gauss_func(const alglib::real_1d_array &u, const alglib::real_1d_array &x, double &func, void *ptr) 
{
    // this callback calculates f(x,un)=A*exp(-0.5*(x-mu)^2/sigma^2)
    // where x is a position on X-axis and un is adjustable parameter
    func=u[0]*exp(-0.5*pow((x[0]-u[1])/u[2],2));
}
void gauss_grad(const alglib::real_1d_array &u, const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr) 
{
    // this callback calculates f(x,un)=A*exp(-0.5*(x-mu)^2/sigma^2) and gradient G={df/dc[i]}
    // where x is a position on X-axis and u is adjustable parameter.
    func=u[0]*exp(-0.5*pow((x[0]-u[1])/u[2],2));
    grad[0]=func/u[0];
    grad[1]=func*(x[0]-u[1])/pow(u[2],2);
    grad[2]=func*pow((x[0]-u[1]),2)/pow(u[2],3);
}

int gaussian_fitting::fit_lev_marq(std::vector<double> x,std::vector<double> y,std::vector<double> u0,std::vector<double>& params,alglib::lsfitreport& rep,int maxit,double tol) {
    alglib::real_2d_array x_;
    x_.setlength(x.size(),1);
    x_.setcontent(x.size(),1,x.data());
    alglib::real_1d_array f;
    f.setlength(y.size());
    f.setcontent(y.size(),y.data());
    alglib::real_1d_array u;
    u.setlength(u0.size());
    u.setcontent(u0.size(),u0.data());    
    alglib::ae_int_t maxits=maxit;
    alglib::ae_int_t info;
    alglib::lsfitstate state;
    // alglib::lsfitreport rep;
    lsfitcreatefg(x_, f, u, true, state);
    lsfitsetcond(state, tol, maxits);
    alglib::lsfitfit(state, gauss_func, gauss_grad);
    lsfitresults(state, info, u, rep);
    params.push_back(u[0]);
    params.push_back(u[1]);
    params.push_back(u[2]);
    
    return int(info);
}
