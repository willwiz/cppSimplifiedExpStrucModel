//
//  ModelCollagenEnsemble.hpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/11/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#ifndef ModelCollagenEnsemble_hpp
#define ModelCollagenEnsemble_hpp

#include "math.h"
#include "Gauss_Quad_Points.hpp"
#include "BetaDistribution.hpp"

class model_col_fiber_PF
{
    
private:
    double m_1_lf; // 1/lambda_fiber, which is the applied stretch
    
public:
    model_col_fiber_PF():m_1_lf(1.0){};
    
    model_col_fiber_PF(double lf):m_1_lf(1.0/lf){
    }
    
    double stress(double ls)
    {
        double m_1_ls = 1.0/ls;
        return m_1_ls * (m_1_ls - m_1_lf);
    }
    //  Given the stretch of the fibers, this gives the stress for any slack stretch
    //      This is done to save 1 operation when integrating over the fiber model.
    
};


class model_col_ens_struc_PF
{
    //PF = Linear in 1st PK and stretch
    
    //  The parameters of this model are:
    //      mean
    //      stdev
    //      lb
    //      ub
    //
    //  The arguments are
    //      lf
    //      l_pushforward
    
    
private:
    const int gq_n = 21;    // The number of gauss quadrature points useed
    const double *gq_abs = C_01_A21;    // The non transformed abscissas
    const double *gq_weight = C_01_W21; // The non-transformed weights
    
public:
    
    int is_deltafunction;
    //  Check if the recruitment function is too narrow for integration
    //      In this case, the recruitment distribution is assumed to behaves as a delta function and proceed
    //      accordingly
    
    double lower_limit;
    double upper_limit;
    
    double mean;
    double stdev;
    
    model_col_fiber_PF fibermodel;
    //
    
    beta_PDF Dx_PDF;
    //The collagen fiber recruitment distribution
    //  Which has the parameters:
    //      mean,   mean
    //      stdev,  standard deviation
    //      lb,     lowerbound
    //      ub,     upperbound
    
    void set_parameters(double mu, double sigma, double lb, double ub);
    
    model_col_ens_struc_PF(): mean(1.2), stdev(0.015), lower_limit(1.0), upper_limit(1.24), Dx_PDF(1.2, 0.12, 1.12, 1.24), is_deltafunction(0)
    {
    };
    // default constructor with common parameters for valve tissues
    
    model_col_ens_struc_PF(double mu, double sigma, double lb, double ub):
    mean(mu), stdev(sigma), lower_limit(lb), upper_limit(ub), Dx_PDF(mu, sigma, lb, ub), is_deltafunction(0)
    {
    }
    
    
    
    double stress(double lambda_fiber, double lambda_pushforward);
};

#endif /* ModelCollagenEnsemble_hpp */
