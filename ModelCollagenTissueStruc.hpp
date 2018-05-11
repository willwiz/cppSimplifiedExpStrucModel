//
//  ModelCollagenTissueStruc.hpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/11/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#ifndef ModelCollagenTissueStruc_hpp
#define ModelCollagenTissueStruc_hpp

#include "kinematics.hpp"
#include "BetaDistribution.hpp"
#include "ModelCollagenEnsemble.hpp"

void calc_C(double vT[4], double res[4]);


class model_col_tissue_struc_PF
{
    //  The parameters of this model are
    //      kC = arg[0]
    //      ODFmean = arg[]
    //      ODFstdev = arg[]
    //      Dxmean = arg[]
    //      Dxstdev = arg[]
    //      Dxlb = arg[]
    //      DXub = arg[]
    //      Dxaniso = arg[]
    //      kI = arg[]
    //
    //  The constants of this model is
    //      PF11 = arg[]
    //      PF12 = arg[]
    //      PF21 = arg[]
    //      PF22 = arg[]
    //
    //  The argument of this model are
    //      F11, F12, F21, F22
    
private:
    double kC;
    double ODFmean;
    double ODFstdev;
    double Dxmean;
    double Dxstdev;
    double Dxlb;
    double Dxub;
    double kI;
    double Dx_aniso;
    
    const int gq_n = 31;    // The number of gauss quadrature points useed
    const double *gq_abs = C_11_A31;    // The non transformed abscissas
    const double *gq_weight = C_11_W31; // The non-transformed weights
    
    double cos2_theta[31];
    double sin2_theta[31];
    double cossin_theta[31];
    double odfweight[31];
    
    double int_range_theta; //range of integration over the splay
    
    void compute_range(double stdev){int_range_theta = fmin(M_PI_2, 4.0 * stdev);}
    
public:
    
    beta_ODF fiber_odf;
    model_col_ens_struc_PF ensemblemodel;
    
    void set_parameters();
    
    model_col_tissue_struc_PF(): kC{100000.0}, kI{2000.0}, fiber_odf(0.0, 0.3),
        ensemblemodel(1.2, 0.015, 1.0, 1.24), Dx_aniso{1.0}, int_range_theta(M_PI_2)
    {
        double theta; //    orientation in the current reference state
        double omega; //    orientation in the original reference state
        double dummie_cos;
        double dummie_sin;
        for (int i = 0; i < 31; i++) {
            theta = int_range_theta * C_11_A31[i];
            dummie_cos = cos(theta);
            dummie_sin = sin(theta);
            omega = theta;
            cos2_theta[i] = dummie_cos*dummie_cos;
            sin2_theta[i] = dummie_sin*dummie_sin;
            cossin_theta[i] = dummie_cos*dummie_sin;
            odfweight[i] = int_range_theta * C_11_W31[i] * fiber_odf.ODF(theta);
        }
    };
    
    model_col_tissue_struc_PF(double arg[13]):
    kC{arg[0]},
    ODFmean{arg[1]},
    ODFstdev{arg[2]},
    kI{arg[8]},
    Dx_aniso{arg[7]},
    fiber_odf(arg[1], arg[2]),
    ensemblemodel(arg[3], arg[4], arg[5], arg[6]),
    int_range_theta(M_PI_2)
    {
        compute_range(arg[2]);
        double theta;
        double hereonly_cos;
        double hereonly_sin;
        for (int i = 0; i < 31; i++) {
            theta = int_range_theta * C_11_A31[i] + ODFmean;
            hereonly_cos = cos(theta);
            hereonly_sin = sin(theta);
            cos2_theta[i] = hereonly_cos*hereonly_cos;
            sin2_theta[i] = hereonly_sin*hereonly_sin;
            cossin_theta[i] = hereonly_cos*hereonly_sin;
            odfweight[i] = int_range_theta * C_11_W31[i] * fiber_odf.ODF(theta);
        }
    };
    
    int stress(double vC[4], double res[4]);
    
    
    
};


#endif /* ModelCollagenTissueStruc_hpp */
