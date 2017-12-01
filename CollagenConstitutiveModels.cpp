//
//  CollagenConstitutiveModels.cpp
//  CompleteModelCollectionTest
//
//  Created by Will Zhang on 9/1/17.
//  Copyright © 2017 Will Zhang. All rights reserved.
//

#include "CollagenConstitutiveModels.hpp"

/*****          Collagen fiber models           *****/

//Fiber model that's linear in PF
double col_fiber_PF (double Ef, double Es){
    double lf = 1.0 / Ef;
    double ls = 1.0 / Es;
    
    return ls * (ls - lf);
};




//Fiber Model that's linear in SE
double col_fiber_SE (double lf, double ls)
{
    double Ef = lf*lf;
    double Es = ls*ls;
    double denom = 1.0/Es;
    
    return 0.5*(Ef - Es) * denom * denom;
};




/*****          Exponential type ensemble models           *****/

double col_ens_PF_exp::stress(double lf)
{
    return (exp(parB * (lf - lmax)) - scaling);
}

double col_ens_PF_exp::strainenergy(double lf)
{
    double lfm1 = lf - 1.0;
    return scaling*( (exp(parB * lfm1) - 1.0)/parB - lfm1);
}


/*****          Exponential type tissue models           *****/



