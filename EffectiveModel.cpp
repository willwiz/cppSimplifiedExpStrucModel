//
//  EffectiveModel.cpp
//  CompleteModelCollection
//
//  Created by Will Zhang on 9/26/17.
//  Copyright © 2017 Will Zhang. All rights reserved.
//

#include "EffectiveModel.hpp"

void smpfl_model::stress(double *parm, double *strain, int n)
{
//    double gamma1 = strain[0];
//    double gamma2 = strain[1];
//    double gamma3 = strain[2];
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += strain[i];
    }
    return;
}

void smpfl_model::degrees(int deg, int num)
{
    return;
}
