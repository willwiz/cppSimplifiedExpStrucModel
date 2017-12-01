//
//  EffectiveModel.hpp
//  CompleteModelCollection
//
//  Created by Will Zhang on 9/26/17.
//  Copyright © 2017 Will Zhang. All rights reserved.
//

#ifndef EffectiveModel_hpp
#define EffectiveModel_hpp

#include <stdio.h>
#include <vector>
#include <cmath>

#endif /* EffectiveModel_hpp */






class smpfl_model
{
private:
    int ex_i;
    int ex_j;
    int ex_k;
    
    int deg[3];
    
    double vec_M[2];
    double vec_S[2];
    
    double gamma_1;
    double gamma_2;
    double gamma_3;
    
    double energy;
    
public:
    void degrees(int degrees);
    
    void calc_vec_M(double phi);
    
    void calc_vec_S(double phi);
    
    void calc_lamdbda_M2();
    void calc_lamdbda_S2();
    void calc_phi();
    
    void calc_gammas();
    
    void stress(double parm[], double strain[], int n);
    
    double partial_i(double w, double gamma1, double i);
    
    smpfl_model(double phi, double rank, double degrees);
    
    smpfl_model(): ex_i{2}, ex_j{2}, ex_k{2}, deg{10,10,5}, vec_M{1,0}, vec_S{0,1} {};
};

void smpfl_model::degrees(int degrees)
{
    deg[0] = degrees;
    deg[1] = degrees;
    deg[2] = degrees/2;
    return;
}

smpfl_model::smpfl_model(double theta, int rank, int exponent)
{
    degrees(rank);
    ex_i = exponent/deg[0]/deg[1];
    ex_j = (exponent % (deg[0]/deg[1])) / deg[2];
    ex_k = (exponent % (deg[0]/deg[1])) % deg[2];
    calc_vec_M(theta);
    calc_vec_S(theta);
    return;
}



void smpfl_model::calc_vec_M(double theta)
{
    vec_M[0] = cos(theta);
    vec_M[1] = sin(theta);
    return;
}

void smpfl_model::calc_vec_S(double theta)
{
    vec_S[0] = -sin(theta);
    vec_S[1] = cos(theta);
    return;
}



double mat_multiply(double tsr[4], double vec_M[2], double vec_S[2])
{
    double result = 0.0;
    
    for (int i = 0; i<2; i++) {
        for (int j = 0; j<2; j++) {
            result += vec_M[i]*tsr[2*i + j]*vec_S[j];
        }
    }
    return result;
}




double strainenergy(double gamma1, double gamma2, double gamma3, int i, int j, int k)
{
    return pow(gamma1, i) * pow(gamma2, j) * pow(gamma3, k);
}



double smpfl_model::partial_i(double w, double gamma1, double i)
{
    return (i - 1) * w / gamma1;
}


