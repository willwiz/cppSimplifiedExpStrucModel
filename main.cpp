//
//  main.cpp
//  testmain
//
//  Created by Will Zhang on 11/29/17.
//  Copyright © 2017 Will Zhang. All rights reserved.
//

#include <iostream>
#include "BetaDistribution.hpp"
#include "ModelCollagenTissueExp.hpp"
#include "ModelPSstruc.hpp"

int main(int argc, const char * argv[]) {
    
//    col_ens_PF_exp model(30,1.2);
    
//    col_tissue_PF_exp model(1000, 0.0, 0.8, 30, 1.0);
    double parameters[23] = {100, 0.0, 1.0 , 1.0, 100000, 0.0, 0.3, 1.2, 0.015, 1.0, 1.24, 1.0, 1000, 1.0, 0.0, 0.0, 1.0, 0.0001, 1000, 1.05, 0.0, 0.0, 1.05};
    model_ps_struc model(parameters);
    double strain[4] = {1.2, 0.0, 0.0, 1.4};
    double res[4];
//    res = model.strainenergy(strain);

//    beta_PDF pdf(1.2, 0.02, 1.12, 1.25);
    
    beta_ODF odf(0.0,0.8);
    
    int err;
    
    std::clock_t    start;
    
    start = std::clock();
    err = model.stress(strain, res);
    std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    
    std::cout << res[0] << ",\n" << res[1] << ",\n" << res[2] <<",\n" << res[3] <<",\n";
    
    std::cout << "\n" << model.fraction_old << ",\n" << model.fraction_new << "\n";
    
//    for (int i = 0; i < 10; i++) {
//        std::cout << odf.ODF(0.1*i*M_PI) << ",\n";
//    }
//
//    std::cout << odf.is_isotropic << ",\n" << odf.has_baseline << ",\n";
//    std::cin >> x;
//    std::cout << ;
//    std::cout << pdf.PDF(x);
//    for (int i = 0; i < 20; i++) {
//        std::cout << pdf.PDF(1.0 + i*0.3/20.0) <<"\n";
//    }
//    std::cout << pdf.alpham1 <<"\n" << pdf.betam1 <<"\n" << pdf.normalizing_constant <<"\n" << pdf.betaexponent <<"\n";
    return 0;
}
