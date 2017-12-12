//
//  MathLibLink.cpp
//  EXL_Model
//
//  Created by Will Zhang on 7/22/15.
//  Copyright (c) 2015 Will Zhang. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include "WolframLibrary.h"
#include "BetaDistribution.hpp"
#include "CollagenConstitutiveModels.hpp"

/* Return the version of Library Link */

EXTERN_C DLLEXPORT mint WolframLibrary_getVersion( ) {
    return WolframLibraryVersion;
}

/* Initialize Library */
EXTERN_C DLLEXPORT int WolframLibrary_initialize( WolframLibraryData libData) {
    return LIBRARY_NO_ERROR;
}

/* Uninitialize Library */
EXTERN_C DLLEXPORT void WolframLibrary_uninitialize( WolframLibraryData libData) {
    return;
}





/* */
EXTERN_C DLLEXPORT int BetaODF( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    
    double result;
    //Start of Mathematica to C Communication
    
    //get the inputs
    mreal T1 = MArgument_getReal(Args[0]);
    mreal T2 = MArgument_getReal(Args[1]);
    mreal T3 = MArgument_getReal(Args[2]);
    
    beta_ODF fiberODF(T1, T2);
    
    result = fiberODF.ODF(T3);
    
    MArgument_setReal(Res, result); //Set results to be the actually output
    
    return LIBRARY_NO_ERROR;
}





/* */
EXTERN_C DLLEXPORT int BetaPDF( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    
    //Start of Mathematica to C Communication
    
    //get the inputs
    mreal mu = MArgument_getReal(Args[0]);
    mreal sd = MArgument_getReal(Args[1]);
    mreal lb = MArgument_getReal(Args[2]);
    mreal ub = MArgument_getReal(Args[3]);
    mreal x = MArgument_getReal(Args[4]);
    

    
    beta_PDF fiberPDF(mu, sd, lb, ub);
    
    double result = fiberPDF.PDF(x);
    
    MArgument_setReal(Res, result); //Set results to be the actually output
    return LIBRARY_NO_ERROR;
}




EXTERN_C DLLEXPORT int ColExpModel( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    
    //Start of Mathematica to C Communication
    
    int err;
    
    MTensor T0;
    MTensor T1;
    MTensor Tout;
    mint dims[1];
    
    double* vC;
    double* vparm;
    double* results;
    
    dims[0] = 4;
    
    T0 = MArgument_getMTensor(Args[0]);
    T1 = MArgument_getMTensor(Args[1]);
    
    vparm = libData->MTensor_getRealData(T0);
    vC = libData->MTensor_getRealData(T1);
    
    double kc = vparm[0];
    double mu = vparm[1];
    double sd = vparm[2];
    double expB = vparm[3];
    double lfmax = vparm[4];
    
    err = libData->MTensor_new(MType_Real, 1, dims, &Tout);
    
    results = libData->MTensor_getRealData(Tout);
    
    col_tissue_PF_exp model(kc, mu, sd, expB, lfmax);
    
    err = model.stress(vC, results); //Compute results
    
    MArgument_setMTensor(Res, Tout); //Set results to be the actually output
    
    return LIBRARY_NO_ERROR;
}



EXTERN_C DLLEXPORT int ColExpModelPsi( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    
    //Start of Mathematica to C Communication
    
    //Mathematica type for the inputs
    MTensor T0;
    MTensor T1;
    
    //The inputs from mathematica
    double* vC;
    double* vparm;
    
    //The type for the result
    double results;
    
    //Take input from mathematica
    T0 = MArgument_getMTensor(Args[0]);
    T1 = MArgument_getMTensor(Args[1]);
    
    //Convert it to C types
    vparm = libData->MTensor_getRealData(T0);
    vC = libData->MTensor_getRealData(T1);
    
    //Define the parameters from the input vector
    double kc = vparm[0];
    double mu = vparm[1];
    double sd = vparm[2];
    double expB = vparm[3];
    double lfmax = vparm[4];
    
    //Create the model object
    col_tissue_PF_exp model(kc, mu, sd, expB, lfmax);
    
    //Compute the strain energy and set it to the results
    results = model.strainenergy(vC); //Compute results
    
    //Send the results to mathematica
    MArgument_setReal(Res, results); //Set results to be the actually output
    
    return LIBRARY_NO_ERROR;
}
