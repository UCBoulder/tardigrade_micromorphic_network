/**
  ******************************************************************************
  * \file micromorphic_network.h
  ******************************************************************************
  * The header file for the micromorphic network constitutive model.
  ******************************************************************************
  */

#ifndef VIPOR_H
#define VIPOR_H

#define USE_EIGEN
#include<vector_tools.h>

namespace micronet{

    errorOut compute_energy( const floatVector &deformationGradient, const floatVector &microDeformation,
                             const floatVector &gradientMicroDeformation, const floatMatrix &params,
                             floatType &energy );

}

#endif
