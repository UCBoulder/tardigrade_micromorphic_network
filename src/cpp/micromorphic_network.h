/**
  ******************************************************************************
  * \file micromorphic_network.h
  ******************************************************************************
  * The header file for the micromorphic network constitutive model.
  ******************************************************************************
  */

#ifndef VIPOR_H
#define VIPOR_H

#include<error_tools.h>
#define USE_EIGEN
#include<vector_tools.h>
#include<constitutive_tools.h>

namespace microNet{

    typedef errorTools::Node errorNode; //!< Define a error Node (node in a linked list) type
    typedef errorNode* errorOut; //!< Define the error output type
    typedef double floatType; //!< Define the float values type.
    typedef std::vector< floatType > floatVector; //!< Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!< Define a matrix of floats

    const unsigned int _spatialDimension = 3; //!< Define the spatial dimension

    errorOut computeEnergy( const floatVector &deformationGradient, const floatVector &microDeformation,
                            const floatVector &gradientMicroDeformation, const floatMatrix &params,
                            floatType &energy );

    errorOut assembleStiffness( const floatMatrix &params, floatVector &stiffnessMatrix );

    errorOut computeStress( const floatVector &deformationGradient, const floatVector &microDeformation,
                            const floatVector &gradientMicroDeformation, const floatMatrix &params,
                            floatVector &PK2Stress, floatVector &SIGMA, floatVector &higherOrderStress );

}

#endif
