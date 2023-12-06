/**
  ******************************************************************************
  * \file micromorphic_network.cpp
  ******************************************************************************
  * The header file for the micromorphic network constitutive model.
  ******************************************************************************
  */

#include<micromorphic_network.h>

namespace microNet{

    errorOut computeEnergy( const floatVector &deformationGradient, const floatVector &microDeformation,
                            const floatVector &gradientMicroDeformation, const floatMatrix &params,
                            floatType &energy ){
        /*!
         * Compute the energy of the micromorphic network.
         * 
         * \f$\left(\psi \rho\right) = \frac{1}{2} E_{IJ}^{\chi} C_{IJKL} E_{KL}^{\chi} + \sum_{n=1}^N \frac{1}{4} k^n \left[ \sqrt{c_i c_i} - \sqrt{C_I C_I} \right]^2\f$
         * 
         * Where $C$ is the stiffness matrix of the central mass, and the vectors which describe the
         * spring in the current and reference configurations are
         * 
         * \f$c_i = F_{iI} dX_I + \left( \chi_{iI} + \chi_{iI,J} dX_J \right) \Xi_I^2 - \chi_{iI} \Xi_I^1 \f$
         * 
         * \f$C_I = dX_I + \Xi_I^2 - \Xi_I^1 \f$
         * 
         * \f$k^n\f$ is the stiffness of the connecting linear spring, \f$dX_I\f$ are the components of the vector
         * describing the center to center spacing of the network, \f$\Xi_I^1\f$ and \f$\Xi_I^2\f$ are the endpoints
         * of the spring on the local (1) and neighboring (2) network object.
         * 
         * \param &deformationGradient: The total deformation gradient at the current increment.
         * \param &microDeformation: The total micro-deformation at the current increment.
         * \param &gradientMicroDeformation: The gradient of the micro-deformation w.r.t. the reference
         *     configuration.
         * \param &params: The parameter matrix. The first index are the parameters for the internal
         *     stiffness matrix and the remaining indices describe the springs.
         * \param &energy: The computed strain energy.
         */

        // Initialize the energy
        energy = 0;

        // Compute the inverse of the deformation gradient
        floatVector inverseDeformationGradient = vectorTools::inverse( deformationGradient,
                                                                       _spatialDimension, _spatialDimension );

        // Compute the stiffness matrix of the central portion
        floatVector C_center;
        errorOut error = assembleStiffness( params, C_center );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the assembly of the central portion's stiffness matrix." );
            result->addNext( error );
            return result;

        }

        floatVector microStrain;
        error = constitutiveTools::computeGreenLagrangeStrain( microDeformation, microStrain );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the micro Green-Lagrange strain" );
            result->addNext( error );
            return result;

        }

        // Add the energy of the central section
        for ( unsigned int I = 0; I < _spatialDimension; I++ ){

            for ( unsigned int J = 0; J < _spatialDimension; J++ ){

                for ( unsigned int K = 0; K < _spatialDimension; K++ ){

                    for ( unsigned int L = 0; L < _spatialDimension; L++ ){

                        energy += 0.5 * microStrain[ 3 * I + J ] * C_center[ 27 * I + 9 * J + 3 * K + L ] * microStrain[ 3 * K + L ];

                    }

                }

            }

        }

        // Compute the energy from the springs
        for ( auto p = params.begin( ) + 1; p != params.end( ); p++ ){

            if ( p->size( ) != 10 ){

                std::ostringstream message;
                message << "The parameters for the " << ( int )( p - params.begin( ) ) + 1
                        << "th spring has a size of " << p->size( ) << " and 10 values are required.";

                return new errorNode( __func__, message.str( ) );

            }

            // Extract the parameters
            floatType kn     = ( *p )[ 0 ];
            floatVector dX   = { ( *p )[ 1 ], ( *p )[ 2 ], ( *p )[ 3 ] };
            floatVector Xi_1 = { ( *p )[ 4 ], ( *p )[ 5 ], ( *p )[ 6 ] };
            floatVector Xi_2 = { ( *p )[ 7 ], ( *p )[ 8 ], ( *p )[ 9 ] };

            floatVector dx( _spatialDimension, 0 );
            floatVector xi_1( _spatialDimension, 0 );
            floatVector xi_2( _spatialDimension, 0 );

            for ( unsigned int i = 0; i < _spatialDimension; i++ ){

                for ( unsigned int I = 0; I < _spatialDimension; I++ ){

                    dx[ i ]   += deformationGradient[ 3 * i + I ] * dX[ I ];
                    xi_1[ i ] += microDeformation[ 3 * i + I ] * Xi_1[ I ];
                    xi_2[ i ] += microDeformation[ 3 * i + I ] * Xi_2[ I ];

                    for ( unsigned int J = 0; J < _spatialDimension; J++ ){

                        xi_2[ i ] += gradientMicroDeformation[ 9 * i + 3 * I + J ] * Xi_2[ I ] * dX[ J ];

                    }

                }

            }

            // Assemble the linear spring lengths
            floatVector C = dX + Xi_2 - Xi_1;
            floatVector c = dx + xi_2 - xi_1;

            // Compute the energy contribution
            energy += 0.25 * kn * std::pow(std::sqrt( vectorTools::dot( c, c ) ) - std::sqrt( vectorTools::dot( C, C ) ), 2 );

        }

        return NULL;

    }

    errorOut assembleStiffness( const floatMatrix &params, floatVector &stiffnessMatrix ){
        /*!
         * Assemble the stiffness matrix of the central portion. Current assumption is isotropy.
         * 
         * \f$ C_{IJKL} = \lambda \delta_{IJ} \delta_{KL} + \mu \left( \delta_{IK} \delta_{JL} + \delta_{IL} \delta_{JK} \right) \f$
         * 
         * \param &params: The parameter value matrix.
         * \param &stiffnessMatrix: The resulting stiffness matrix of the central section.
         */

        if ( params.size( ) < 1 ){

            return new errorNode( __func__, "The parameter matrix must at least define the central section properties." );

        }

        if ( params[ 0 ].size( ) != 2 ){

            return new errorNode( __func__, "The central section must be a linear elastic material" );

        }

        floatType lambda = params[ 0 ][ 0 ];
        floatType mu     = params[ 0 ][ 1 ];
        floatVector eye( _spatialDimension * _spatialDimension, 0 );
        vectorTools::eye( eye );

        stiffnessMatrix = floatVector( _spatialDimension * _spatialDimension * _spatialDimension * _spatialDimension, 0 );

        const unsigned int _sd3 = _spatialDimension * _spatialDimension * _spatialDimension;
        const unsigned int _sd2 = _spatialDimension * _spatialDimension;
        const unsigned int _sd1 = _spatialDimension;

        for ( unsigned int I = 0; I < _spatialDimension; I++ ){

            for ( unsigned int J = 0; J < _spatialDimension; J++ ){

                for ( unsigned int K = 0; K < _spatialDimension; K++ ){

                    for ( unsigned int L = 0; L < _spatialDimension; L++ ){

                        stiffnessMatrix[ _sd3 * I + _sd2 * J + _sd1 * K + L ] = lambda * eye[ _sd1 * I + J ] * eye[ _sd1 * K + L ]
                                                                              + mu * ( eye[ _sd1 * I + K ] * eye[ _sd1 * J + L ] + eye[ _sd1 * I + L ] * eye[ _sd1 * J + K ] );

                    }

                }

            }

        }

        return NULL;

    }

    errorOut computeStress( const floatVector &deformationGradient, const floatVector &microDeformation,
                            const floatVector &gradientMicroDeformation, const floatMatrix &params,
                            floatVector &PK2, floatVector &SIGMA, floatVector &M ){
        /*!
         * Compute the stresses of the micromorphic network with a strain energy function of
         * 
         * \f$\left(\psi \rho\right) = \frac{1}{2} E_{IJ}^{\chi} C_{IJKL} E_{KL}^{\chi} + \sum_{n=1}^N \frac{1}{4} k^n \left[ \sqrt{c_i c_i} - \sqrt{C_I C_I} \right]^2\f$
         * 
         * Where $C$ is the stiffness matrix of the central mass, and the vectors which describe the
         * spring in the current and reference configurations are
         * 
         * \f$c_i = F_{iI} dX_I + \left( \chi_{iI} + \chi_{iI,J} dX_J \right) \Xi_I^2 - \chi_{iI} \Xi_I^1 \f$
         * 
         * \f$C_I = dX_I + \Xi_I^2 - \Xi_I^1 \f$
         * 
         * \f$k^n\f$ is the stiffness of the connecting linear spring, \f$dX_I\f$ are the components of the vector
         * describing the center to center spacing of the network, \f$\Xi_I^1\f$ and \f$\Xi_I^2\f$ are the endpoints
         * of the spring on the local (1) and neighboring (2) network object. This strain energy function
         * leads to stresses in the reference configuration of
         * 
         * \f$ S_{AB} = \sum_{n}^N \left\{f^n c_a dX_A F_{Ba}^{-1}\right\}\f$
         * 
         * \f$ \Sigma_{AB} = C_{ABIJ} E_{IJ}^{\chi} + \sum_{n}^N \left\{ f^n F_{Aa}^{-1} c_a F_{Bb}^{-1} c_b \right\}\f$
         * 
         * \f$ M_{ABC} = \sum_{n}^N \left\{ f^n c_a dX_A F_{Ba}^{-1} \Xi_C^2 \right\} \f$
         * 
         * where \f$\bf{S}\f$ is the second Piola-Kirchhoff stress, \f$\bf{\Sigma}\f$ is the reference symmetric micro stress,
         * and \f$\bf{M}\f$ is the reference higher order stress.
         * 
         * \param &deformationGradient: The total deformation gradient at the current increment.
         * \param &microDeformation: The total micro-deformation at the current increment.
         * \param &gradientMicroDeformation: The gradient of the micro-deformation w.r.t. the reference
         *     configuration.
         * \param &params: The parameter matrix. The first index are the parameters for the internal
         *     stiffness matrix and the remaining indices describe the springs.
         * \param &PK2: The second Piola-Kirchhoff stress
         * \param &SIGMA: The reference symmetric micro-stress
         * \param &M: The reference higher order stress.
         */

        floatMatrix dPK2dF,   dPK2dchi,   dPK2dgrad_chi;
        floatMatrix dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi;
        floatMatrix dMdF,     dMdchi,     dMdgrad_chi;

        return computeStress( deformationGradient, microDeformation, gradientMicroDeformation, params,
                              PK2,       SIGMA,      M,
                              dPK2dF,    dPK2dchi,   dPK2dgrad_chi,
                              dSIGMAdF,  dSIGMAdchi, dSIGMAdgrad_chi,
                              dMdF,      dMdchi,     dMdgrad_chi,
                              false );

    }

    errorOut computeStress( const floatVector &deformationGradient, const floatVector &microDeformation,
                            const floatVector &gradientMicroDeformation, const floatMatrix &params,
                            floatVector &PK2,      floatVector &SIGMA,      floatVector &M,
                            floatMatrix &dPK2dF,   floatMatrix &dPK2dchi,   floatMatrix &dPK2dgrad_chi,
                            floatMatrix &dSIGMAdF, floatMatrix &dSIGMAdchi, floatMatrix &dSIGMAdgrad_chi,
                            floatMatrix &dMdF,     floatMatrix &dMdchi,     floatMatrix &dMdgrad_chi ){
        /*!
         * Compute the stresses of the micromorphic network with a strain energy function of
         * 
         * \f$\left(\psi \rho\right) = \frac{1}{2} E_{IJ}^{\chi} C_{IJKL} E_{KL}^{\chi} + \sum_{n=1}^N \frac{1}{4} k^n \left[ \sqrt{c_i c_i} - \sqrt{C_I C_I} \right]^2\f$
         * 
         * Where $C$ is the stiffness matrix of the central mass, and the vectors which describe the
         * spring in the current and reference configurations are
         * 
         * \f$c_i = F_{iI} dX_I + \left( \chi_{iI} + \chi_{iI,J} dX_J \right) \Xi_I^2 - \chi_{iI} \Xi_I^1 \f$
         * 
         * \f$C_I = dX_I + \Xi_I^2 - \Xi_I^1 \f$
         * 
         * \f$k^n\f$ is the stiffness of the connecting linear spring, \f$dX_I\f$ are the components of the vector
         * describing the center to center spacing of the network, \f$\Xi_I^1\f$ and \f$\Xi_I^2\f$ are the endpoints
         * of the spring on the local (1) and neighboring (2) network object. This strain energy function
         * leads to stresses in the reference configuration of
         * 
         * \f$ S_{AB} = \sum_{n}^N \left\{f^n c_a dX_A F_{Ba}^{-1}\right\}\f$
         * 
         * \f$ \Sigma_{AB} = C_{ABIJ} E_{IJ}^{\chi} + \sum_{n}^N \left\{ f^n F_{Aa}^{-1} c_a F_{Bb}^{-1} c_b \right\}\f$
         * 
         * \f$ M_{ABC} = \sum_{n}^N \left\{ f^n c_a dX_A F_{Ba}^{-1} \Xi_C^2 \right\} \f$
         * 
         * where \f$\bf{S}\f$ is the second Piola-Kirchhoff stress, \f$\bf{\Sigma}\f$ is the reference symmetric micro stress,
         * and \f$\bf{M}\f$ is the reference higher order stress.
         * 
         * \param &deformationGradient: The total deformation gradient at the current increment.
         * \param &microDeformation: The total micro-deformation at the current increment.
         * \param &gradientMicroDeformation: The gradient of the micro-deformation w.r.t. the reference
         *     configuration.
         * \param &params: The parameter matrix. The first index are the parameters for the internal
         *     stiffness matrix and the remaining indices describe the springs.
         * \param &PK2: The second Piola-Kirchhoff stress
         * \param &SIGMA: The reference symmetric micro-stress
         * \param &M: The reference higher order stress.
         * \param &dPK2dF: The Jacobian of the PK2 stress w.r.t. the deformation gradient.
         * \param &dPK2dchi: The Jacobian of the PK2 stress w.r.t. the micro deformation.
         * \param &dPK2dgrad_chi: The Jacobian of the PK2 stress w.r.t. the gradient of the micro deformation.
         * \param &dSIGMAdF: The Jacobian of the reference micro stress w.r.t. the deformation gradient.
         * \param &dSIGMAdchi: The Jacobian of the reference micro stress w.r.t. the micro deformation.
         * \param &dSIGMAdgrad_chi: The Jacobian of the reference micro stress w.r.t. the gradient of the micro deformation.
         * \param &dMdF: The Jacobian of the higher order stress w.r.t. the deformation gradient.
         * \param &dMdchi: The Jacobian of the higher order stress w.r.t. the micro deformation.
         * \param &dMdgrad_chi: The Jacobian of the higher order stress w.r.t. the gradient of the micro deformation.
         */

        return computeStress( deformationGradient, microDeformation, gradientMicroDeformation, params,
                              PK2,       SIGMA,      M,
                              dPK2dF,    dPK2dchi,   dPK2dgrad_chi,
                              dSIGMAdF,  dSIGMAdchi, dSIGMAdgrad_chi,
                              dMdF,      dMdchi,     dMdgrad_chi,
                              true );

    }

    errorOut computeStress( const floatVector &deformationGradient, const floatVector &microDeformation,
                            const floatVector &gradientMicroDeformation, const floatMatrix &params,
                            floatVector &PK2,      floatVector &SIGMA, floatVector &M,
                            floatMatrix &dPK2dF,   floatMatrix &dPK2dchi,   floatMatrix &dPK2dgrad_chi,
                            floatMatrix &dSIGMAdF, floatMatrix &dSIGMAdchi, floatMatrix &dSIGMAdgrad_chi,
                            floatMatrix &dMdF,     floatMatrix &dMdchi,     floatMatrix &dMdgrad_chi,
                            bool computeJacobians ){
        /*!
         * Compute the stresses of the micromorphic network with a strain energy function of
         * 
         * \f$\left(\psi \rho\right) = \frac{1}{2} E_{IJ}^{\chi} C_{IJKL} E_{KL}^{\chi} + \sum_{n=1}^N \frac{1}{4} k^n \left[ \sqrt{c_i c_i} - \sqrt{C_I C_I} \right]^2\f$
         * 
         * Where $C$ is the stiffness matrix of the central mass, and the vectors which describe the
         * spring in the current and reference configurations are
         * 
         * \f$c_i = F_{iI} dX_I + \left( \chi_{iI} + \chi_{iI,J} dX_J \right) \Xi_I^2 - \chi_{iI} \Xi_I^1 \f$
         * 
         * \f$C_I = dX_I + \Xi_I^2 - \Xi_I^1 \f$
         * 
         * \f$k^n\f$ is the stiffness of the connecting linear spring, \f$dX_I\f$ are the components of the vector
         * describing the center to center spacing of the network, \f$\Xi_I^1\f$ and \f$\Xi_I^2\f$ are the endpoints
         * of the spring on the local (1) and neighboring (2) network object. This strain energy function
         * leads to stresses in the reference configuration of
         * 
         * \f$ S_{AB} = \sum_{n}^N \left\{f^n c_a dX_A F_{Ba}^{-1}\right\}\f$
         * 
         * \f$ \Sigma_{AB} = C_{ABIJ} E_{IJ}^{\chi} + \sum_{n}^N \left\{ f^n F_{Aa}^{-1} c_a F_{Bb}^{-1} c_b \right\}\f$
         * 
         * \f$ M_{ABC} = \sum_{n}^N \left\{ f^n c_a dX_A F_{Ba}^{-1} \Xi_C^2 \right\} \f$
         * 
         * where \f$\bf{S}\f$ is the second Piola-Kirchhoff stress, \f$\bf{\Sigma}\f$ is the reference symmetric micro stress,
         * and \f$\bf{M}\f$ is the reference higher order stress.
         * 
         * \param &deformationGradient: The total deformation gradient at the current increment.
         * \param &microDeformation: The total micro-deformation at the current increment.
         * \param &gradientMicroDeformation: The gradient of the micro-deformation w.r.t. the reference
         *     configuration.
         * \param &params: The parameter matrix. The first index are the parameters for the internal
         *     stiffness matrix and the remaining indices describe the springs.
         * \param &PK2: The second Piola-Kirchhoff stress
         * \param &SIGMA: The reference symmetric micro-stress
         * \param &M: The reference higher order stress.
         * \param &dPK2dF: The Jacobian of the PK2 stress w.r.t. the deformation gradient.
         * \param &dPK2dchi: The Jacobian of the PK2 stress w.r.t. the micro deformation.
         * \param &dPK2dgrad_chi: The Jacobian of the PK2 stress w.r.t. the gradient of the micro deformation.
         * \param &dSIGMAdF: The Jacobian of the reference micro stress w.r.t. the deformation gradient.
         * \param &dSIGMAdchi: The Jacobian of the reference micro stress w.r.t. the micro deformation.
         * \param &dSIGMAdgrad_chi: The Jacobian of the reference micro stress w.r.t. the gradient of the micro deformation.
         * \param &dMdF: The Jacobian of the higher order stress w.r.t. the deformation gradient.
         * \param &dMdchi: The Jacobian of the higher order stress w.r.t. the micro deformation.
         * \param &dMdgrad_chi: The Jacobian of the higher order stress w.r.t. the gradient of the micro deformation.
         * \param &computeJacobians: Flag for whether the Jacobians should be computed or not.
         */

        // Compute the inverse of the deformation gradient
        floatVector inverseDeformationGradient = vectorTools::inverse( deformationGradient,
                                                                       _spatialDimension, _spatialDimension );

        // Compute the stiffness matrix of the central portion
        floatVector C_center;
        errorOut error = assembleStiffness( params, C_center );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the assembly of the central portion's stiffness matrix." );
            result->addNext( error );
            return result;

        }

        floatVector microStrain;
        floatMatrix dMicroStraindchi;

        if ( computeJacobians ){
            error = constitutiveTools::computeGreenLagrangeStrain( microDeformation, microStrain, dMicroStraindchi );
        }
        else{
            error = constitutiveTools::computeGreenLagrangeStrain( microDeformation, microStrain );
        }

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in the computation of the micro Green-Lagrange strain" );
            result->addNext( error );
            return result;

        }
       
        floatVector eye( _spatialDimension * _spatialDimension, 0 );
        vectorTools::eye( eye );

        const unsigned int _sd3 = _spatialDimension * _spatialDimension * _spatialDimension;
        const unsigned int _sd2 = _spatialDimension * _spatialDimension;
        const unsigned int _sd1 = _spatialDimension;

        // Initialize the stresses
        PK2   = floatVector( _spatialDimension * _spatialDimension, 0 );
        SIGMA = floatVector( _spatialDimension * _spatialDimension, 0 );
        M     = floatVector( _spatialDimension * _spatialDimension * _spatialDimension, 0 );

        floatMatrix dFinvdF;

        if ( computeJacobians ){

            // Initialize the Jacobians

            dFinvdF = floatMatrix( _sd2, floatVector( _sd2, 0 ) );

            for ( unsigned int A = 0; A < _spatialDimension; A++ ){

                for ( unsigned int a = 0; a < _spatialDimension; a++ ){

                    for ( unsigned int K = 0; K < _spatialDimension; K++ ){

                        for ( unsigned int k = 0; k < _spatialDimension; k++ ){

                            dFinvdF[ _sd1 * A + a ][ _sd1 * k + K ] = -inverseDeformationGradient[ _sd1 * A + k ] * inverseDeformationGradient[ _sd1 * K + a ];

                        }

                    }

                }

            }

            dPK2dF          = floatMatrix( _sd2, floatVector( _sd2, 0 ) );
            dPK2dchi        = floatMatrix( _sd2, floatVector( _sd2, 0 ) );
            dPK2dgrad_chi   = floatMatrix( _sd2, floatVector( _sd3, 0 ) );

            dSIGMAdF        = floatMatrix( _sd2, floatVector( _sd2, 0 ) );
            dSIGMAdchi      = floatMatrix( _sd2, floatVector( _sd2, 0 ) );
            dSIGMAdgrad_chi = floatMatrix( _sd2, floatVector( _sd3, 0 ) );

            dMdF            = floatMatrix( _sd3, floatVector( _sd2, 0 ) );
            dMdchi          = floatMatrix( _sd3, floatVector( _sd2, 0 ) );
            dMdgrad_chi     = floatMatrix( _sd3, floatVector( _sd3, 0 ) );

        }

        // Compute the stress in the central section

        for ( unsigned int i = 0; i < _spatialDimension; i++ ){

            for ( unsigned int c = 0; c < _spatialDimension; c++ ){

                for ( unsigned int K = 0; K < _spatialDimension; K++ ){

                    for ( unsigned int L = 0; L < _spatialDimension; L++ ){

                        for ( unsigned int _M = 0; _M < _spatialDimension; _M++ ){

                            for ( unsigned int N = 0; N < _spatialDimension; N++ ){

                                for ( unsigned int A = 0; A < _spatialDimension; A++ ){

                                    for ( unsigned int B = 0; B < _spatialDimension; B++ ){

                                        SIGMA[ _sd1 * A + B ] += inverseDeformationGradient[ _sd1 * A + c ] * ( 0.5 * ( microDeformation[ _sd1 * i + L ] * microDeformation[ _sd1 * c + K ] + microDeformation[ _sd1 * i + K ] * microDeformation[ _sd1 * c + L ] ) * C_center[ _sd3 * K + _sd2 * L + _sd1 * _M + N ] * microStrain[ _sd1 * _M + N ] ) * inverseDeformationGradient[ _sd1 * B + i ];

                                        if ( computeJacobians ){

                                            for ( unsigned int r = 0; r < _spatialDimension; r++ ){

                                                for ( unsigned int R = 0; R < _spatialDimension; R++ ){

                                                    dSIGMAdF[ _sd1 * A + B ][ _sd1 * r + R ] += dFinvdF[ _sd1 * A + c ][ _sd1 * r + R ] * ( 0.5 * ( microDeformation[ _sd1 * i + L ] * microDeformation[ _sd1 * c + K ] + microDeformation[ _sd1 * i + K ] * microDeformation[ _sd1 * c + L ] ) * C_center[ _sd3 * K + _sd2 * L + _sd1 * _M + N ] * microStrain[ _sd1 * _M + N ] ) * inverseDeformationGradient[ _sd1 * B + i ]
                                                                                              + inverseDeformationGradient[ _sd1 * A + c ] * ( 0.5 * ( microDeformation[ _sd1 * i + L ] * microDeformation[ _sd1 * c + K ] + microDeformation[ _sd1 * i + K ] * microDeformation[ _sd1 * c + L ] ) * C_center[ _sd3 * K + _sd2 * L + _sd1 * _M + N ] * microStrain[ _sd1 * _M + N ] ) * dFinvdF[ _sd1 * B + i ][ _sd1 * r + R ];

                                                    dSIGMAdchi[ _sd1 * A + B ][ _sd1 * r + R ] += inverseDeformationGradient[ _sd1 * A + c ] * ( 0.5 * ( eye[ _sd1 * i + r] * eye[ _sd1 * L + R ] * microDeformation[ _sd1 * c + K ] + microDeformation[ _sd1 * i + L ] * eye[ _sd1 * c + r ] * eye[ _sd1 * K + R ] + eye[ _sd1 * i + r ] * eye[ _sd1 * K + R ] * microDeformation[ _sd1 * c + L ] + microDeformation[ _sd1 * i + K ] * eye[ _sd1 * c + r ] * eye[ _sd1 * L + R ] ) * C_center[ _sd3 * K + _sd2 * L + _sd1 * _M + N ] * microStrain[ _sd1 * _M + N ] ) * inverseDeformationGradient[ _sd1 * B + i ]
                                                                                                + inverseDeformationGradient[ _sd1 * A + c ] * ( 0.5 * ( microDeformation[ _sd1 * i + L ] * microDeformation[ _sd1 * c + K ] + microDeformation[ _sd1 * i + K ] * microDeformation[ _sd1 * c + L ] ) * C_center[ _sd3 * K + _sd2 * L + _sd1 * _M + N ] * dMicroStraindchi[ _sd1 * _M + N ][ _sd1 * r + R ] ) * inverseDeformationGradient[ _sd1 * B + i ];

                                                }

                                            }

                                        }

                                    }

                                }

                            }

                        }

                    }

                }

            }

        }

        // Compute the stress because of the connecting springs

        for ( auto p = params.begin( ) + 1; p != params.end( ); p++ ){

            if ( p->size( ) != 10 ){

                std::ostringstream message;
                message << "The parameters for the " << ( int )( p - params.begin( ) ) + 1
                        << "th spring has a size of " << p->size( ) << " and 10 values are required.";

                return new errorNode( __func__, message.str( ) );

            }

            // Extract the parameters
            floatType kn     = ( *p )[ 0 ];
            floatVector dX   = { ( *p )[ 1 ], ( *p )[ 2 ], ( *p )[ 3 ] };
            floatVector Xi_1 = { ( *p )[ 4 ], ( *p )[ 5 ], ( *p )[ 6 ] };
            floatVector Xi_2 = { ( *p )[ 7 ], ( *p )[ 8 ], ( *p )[ 9 ] };

            floatVector dx( _spatialDimension, 0 );
            floatVector xi_1( _spatialDimension, 0 );
            floatVector xi_2( _spatialDimension, 0 );

            for ( unsigned int i = 0; i < _spatialDimension; i++ ){

                for ( unsigned int I = 0; I < _spatialDimension; I++ ){

                    dx[ i ]   += deformationGradient[ 3 * i + I ] * dX[ I ];
                    xi_1[ i ] += microDeformation[ 3 * i + I ] * Xi_1[ I ];
                    xi_2[ i ] += microDeformation[ 3 * i + I ] * Xi_2[ I ];

                    for ( unsigned int J = 0; J < _spatialDimension; J++ ){

                        xi_2[ i ] += gradientMicroDeformation[ 9 * i + 3 * I + J ] * Xi_2[ I ] * dX[ J ];

                    }

                }

            }

            // Assemble the linear spring lengths
            floatVector C = dX + Xi_2 - Xi_1;
            floatVector c = dx + xi_2 - xi_1;

            // Compute the force value
            floatType fn = 0.5 * kn * ( 1 - std::sqrt( vectorTools::dot( C, C ) ) / std::sqrt( vectorTools::dot( c, c ) ) );

            floatMatrix ddxdF;
            floatMatrix dxi_1dchi;
            floatMatrix dxi_2dchi, dxi_2dgrad_chi;
            floatMatrix dcdF, dcdchi, dcdgrad_chi;
            floatVector dfdF, dfdchi, dfdgrad_chi;

            if ( computeJacobians ){

                ddxdF          = floatMatrix( _spatialDimension, floatVector( _sd2, 0 ) );
                dxi_1dchi      = floatMatrix( _spatialDimension, floatVector( _sd2, 0 ) );
                dxi_2dchi      = floatMatrix( _spatialDimension, floatVector( _sd2, 0 ) );
                dxi_2dgrad_chi = floatMatrix( _spatialDimension, floatVector( _sd3, 0 ) );

                for ( unsigned int i = 0; i < _spatialDimension; i++ ){

                    for ( unsigned int j = 0; j < _spatialDimension; j++ ){

                        for ( unsigned int J = 0; J < _spatialDimension; J++ ){

                            ddxdF[ i ][ _sd1 * j + J ] = eye[ _sd1 * i + j ] * dX[ J ];

                            dxi_1dchi[ i ][ _sd1 * j + J ] = eye[ _sd1 * i + j ] * Xi_1[ J ];

                            dxi_2dchi[ i ][ _sd1 * j + J ] = eye[ _sd1 * i + j ] * Xi_2[ J ];

                            for ( unsigned int K = 0; K < _spatialDimension; K++ ){

                                dxi_2dgrad_chi[ i ][ _sd2 * j + _sd1 * J + K ] = eye[ _sd1 * i + j ] * Xi_2[ J ] * dX[ K ];

                            }

                        }

                    }

                }

                // Assemble the Jacobians of the spring vectors

                dcdF        = ddxdF;
                dcdchi      = dxi_2dchi - dxi_1dchi;
                dcdgrad_chi = dxi_2dgrad_chi;

                // Assemble the Jacobians of the force values

                dfdF        = floatVector( _sd2, 0 );
                dfdchi      = floatVector( _sd2, 0 );
                dfdgrad_chi = floatVector( _sd3, 0 );

                floatType factor = 0.5 * kn * std::sqrt( vectorTools::dot( C, C ) ) * std::pow( vectorTools::dot( c, c ), -1.5 );

                for ( unsigned int k = 0; k < _spatialDimension; k++ ){

                    for ( unsigned int K = 0; K < _spatialDimension; K++ ){

                        dfdF[ _sd1 * k + K ] += factor * c[ k ] * dX[ K ];

                        dfdchi[ _sd1 * k + K ] += factor * c[ k ] * ( Xi_2[ K ] - Xi_1[ K ] );

                        for ( unsigned int L = 0; L < _spatialDimension; L++ ){

                            dfdgrad_chi[ _sd2 * k + _sd1 * K + L ] = factor * c[ k ] * Xi_2[ K ] * dX[ L ];

                        }

                    }

                }

            }

            // Add the contributions of the current spring to the stresses

            for ( unsigned int a = 0; a < _sd1; a++ ){

                for ( unsigned int A = 0; A < _sd1; A++ ){

                    for ( unsigned int B = 0; B < _sd1; B++ ){

                        PK2[ _sd1 * A + B ] += fn * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * c[ a ];

                        if ( computeJacobians ){

                            for ( unsigned int r = 0; r < _sd1; r++ ){

                                for ( unsigned int R = 0; R < _sd1; R++ ){

                                    dPK2dF[ _sd1 * A + B ][ _sd1 * r + R ] += dfdF[ _sd1 * r + R ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * c[ a ]
                                                                            + fn * dX[ A ] * dFinvdF[ _sd1 * B + a ][ _sd1 * r + R ] * c[ a ]
                                                                            + fn * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * dcdF[ a ][ _sd1 * r + R ];

                                    dPK2dchi[ _sd1 * A + B ][ _sd1 * r + R ] += dfdchi[ _sd1 * r + R ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * c[ a ]
                                                                              + fn * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * dcdchi[ a ][ _sd1 * r + R ];

                                    for ( unsigned int S = 0; S < _sd1; S++ ){

                                        dPK2dgrad_chi[ _sd1 * A + B ][ _sd2 * r + _sd1 * R + S ] += dfdgrad_chi[ _sd2 * r + _sd1 * R + S ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * c[ a ]
                                                                                                  + fn * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * dcdgrad_chi[ a ][ _sd2 * r + _sd1 * R + S ];

                                    }

                                }

                            }

                        }

                        for ( unsigned int b = 0; b < _sd1; b++ ){

                            SIGMA[ _sd1 * A + B ] += fn * c[ a ] * c[ b ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ];

                            M[ _sd2 * A + _sd1 * B + b ] += fn * c[ a ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * Xi_2[ b ];

                            if ( computeJacobians ){

                                for ( unsigned int r = 0; r < _sd1; r++ ){

                                    for ( unsigned int R = 0; R < _sd1; R++ ){

                                        dSIGMAdF[ _sd1 * A + B ][ _sd1 * r + R ] += dfdF[ _sd1 * r + R ] * c[ a ] * c[ b ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ]
                                                                                  + fn * dcdF[ a ][ _sd1 * r + R ] * c[ b ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ]
                                                                                  + fn * c[ a ] * dcdF[ b ][ _sd1 * r + R ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ]
                                                                                  + fn * c[ a ] * c[ b ] * dFinvdF[ _sd1 * B + b ][ _sd1 * r + R ] * inverseDeformationGradient[ _sd1 * A + a ]
                                                                                  + fn * c[ a ] * c[ b ] * inverseDeformationGradient[ _sd1 * B + b ] * dFinvdF[ _sd1 * A + a ][ _sd1 * r + R ];

                                        dMdF[ _sd2 * A + _sd1 * B + b ][ _sd1 * r + R ] += dfdF[ _sd1 * r + R ] * c[ a ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * Xi_2[ b ]
                                                                                         + fn * dcdF[ a ][ _sd1 * r + R ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * Xi_2[ b ]
                                                                                         + fn * c[ a ] * dX[ A ] * dFinvdF[ _sd1 * B + a ][ _sd1 * r + R ] * Xi_2[ b ];

                                        dSIGMAdchi[ _sd1 * A + B ][ _sd1 * r + R ] += dfdchi[ _sd1 * r + R ] * c[ a ] * c[ b ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ]
                                                                                    + fn * dcdchi[ a ][ _sd1 * r + R ] * c[ b ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ]
                                                                                    + fn * c[ a ] * dcdchi[ b ][ _sd1 * r + R ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ];

                                        dMdchi[ _sd2 * A + _sd1 * B + b ][ _sd1 * r + R ] += dfdchi[ _sd1 * r + R ] * c[ a ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * Xi_2[ b ]
                                                                                           + fn * dcdchi[ a ][ _sd1 * r + R ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * Xi_2[ b ];


                                        for ( unsigned int S = 0; S < _sd1; S++ ){

                                            dSIGMAdgrad_chi[ _sd1 * A + B ][ _sd2 * r + _sd1 * R + S ] += dfdgrad_chi[ _sd2 * r + _sd1 * R + S ] * c[ a ] * c[ b ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ]
                                                                                    + fn * dcdgrad_chi[ a ][ _sd2 * r + _sd1 * R + S ] * c[ b ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ]
                                                                                    + fn * c[ a ] * dcdgrad_chi[ b ][ _sd2 * r + _sd1 * R + S ] * inverseDeformationGradient[ _sd1 * B + b ] * inverseDeformationGradient[ _sd1 * A + a ];

                                            dMdgrad_chi[ _sd2 * A + _sd1 * B + b ][ _sd2 * r + _sd1 * R + S ] += dfdgrad_chi[ _sd2 * r + _sd1 * R + S] * c[ a ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * Xi_2[ b ]
                                                                                                               + fn * dcdgrad_chi[ a ][ _sd2 * r + _sd1 * R + S ] * dX[ A ] * inverseDeformationGradient[ _sd1 * B + a ] * Xi_2[ b ];

                                        }

                                    }

                                }

                            }

                        }

                    }

                }

            }

        }

        return NULL;

    }

}
