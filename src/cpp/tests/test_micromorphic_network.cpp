//Tests for micromorphic_network

#include<micromorphic_network.h>
#include<sstream>
#include<fstream>
#include<iostream>

#define BOOST_TEST_MODULE test_micromorphic_network
#include <boost/test/included/unit_test.hpp>

typedef microNet::errorOut errorOut;
typedef microNet::errorNode errorNode;
typedef microNet::floatType floatType;
typedef microNet::floatVector floatVector;
typedef microNet::floatMatrix floatMatrix;

BOOST_AUTO_TEST_CASE( test_assembleStiffness ){
    /*!
     * Test of the assembly of the central mass' stiffness matrix
     */

    floatType lambda = 7.29765590446705;
    floatType mu     = 2.506265664160401;

    floatMatrix params = {
                             { lambda, mu },
                             { 2.,  1.0,  0.0,  0.0,  0.2,  0.2,  0.0, -0.2, -0.2,  0.0 },
                             { 2., -1.0,  0.0,  0.0, -0.2, -0.2,  0.0,  0.2,  0.2,  0.0 },
                             { 2.,  0.0,  1.0,  0.0, -0.2,  0.2,  0.0,  0.2, -0.2,  0.0 },
                             { 2.,  0.0, -1.0,  0.0,  0.2, -0.2,  0.0, -0.2,  0.2,  0.0 },
                             { 2.,  0.0,  0.0,  1.0,  0.0,  0.0,  0.2,  0.0,  0.0, -0.2 },
                             { 2.,  0.0,  0.0, -1.0,  0.0,  0.0, -0.2,  0.0,  0.0,  0.2 }
                         };

    floatVector answer = { 12.31018723,  0.        ,  0.        ,  0.        ,  7.2976559 ,
                            0.        ,  0.        ,  0.        ,  7.2976559 ,  0.        ,
                            2.50626566,  0.        ,  2.50626566,  0.        ,  0.        ,
                            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                            2.50626566,  0.        ,  0.        ,  0.        ,  2.50626566,
                            0.        ,  0.        ,  0.        ,  2.50626566,  0.        ,
                            2.50626566,  0.        ,  0.        ,  0.        ,  0.        ,
                            0.        ,  7.2976559 ,  0.        ,  0.        ,  0.        ,
                           12.31018723,  0.        ,  0.        ,  0.        ,  7.2976559 ,
                            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                            2.50626566,  0.        ,  2.50626566,  0.        ,  0.        ,
                            0.        ,  2.50626566,  0.        ,  0.        ,  0.        ,
                            2.50626566,  0.        ,  0.        ,  0.        ,  0.        ,
                            0.        ,  0.        ,  0.        ,  2.50626566,  0.        ,
                            2.50626566,  0.        ,  7.2976559 ,  0.        ,  0.        ,
                            0.        ,  7.2976559 ,  0.        ,  0.        ,  0.        ,
                           12.31018723 };

    floatVector result;

    errorOut error = microNet::assembleStiffness( params, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( test_computeEnergy ){
    /*!
     * Test the computation of the Helmholtz free energy.
     */

    floatType lambda = 7.29765590446705;
    floatType mu     = 2.506265664160401;

    floatMatrix params = {
                             { lambda, mu },
                             { 2.,  1.0,  0.0,  0.0,  0.2,  0.2,  0.0, -0.2, -0.2,  0.0 },
                             { 2., -1.0,  0.0,  0.0, -0.2, -0.2,  0.0,  0.2,  0.2,  0.0 },
                             { 2.,  0.0,  1.0,  0.0, -0.2,  0.2,  0.0,  0.2, -0.2,  0.0 },
                             { 2.,  0.0, -1.0,  0.0,  0.2, -0.2,  0.0, -0.2,  0.2,  0.0 },
                             { 2.,  0.0,  0.0,  1.0,  0.0,  0.0,  0.2,  0.0,  0.0, -0.2 },
                             { 2.,  0.0,  0.0, -1.0,  0.0,  0.0, -0.2,  0.0,  0.0,  0.2 }
                         };

    floatVector F = { 1.04571166, 0.01456296, 0.02111255, 0.03342643, 1.09937395,
                      0.0060523 , 0.02281969, 0.09017899, 1.01420422 };

    floatVector chi = { 1.06412257, 0.053143  , 0.01451306, 0.06457342, 1.04782925,
                        0.05113791, 0.04144923, 0.0172515 , 1.09107831 };

    floatVector grad_chi = { 0.35617972, 0.5451236 , 0.98177534, 0.33615045, 0.760854  ,
                             0.47847324, 0.28144915, 0.24489963, 0.89853014, 0.24501212,
                             0.41851747, 0.09278363, 0.63127664, 0.19385062, 0.39255271,
                             0.99951862, 0.96871357, 0.10916816, 0.03280607, 0.99890606,
                             0.30816764, 0.8519894 , 0.05398674, 0.99729743, 0.2978061 ,
                             0.67015276, 0.86704963 };

    floatType answer = 0.2942345813887143;

    floatType result;

    errorOut error = microNet::computeEnergy( F, chi, grad_chi, params, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( test_computeStress ){
    /*!
     * Test the computation of the stresses
     */

    floatType lambda = 7.29765590446705;
    floatType mu     = 2.506265664160401;

    floatMatrix params = {
                             { lambda, mu },
                             { 2.,  1.0,  0.0,  0.0,  0.2,  0.2,  0.0, -0.2, -0.2,  0.0 },
                             { 2., -1.0,  0.0,  0.0, -0.2, -0.2,  0.0,  0.2,  0.2,  0.0 },
                             { 2.,  0.0,  1.0,  0.0, -0.2,  0.2,  0.0,  0.2, -0.2,  0.0 },
                             { 2.,  0.0, -1.0,  0.0,  0.2, -0.2,  0.0, -0.2,  0.2,  0.0 },
                             { 2.,  0.0,  0.0,  1.0,  0.0,  0.0,  0.2,  0.0,  0.0, -0.2 },
                             { 2.,  0.0,  0.0, -1.0,  0.0,  0.0, -0.2,  0.0,  0.0,  0.2 }
                         };

    floatVector F = { 1.04571166, 0.01456296, 0.02111255, 0.03342643, 1.09937395,
                      0.0060523 , 0.02281969, 0.09017899, 1.01420422 };

    floatVector chi = { 1.06412257, 0.053143  , 0.01451306, 0.06457342, 1.04782925,
                        0.05113791, 0.04144923, 0.0172515 , 1.09107831 };

    floatVector grad_chi = { 0.35617972, 0.5451236 , 0.98177534, 0.33615045, 0.760854  ,
                             0.47847324, 0.28144915, 0.24489963, 0.89853014, 0.24501212,
                             0.41851747, 0.09278363, 0.63127664, 0.19385062, 0.39255271,
                             0.99951862, 0.96871357, 0.10916816, 0.03280607, 0.99890606,
                             0.30816764, 0.8519894 , 0.05398674, 0.99729743, 0.2978061 ,
                             0.67015276, 0.86704963 };

    floatVector PK2_answer = { 0.09149356, -0.05428376,  0.00832308,  0.1136431 ,  0.19184733,
                               0.02413619,  0.09996321,  0.01058891,  0.01742596 };

    floatVector SIGMA_answer = { 2.11344778, 0.47086981, 0.23615829, 0.47086981, 1.84329275,
                                 0.18956751, 0.23615829, 0.18956751, 2.43694791 };

    floatVector M_answer = { 0.00758761,  0.00758761,  0.        ,  0.00213902,  0.00213902,
                             0.        ,  0.00491321,  0.00491321,  0.        ,  0.00326427,
                            -0.00326427,  0.        ,  0.012122  , -0.012122  ,  0.        ,
                             0.01151558, -0.01151558,  0.        ,  0.        ,  0.        ,
                            -0.00441271,  0.        ,  0.        , -0.00234627,  0.        ,
                             0.        ,  0.06341387 };

    floatVector PK2_result, SIGMA_result, M_result;

    errorOut error = microNet::computeStress( F, chi, grad_chi, params, PK2_result, SIGMA_result, M_result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2_result, PK2_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( M_result, M_answer ) );

}
