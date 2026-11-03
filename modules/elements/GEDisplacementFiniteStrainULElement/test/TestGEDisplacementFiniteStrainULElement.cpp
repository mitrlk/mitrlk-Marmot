#include "Marmot/GEDisplacementFiniteStrainULElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotNumericalDifferentiation.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Elements;
using namespace Marmot::Testing;
namespace NumDiff = Marmot::NumericalAlgorithms::Differentiation;

// One increment of DufourModel from rest, far in the plastic range. In plastic flow dL_dF is not zero, so the whole
// stiffness, the nonlocal rows and columns included, must match central differences of the residual.
// computeYourself() returns the negative residual Pe = -P_int and K = dP_int / dQ = -dPe / dQ.
template < class Element >
void checkTangentAgainstNumericalDifferentiation( const std::string&            label,
                                                  const std::vector< double >&  nodeCoordsVec,
                                                  typename Element::SectionType sectionType )
{
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration;

  // SikaPower-498 card of the DufourModel test (K, G, ft, fc, Q1, Q2, b1, b2, b3, b4, b5, eta_VP, n, nuP+, nuP-,
  // epsF, omegaMax, ld, m, density), failure driver off
  const std::vector< double > matProps   = { 1793.9, 828.0, 12.0,   21.6, 4051.0, 15.8, 552.0, 243.0, 0.0, 0.0,
                                             15.4,   1e-6,  0.0435, 0.3,  0.5,    1.75, 0.99,  0.12,  1.0, 1.0 };
  const std::vector< double > elPropsVec = { 1.0 }; // thickness, read by initializeYourself() for nDim = 2

  auto element = std::make_unique< Element >( 1, intType, sectionType );
  element->assignNodeCoordinates( nodeCoordsVec.data() );
  element->assignProperty(
    MarmotMaterialSection( "DUFOURMODEL", matProps.data(), static_cast< int >( matProps.size() ) ) );
  element->assignProperty( ElementProperties( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) ) );

  const int             nStateVars = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVars, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVars );
  element->initializeYourself();
  element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  const std::vector< double > stateVarsAtRest = stateVars;

  const int    nDof   = element->getNDofPerElement();
  const int    nU     = Element::bsU;
  const double time[] = { 0.0, 0.0 };
  const double dT     = 1.0;

  // displacements of up to 0.04 (strains of several times the yield strain of about 0.015) and a small positive
  // nonlocal field
  Eigen::VectorXd Q0( nDof );
  for ( int i = 0; i < nU; i++ )
    Q0[i] = 0.04 * std::sin( 1.7 * ( i + 1 ) );
  for ( int i = nU; i < nDof; i++ )
    Q0[i] = 1e-3 * ( 1.0 + 0.1 * ( i - nU ) );

  MarmotElement& asElement = *element; // computeYourself() is private in the axisymmetric class

  // negative residual and stiffness of one increment from rest, at the displacement and nonlocal field Q
  auto evaluate = [&]( const Eigen::VectorXd& Q, Eigen::MatrixXd& K ) -> Eigen::VectorXd {
    stateVars          = stateVarsAtRest;
    Eigen::VectorXd Pe = Eigen::VectorXd::Zero( nDof );
    K.setZero( nDof, nDof );
    double pNewdT = 1.0;
    asElement.computeYourself( Q.data(), Q.data(), Pe.data(), K.data(), time, dT, pNewdT );
    throwExceptionOnFailure( pNewdT >= 1.0, label + ": computeYourself() requested a cutback." );
    return Pe;
  };

  Eigen::MatrixXd K( nDof, nDof );
  evaluate( Q0, K );

  // the increment must be plastic: in an elastic one dL_dF is zero, and with it the block AU
  const double maxKAU = K.block( nU, 0, nDof - nU, nU ).cwiseAbs().maxCoeff();
  throwExceptionOnFailure( maxKAU > 1e-6,
                           label + ": the increment is not plastic, max |K_AU| " + std::to_string( maxKAU ) );

  NumDiff::vector_to_vector_function_type negativeResidual = [&]( const Eigen::VectorXd& Q ) -> Eigen::VectorXd {
    Eigen::MatrixXd Ktmp( nDof, nDof );
    return -evaluate( Q, Ktmp );
  };
  const Eigen::MatrixXd numK = NumDiff::centralDifference( negativeResidual, Q0 );

  // block by block, relative to the largest entry of each block: UU, UA, AU, AA
  const int   blockStart[2] = { 0, nU }, blockSize[2] = { nU, nDof - nU };
  std::string failures;
  for ( int a = 0; a < 2; a++ )
    for ( int b = 0; b < 2; b++ ) {
      const Eigen::MatrixXd Kb    = K.block( blockStart[a], blockStart[b], blockSize[a], blockSize[b] );
      const Eigen::MatrixXd numKb = numK.block( blockStart[a], blockStart[b], blockSize[a], blockSize[b] );
      const double          scale = std::max( numKb.cwiseAbs().maxCoeff(), 1e-12 );
      const double          error = ( Kb - numKb ).cwiseAbs().maxCoeff() / scale;
      if ( !( error < 1e-5 ) )
        failures += std::string( " " ) + ( a == 0 ? "U" : "A" ) + ( b == 0 ? "U" : "A" ) + " " +
                    std::to_string( error );
    }
  throwExceptionOnFailure( failures.empty(),
                           label +
                             ": stiffness blocks differ from central differences (block, relative error):" + failures );
}

// quad8, corners then the midsides of the edges 0-1, 1-2, 2-3, 3-0; for the axisymmetric element the ring r in [1, 2],
// z in [0, 1], away from the axis since the hoop terms divide by r
const std::vector< double > quad8Ring =
  { 1.0, 0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 1.0, 1.5, 0.0, 2.0, 0.5, 1.5, 1.0, 1.0, 0.5 };

// GCX8RUL. The hoop part of the coupling block dL / dq_U is dL_dF( 2, 2 ) N / r, once per node, and the hoop force
// depends on the nonlocal field through dTau_dA( 2, 2 ).
void testAxiSymmetricGETangentMatchesNumericalDifferentiation()
{
  using Element = AxiSymmetricGEDisplacementFiniteStrainULElement< 8 >;
  checkTangentAgainstNumericalDifferentiation< Element >( "GCX8RUL", quad8Ring, Element::PlaneStrain );
}

// GCPE8RUL, the same element without the hoop terms
void testPlaneStrainGETangentMatchesNumericalDifferentiation()
{
  using Element = GEDisplacementFiniteStrainULElement< 2, 8 >;
  checkTangentAgainstNumericalDifferentiation< Element >( "GCPE8RUL", quad8Ring, Element::PlaneStrain );
}

// GC3D20RUL on the unit cube
void testSolidGETangentMatchesNumericalDifferentiation()
{
  using Element                      = GEDisplacementFiniteStrainULElement< 3, 20 >;
  const std::vector< double > hexa20 = { 0,  0, 0,  1, 0,  0, 1,  1, 0,  0, 1,  0, 0, 0,  1,  1, 0,  1, 1, 1,
                                         1,  0, 1,  1, .5, 0, 0,  1, .5, 0, .5, 1, 0, 0,  .5, 0, .5, 0, 1, 1,
                                         .5, 1, .5, 1, 1,  0, .5, 1, 0,  0, .5, 1, 0, .5, 1,  1, .5, 0, 1, .5 };
  checkTangentAgainstNumericalDifferentiation< Element >( "GC3D20RUL", hexa20, Element::Solid );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testAxiSymmetricGETangentMatchesNumericalDifferentiation,
    testPlaneStrainGETangentMatchesNumericalDifferentiation,
    testSolidGETangentMatchesNumericalDifferentiation,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
