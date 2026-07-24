#include "Marmot/DufourModel.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <array>
#include <functional>
#include <string>
#include <vector>

using namespace Marmot::Testing;
using namespace Marmot::Materials;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

namespace {

  // SikaPower-498 parameter card (Dufour), order as read by DufourModel:
  // K, G, ft, fc, Q1, Q2, b1, b2, b3, b4, b5, eta_VP, n, nuP+, nuP-, epsF, omegaMax, ld, m, density
  const std::array< double, 20 > props = { 1793.9, 828.0, 12.0, 21.6, 4051.0, 15.8, 552.0, 243.0, 0.0, 0.0,
                                           15.4,   1e-6,  0.0435, 0.3, 0.5,    1.75, 0.99,  0.12,  1.0, 1.0 };
  const int                      elLabel = 1;

  // Evaluate the Kirchhoff-stress response (and, optionally, the analytic tangent) for a given
  // deformation gradient F and nonlocal field A, always starting from a fresh (undeformed) state.
  DufourModel::ConstitutiveResponse< 3 > evaluate( const Tensor33d&                    F,
                                                   double                              A,
                                                   DufourModel::AlgorithmicModuli< 3 >* tangentOut = nullptr )
  {
    DufourModel         mat( props.data(), static_cast< int >( props.size() ), elLabel );
    std::vector< double > sv( mat.getNumberOfRequiredStateVars(), 0.0 );
    mat.assignStateVars( sv.data(), static_cast< int >( sv.size() ) );
    mat.initializeYourself();

    DufourModel::Deformation< 3 >          def{ F, A };
    DufourModel::TimeIncrement             timeInc{ 0.0, 1.0 };
    DufourModel::ConstitutiveResponse< 3 > response{ Tensor33d( 0.0 ), 0.0, 0.0, 0.0, 0.0 };
    DufourModel::AlgorithmicModuli< 3 >    tangent{ Tensor3333d( 0.0 ), Tensor33d( 0.0 ), Tensor33d( 0.0 ) };

    mat.computeStress( response, tangent, def, timeInc );
    if ( tangentOut )
      *tangentOut = tangent;
    return response;
  }

} // namespace

// I-1: undeformed configuration -> zero Kirchhoff stress
void testUndeformedResponse()
{
  const auto r = evaluate( Spatial3D::I, 0.0 );
  throwExceptionOnFailure( checkIfEqual( r.tau, Tensor33d( 0.0 ), 1e-9 ),
                           "DufourModel I-1: undeformed configuration must yield zero Kirchhoff stress in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// I-2: Kirchhoff stress must be symmetric under an arbitrary (small) deformation
void testStressSymmetry()
{
  Tensor33d F = Spatial3D::I;
  F( 0, 0 ) += 0.004;
  F( 1, 1 ) -= 0.002;
  F( 2, 2 ) += 0.003;
  F( 0, 1 ) += 0.0015;
  F( 1, 0 ) += 0.0015;
  F( 0, 2 ) -= 0.0010;
  F( 2, 0 ) -= 0.0010;

  const auto r = evaluate( F, 0.0 );
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      throwExceptionOnFailure( checkIfEqual( r.tau( i, j ), r.tau( j, i ), 1e-9 ),
                               "DufourModel I-2: Kirchhoff stress must be symmetric in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
}

// I-3: pure rigid-body rotation -> zero stress (frame indifference)
void testPureRotation()
{
  for ( int phi_deg = 0; phi_deg <= 180; phi_deg += 30 ) {
    const double phi = Marmot::Math::degToRad( phi_deg );
    Tensor33d    F( 0.0 );
    F( 0, 0 ) = cos( phi );
    F( 0, 1 ) = -sin( phi );
    F( 1, 0 ) = sin( phi );
    F( 1, 1 ) = cos( phi );
    F( 2, 2 ) = 1.0;

    const auto r = evaluate( F, 0.0 );
    throwExceptionOnFailure( checkIfEqual( r.tau, Tensor33d( 0.0 ), 1e-8 ),
                             "DufourModel I-3: pure rotation (phi_deg=" + std::to_string( phi_deg ) +
                               ") must yield zero stress in " + std::string( __PRETTY_FUNCTION__ ) );
  }
}

// I-4: analytic dTau/dF must match central finite differences (elastic regime)
void testAlgorithmicTangent()
{
  Tensor33d F0 = Spatial3D::I;
  F0( 0, 0 ) += 0.0004;
  F0( 1, 1 ) += 0.0002;
  F0( 2, 2 ) -= 0.0001;
  F0( 0, 1 ) += 0.0003;
  F0( 1, 0 ) += 0.0003;

  DufourModel::AlgorithmicModuli< 3 > tangent;
  evaluate( F0, 0.0, &tangent );

  const double eps = 1e-6;
  Tensor3333d  dTau_dF_fd( 0.0 );
  for ( int k = 0; k < 3; k++ )
    for ( int l = 0; l < 3; l++ ) {
      Tensor33d Fp = F0, Fm = F0;
      Fp( k, l ) += eps;
      Fm( k, l ) -= eps;
      const auto rp = evaluate( Fp, 0.0 );
      const auto rm = evaluate( Fm, 0.0 );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          dTau_dF_fd( i, j, k, l ) = ( rp.tau( i, j ) - rm.tau( i, j ) ) / ( 2.0 * eps );
    }

  throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, dTau_dF_fd, 1e-3 ),
                           "DufourModel I-4: analytic dTau/dF must match central finite differences in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testUndeformedResponse,
    testStressSymmetry,
    testPureRotation,
    testAlgorithmicTangent,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
