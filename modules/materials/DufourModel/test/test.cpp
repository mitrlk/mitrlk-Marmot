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
  const std::array< double, 20 > props   = { 1793.9, 828.0, 12.0,   21.6, 4051.0, 15.8, 552.0, 243.0, 0.0, 0.0,
                                             15.4,   1e-6,  0.0435, 0.3,  0.5,    1.75, 0.99,  0.12,  1.0, 1.0 };
  const int                      elLabel = 1;

  // Evaluate the Kirchhoff-stress response (and, optionally, the analytic tangent) for a given
  // deformation gradient F and nonlocal field A, always starting from a fresh (undeformed) state.
  // SWDFM/paraboloidal card: base 20 entries + cSW, bSW, kSW, dF, volDriver, Xt, Xc.
  // Xt/Xc active -> the stress-driven omega_f branch and its tangent contribution are exercised.
  const std::vector< double > propsStress = { 1793.9, 828.0, 12.0, 21.6,   4051.0, 15.8, 552.0, 243.0, 0.0,
                                              0.0,    15.4,  1e-6, 0.0435, 0.3,    0.5,  1.75,  0.99,  0.12,
                                              1.0,    1.0,   0.0,  1e6,    0.0,    0.1,  0.0,   8.0,   14.4 };

  // ---- viscoelasticity cards -----------------------------------------------------------------
  // Base 20 entries, then the damage-variant slots 20..26 switched OFF (cSW = 0, Xt = 0) so that
  // ONLY the classical omega = 1 - exp(-alphaP_bar/epsF) law is active, then:
  //   [27] nMaxwell, and per branch [28+3i] gammaG_i, [29+3i] gammaK_i, [30+3i] tau_i.
  const std::vector< double > propsBase = { 1793.9, 828.0, 12.0, 21.6, 4051.0, 15.8, 552.0, 243.0, 0.0, 0.0,
                                            15.4,   1e-6,  0.0435, 0.3, 0.5,   1.75, 0.99,  0.12,  1.0, 1.0,
                                            0.0,    1e6,   0.0,  0.1,  0.0,    0.0,  0.0 };

  // two branches, deliberately unequal shear/bulk weights and well-separated relaxation times
  std::vector< double > withMaxwell( std::vector< double > gammaGKTau )
  {
    std::vector< double > card = propsBase;
    card.push_back( static_cast< double >( static_cast< int >( gammaGKTau.size() ) / 3 ) );
    card.insert( card.end(), gammaGKTau.begin(), gammaGKTau.end() );
    return card;
  }

  const std::vector< double > propsVisco = withMaxwell( { 0.30, 0.15, 0.5, 0.20, 0.10, 5.0 } );

  // Advance the material by ONE increment, continuing from (and updating in place) the supplied
  // state vector. An empty vector is sized and initialised first, i.e. starts from rest.
  DufourModel::ConstitutiveResponse< 3 > step( std::vector< double >&               sv,
                                               const Tensor33d&                     F,
                                               double                               A,
                                               double                               dT,
                                               const std::vector< double >&         card,
                                               DufourModel::AlgorithmicModuli< 3 >* tangentOut = nullptr )
  {
    DufourModel mat( card.data(), static_cast< int >( card.size() ), elLabel );

    const bool fresh = sv.empty();
    if ( fresh )
      sv.assign( mat.getNumberOfRequiredStateVars(), 0.0 );
    mat.assignStateVars( sv.data(), static_cast< int >( sv.size() ) );
    if ( fresh )
      mat.initializeYourself();

    DufourModel::Deformation< 3 >          def{ F, A };
    DufourModel::TimeIncrement             timeInc{ 0.0, dT };
    DufourModel::ConstitutiveResponse< 3 > response{ Tensor33d( 0.0 ), 0.0, 0.0, 0.0, 0.0 };
    DufourModel::AlgorithmicModuli< 3 >    tangent{ Tensor3333d( 0.0 ), Tensor33d( 0.0 ), Tensor33d( 0.0 ) };

    mat.computeStress( response, tangent, def, timeInc );
    if ( tangentOut )
      *tangentOut = tangent;
    return response;
  }

  // Central-difference dTau/dF of one increment taken from a FIXED starting history.
  Tensor3333d finiteDifferenceTangent( const std::vector< double >& svStart,
                                       const Tensor33d&             F0,
                                       double                       dT,
                                       const std::vector< double >& card,
                                       double                       eps )
  {
    Tensor3333d fd( 0.0 );
    for ( int k = 0; k < 3; k++ )
      for ( int l = 0; l < 3; l++ ) {
        Tensor33d Fp = F0, Fm = F0;
        Fp( k, l ) += eps;
        Fm( k, l ) -= eps;

        std::vector< double > svp = svStart, svm = svStart;
        const auto            rp = step( svp, Fp, 0.0, dT, card );
        const auto            rm = step( svm, Fm, 0.0, dT, card );

        for ( int i = 0; i < 3; i++ )
          for ( int j = 0; j < 3; j++ )
            fd( i, j, k, l ) = ( rp.tau( i, j ) - rm.tau( i, j ) ) / ( 2.0 * eps );
      }
    return fd;
  }

  DufourModel::ConstitutiveResponse< 3 > evaluate( const Tensor33d&                     F,
                                                   double                               A,
                                                   DufourModel::AlgorithmicModuli< 3 >* tangentOut = nullptr,
                                                   const std::vector< double >&         cardIn     = {} )
  {
    const std::vector< double > card = cardIn.empty() ? std::vector< double >( props.begin(), props.end() ) : cardIn;
    DufourModel                 mat( card.data(), static_cast< int >( card.size() ), elLabel );
    std::vector< double >       sv( mat.getNumberOfRequiredStateVars(), 0.0 );
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

// I-5: with the STRESS-DRIVEN omega_f active (Xt, Xc set low enough that phiBar > 1), the
// analytic dTau/dF must still match central finite differences. This checks the
// -tau_eff (x) ( dOmega/dTau : dTau_eff/dF ) contribution, which is easy to get wrong.
void testAlgorithmicTangentStressDrivenDamage()
{
  Tensor33d F0 = Spatial3D::I;
  F0( 0, 0 ) += 0.010;
  F0( 1, 1 ) -= 0.003;
  F0( 2, 2 ) -= 0.003;
  F0( 0, 1 ) += 0.002;
  F0( 1, 0 ) += 0.002;

  DufourModel::AlgorithmicModuli< 3 > tangent;
  const auto                          r0 = evaluate( F0, 0.0, &tangent, propsStress );
  if ( !( r0.tau( 0, 0 ) > 0.0 ) )
    throw std::runtime_error( "DufourModel I-5: test state carries no stress" );

  // Guard against the test becoming vacuous: omega_f must actually be ACTIVE in this state,
  // i.e. the stress-driven branch must change the response relative to Xt = 0.
  const auto rNoFail = evaluate( F0, 0.0, nullptr );
  if ( std::abs( r0.tau( 0, 0 ) - rNoFail.tau( 0, 0 ) ) < 1e-8 * std::abs( rNoFail.tau( 0, 0 ) ) )
    throw std::runtime_error( "DufourModel I-5: stress-driven omega_f is NOT active - "
                              "the tangent check would be vacuous" );

  const double eps = 1e-7;
  Tensor3333d  dTau_dF_fd( 0.0 );
  for ( int k = 0; k < 3; k++ )
    for ( int l = 0; l < 3; l++ ) {
      Tensor33d Fp = F0, Fm = F0;
      Fp( k, l ) += eps;
      Fm( k, l ) -= eps;
      const auto rp = evaluate( Fp, 0.0, nullptr, propsStress );
      const auto rm = evaluate( Fm, 0.0, nullptr, propsStress );
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          dTau_dF_fd( i, j, k, l ) = ( rp.tau( i, j ) - rm.tau( i, j ) ) / ( 2.0 * eps );
    }

  throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, dTau_dF_fd, 1e-2 ),
                           "DufourModel I-5: analytic dTau/dF with stress-driven damage must match "
                           "central finite differences in " +
                             std::string( __PRETTY_FUNCTION__ ) );
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

// V-1: with the generalized Maxwell chain ACTIVE, the analytic dTau/dF must match central finite
// differences on the FIRST increment (branch states still at rest).
void testViscoelasticTangentFirstIncrement()
{
  Tensor33d F0 = Spatial3D::I;
  F0( 0, 0 ) += 0.0004;
  F0( 1, 1 ) += 0.0002;
  F0( 2, 2 ) -= 0.0001;
  F0( 0, 1 ) += 0.0003;
  F0( 1, 0 ) += 0.0003;
  const double dT = 0.5;

  std::vector< double >               sv;
  DufourModel::AlgorithmicModuli< 3 > tangent;
  const auto                          r0 = step( sv, F0, 0.0, dT, propsVisco, &tangent );

  // Guard against a vacuous check: relaxation must actually move the stress.
  std::vector< double > svNoVE;
  const auto            rNoVE = step( svNoVE, F0, 0.0, dT, propsBase );
  if ( std::abs( r0.tau( 0, 0 ) - rNoVE.tau( 0, 0 ) ) < 1e-6 * std::abs( rNoVE.tau( 0, 0 ) ) )
    throw std::runtime_error( "DufourModel V-1: viscoelasticity is NOT active - "
                              "the tangent check would be vacuous" );

  const std::vector< double > svRest( sv.size(), 0.0 );
  const Tensor3333d           fd = finiteDifferenceTangent( {}, F0, dT, propsVisco, 1e-6 );

  throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, fd, 1e-3 ),
                           "DufourModel V-1: analytic dTau/dF with viscoelasticity must match "
                           "central finite differences in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// V-2: the RECURSIVE path. Take one increment first so the Maxwell branches carry a non-zero
// history, then verify the tangent of the SECOND increment against finite differences taken from
// that same history. This is what catches a mishandled commit/scratch split of the branch states.
void testViscoelasticTangentSecondIncrement()
{
  Tensor33d F1 = Spatial3D::I;
  F1( 0, 0 ) += 0.0004;
  F1( 1, 1 ) += 0.0002;
  F1( 0, 1 ) += 0.0002;
  F1( 1, 0 ) += 0.0002;

  Tensor33d F2 = Spatial3D::I;
  F2( 0, 0 ) += 0.0009;
  F2( 1, 1 ) += 0.0005;
  F2( 2, 2 ) -= 0.0002;
  F2( 0, 1 ) += 0.0006;
  F2( 1, 0 ) += 0.0006;

  const double dT = 0.5;

  // increment 1 -- commits a non-zero Maxwell history
  std::vector< double > svBase;
  step( svBase, F1, 0.0, dT, propsVisco );

  // Guard: the committed history must actually influence increment 2, otherwise this test is
  // no stronger than V-1.
  std::vector< double >               svCont = svBase, svFresh;
  DufourModel::AlgorithmicModuli< 3 > tangent;
  const auto                          rCont  = step( svCont, F2, 0.0, dT, propsVisco, &tangent );
  const auto                          rFresh = step( svFresh, F2, 0.0, dT, propsVisco );
  if ( std::abs( rCont.tau( 0, 0 ) - rFresh.tau( 0, 0 ) ) < 1e-6 * std::abs( rFresh.tau( 0, 0 ) ) )
    throw std::runtime_error( "DufourModel V-2: the Maxwell history does not affect the second "
                              "increment - the recursive tangent check would be vacuous" );

  const Tensor3333d fd = finiteDifferenceTangent( svBase, F2, dT, propsVisco, 1e-6 );

  throwExceptionOnFailure( checkIfEqual( tangent.dTau_dF, fd, 1e-3 ),
                           "DufourModel V-2: analytic dTau/dF on the second increment (non-zero "
                           "Maxwell history) must match central finite differences in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// V-3: the model must reduce EXACTLY to the current one when the chain is switched off, both by
// omitting the entries altogether and by declaring branches with zero relative moduli.
void testViscoelasticReducesToBaseModel()
{
  Tensor33d F = Spatial3D::I;
  F( 0, 0 ) += 0.0006;
  F( 1, 1 ) -= 0.0002;
  F( 0, 1 ) += 0.0004;
  F( 1, 0 ) += 0.0004;
  const double dT = 0.5;

  std::vector< double > svRef;
  const auto            ref = step( svRef, F, 0.0, dT, propsBase );

  // (a) branches declared but with G_i = K_i = 0
  const std::vector< double > propsZeroGamma = withMaxwell( { 0.0, 0.0, 0.5, 0.0, 0.0, 5.0 } );
  std::vector< double >       svZero;
  const auto                  rZero = step( svZero, F, 0.0, dT, propsZeroGamma );

  // (b) nMaxwell = 0
  const std::vector< double > propsNoBranch = withMaxwell( {} );
  std::vector< double >       svNone;
  const auto                  rNone = step( svNone, F, 0.0, dT, propsNoBranch );

  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ ) {
      throwExceptionOnFailure( checkIfEqual( rZero.tau( i, j ), ref.tau( i, j ), 1e-12 ),
                               "DufourModel V-3: G_i = K_i = 0 must reproduce the base model in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
      throwExceptionOnFailure( checkIfEqual( rNone.tau( i, j ), ref.tau( i, j ), 1e-12 ),
                               "DufourModel V-3: nMaxwell = 0 must reproduce the base model in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
    }
}

// V-4: physical limits of the relaxation. With one branch of weight g applied equally to the
// deviatoric and volumetric parts, a step taken in dt << tau must retain essentially the full
// instantaneous stiffness, while dt >> tau must relax onto the equilibrium fraction (1 - g).
void testViscoelasticRelaxationLimits()
{
  Tensor33d F = Spatial3D::I;
  F( 0, 0 ) += 0.0005;
  F( 1, 1 ) -= 0.0002;

  const double                g   = 0.4;
  const double                tau = 1.0;
  const std::vector< double > card = withMaxwell( { g, g, tau } );

  std::vector< double > svRef;
  const auto            ref = step( svRef, F, 0.0, 1e-4, propsBase );

  std::vector< double > svFast, svSlow;
  const auto            fast = step( svFast, F, 0.0, 1e-4 * tau, card ); // dt << tau
  const auto            slow = step( svSlow, F, 0.0, 1e4 * tau, card );  // dt >> tau

  throwExceptionOnFailure( checkIfEqual( fast.tau( 0, 0 ), ref.tau( 0, 0 ), 1e-3 ),
                           "DufourModel V-4: dt << tau must retain the instantaneous stiffness in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( slow.tau( 0, 0 ), ( 1.0 - g ) * ref.tau( 0, 0 ), 1e-3 ),
                           "DufourModel V-4: dt >> tau must relax onto the equilibrium fraction in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// V-5: pin the persistent state layout. The viscoelastic block is sized for nMaxwellMax = 7
// branches unconditionally, so that the layout stays static and independent of the card:
//   base    Fp(9) + alphaP + omega + damageDriver + alphaPBar + alphaD + chiF          = 15
//   visco   PK2Ref(9) + veDev(7 x 9) + veVol(7 x 1)                                    = 79
void testStateVarLayout()
{
  DufourModel mat( propsVisco.data(), static_cast< int >( propsVisco.size() ), elLabel );
  throwExceptionOnFailure( checkIfEqual( static_cast< double >( mat.getNumberOfRequiredStateVars() ), 94.0, 1e-12 ),
                           "DufourModel V-5: unexpected number of required state vars in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testUndeformedResponse,
    testStressSymmetry,
    testPureRotation,
    testAlgorithmicTangent,
    testAlgorithmicTangentStressDrivenDamage,
    testViscoelasticTangentFirstIncrement,
    testViscoelasticTangentSecondIncrement,
    testViscoelasticReducesToBaseModel,
    testViscoelasticRelaxationLimits,
    testStateVarLayout,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
