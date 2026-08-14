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
  const std::vector< double > propsBase = { 1793.9, 828.0, 12.0, 21.6,   4051.0, 15.8, 552.0, 243.0, 0.0,
                                            0.0,    15.4,  1e-6, 0.0435, 0.3,    0.5,  1.75,  0.99,  0.12,
                                            1.0,    1.0,   0.0,  1e6,    0.0,    0.1,  0.0,   0.0,   0.0 };

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

  const double                g    = 0.4;
  const double                tau  = 1.0;
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

// ---------------------------------------------------------------------------------------------
// G-1..G-3: the CUBIC-EXPONENT Rice-Tracey weight of the SWDFM driver (HANDOFF S8.4).
//
//   g(T) = exp( c1 T + c2 T^2 + c3 T^3 ),   ( c1, c2, c3 ) = ( 3.75, -5.75, 3.50 )
//
// The three coefficients sit AFTER the Prony triplets, in the card order ( c2, c3, c1 ). These
// tests exist because the joint re-fit was done in Python (calibration/refit_g.py) and the FE
// runs must integrate the SAME function: the reference values below are the ones printed by that
// script, so a drift between the two shows up here rather than as a mis-scaled FE campaign.
// ---------------------------------------------------------------------------------------------

// Append the ( c2, c3, c1 ) shape entries to a 7-branch Prony card.
std::vector< double > withCubicG( double c2, double c3, double c1 )
{
  std::vector< double > card = withMaxwell( { 0.121756, 0.121756, 100.0,    0.058253, 0.058253, 17.78,    0.054825,
                                              0.054825, 3.1623,   0.038975, 0.038975, 0.5623,   0.030845, 0.030845,
                                              0.1,      0.020811, 0.020811, 0.0178,   0.018788, 0.018788, 0.0032 } );
  card[20]                   = 1.797; // cSW, the centring constant of the same fit
  card.push_back( c2 );
  card.push_back( c3 );
  card.push_back( c1 );
  return card;
}

// G-1: the calibrated cubic reproduces refit_g.py's g(T) at the tabulated triaxialities, and is
// monotone increasing over the whole range the ten runs cover.
void testCubicRiceTraceyWeight()
{
  const std::vector< double > card = withCubicG( -5.75, 3.50, 3.75 );
  DufourModel                 mat( card.data(), static_cast< int >( card.size() ), elLabel );

  // T = 0.0 / 0.2 / 0.4 / 0.6 / 0.8 / 1.0 / 1.2, values from HANDOFF S8.4
  const std::array< double, 7 > T   = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 };
  const std::array< double, 7 > ref = { 1.00, 1.73, 2.23, 2.55, 3.04, 4.48, 9.66 };
  for ( size_t i = 0; i < T.size(); i++ )
    throwExceptionOnFailure( checkIfEqual( std::exp( mat.swdfmLogG( T[i] ) ), ref[i], 5e-3 ),
                             "DufourModel G-1: g(T) does not match the calibrated cubic at T = " +
                               std::to_string( T[i] ) + " in " + std::string( __PRETTY_FUNCTION__ ) );

  // admissibility: a void-growth weight must not decrease with triaxiality
  double prev = mat.swdfmLogG( 0.0 );
  for ( int k = 1; k <= 1300; k++ ) {
    const double cur = mat.swdfmLogG( 1.3 * k / 1300.0 );
    throwExceptionOnFailure( cur > prev,
                             "DufourModel G-1: g(T) is not monotone increasing on T = 0..1.3 in " +
                               std::string( __PRETTY_FUNCTION__ ) );
    prev = cur;
  }
}

// G-2: an absent (or all-zero) shape tail must reproduce the published single-branch exp( 1.3 T )
// bit for bit, so that every deck without the tail is unaffected by this change.
void testCubicReducesToSingleBranch()
{
  const std::vector< double > card = withCubicG( 0.0, 0.0, 0.0 );
  DufourModel                 mat( card.data(), static_cast< int >( card.size() ), elLabel );
  for ( double T = -0.5; T <= 1.3001; T += 0.1 )
    throwExceptionOnFailure( checkIfEqual( mat.swdfmLogG( T ), 1.3 * T, 1e-14 ),
                             "DufourModel G-2: a zero shape tail must reduce to exp( 1.3 T ) in " +
                               std::string( __PRETTY_FUNCTION__ ) );
}

// G-3: a deck written for the SUPERSEDED two-branch law -- ( swdfmT0, swdfmExp2 ) = ( 0.88, 4.10 )
// and no third entry -- must be REFUSED, not silently re-read as the cubic ( c2, c3 ).
void testLegacyTwoBranchCardIsRefused()
{
  std::vector< double > card = withCubicG( -5.75, 3.50, 3.75 );
  card.pop_back();              // drop c1: the legacy card has only two entries
  card[card.size() - 2] = 0.88; // old swdfmT0
  card[card.size() - 1] = 4.10; // old swdfmExp2
  bool threw            = false;
  try {
    DufourModel mat( card.data(), static_cast< int >( card.size() ), elLabel );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "DufourModel G-3: a legacy two-branch card was accepted instead of refused in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// ---------------------------------------------------------------------------------------------
// G-4..G-7: the MONOTONE-BY-CONSTRUCTION exponent and the RATE-DEPENDENT T-sensitivity (S9.10).
//
//   g(T, kdot) = exp( E(T) * (kdot/kdot_ref)^(-s) ),   E' = (b0 + b1 T + b2 T^2)^2
//   ( b0, b1, b2 ) = ( 1.798, -0.702, -3.172 ),  kdot_ref = 20,  s < 0 -> s = n
//
// Reference values come from `calibration/refit_rate.py`, so a drift between the Python fit and
// the C++ the FE actually integrates shows up here and not as a mis-scaled campaign.
// ---------------------------------------------------------------------------------------------
std::vector< double > withMonotoneG( double b0, double b1, double b2, double kdotRef, double s )
{
  std::vector< double > card = withCubicG( 0.0, 0.0, 0.0 ); // Prony tail + three cleared cubic slots
  card[20]                   = 2.006;                       // cSW of the same fit
  card.push_back( b0 );
  card.push_back( b1 );
  card.push_back( b2 );
  card.push_back( kdotRef );
  card.push_back( s );
  return card;
}

// G-4: E(T) matches the expanded quintic, and is monotone for ARBITRARY b -- the whole point of
// writing E' as a perfect square is that no parameter choice can break it.
void testMonotoneExponentIsMonotoneForAnyB()
{
  const double bs[5][3] = { { 1.798, -0.702, -3.172 },
                            { -1.559, -0.514, 4.078 },
                            { 0.0, 3.0, -4.0 },
                            { -5.0, 5.0, 5.0 },
                            { 2.0, -9.0, 1.0 } };
  for ( const auto& b : bs ) {
    const std::vector< double > card = withMonotoneG( b[0], b[1], b[2], 0.0, 0.0 );
    DufourModel                 mat( card.data(), static_cast< int >( card.size() ), elLabel );
    double                      prev = mat.swdfmLogG( -0.35 );
    for ( int k = 1; k <= 1700; k++ ) {
      const double cur = mat.swdfmLogG( -0.35 + 1.70 * k / 1700.0 );
      throwExceptionOnFailure( cur >= prev - 1e-12,
                               "DufourModel G-4: E(T) is not monotone for b = " + std::to_string( b[0] ) + ", " +
                                 std::to_string( b[1] ) + ", " + std::to_string( b[2] ) + " in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
      prev = cur;
    }
  }

  // and the closed form itself, against E(T) = 3.233 T -1.262 T^2 -3.638 T^3 +1.113 T^4 +2.012 T^5
  const std::vector< double > card = withMonotoneG( 1.798, -0.702, -3.172, 0.0, 0.0 );
  DufourModel                 mat( card.data(), static_cast< int >( card.size() ), elLabel );
  for ( double T : { 0.2, 0.6, 1.0, 1.077 } ) {
    const double ref = 3.23280 * T - 1.26220 * T * T - 3.63838 * T * T * T + 1.11337 * T * T * T * T +
                       2.01234 * T * T * T * T * T;
    throwExceptionOnFailure( checkIfEqual( mat.swdfmLogG( T ), ref, 1e-3 ),
                             "DufourModel G-4: E(T) does not match the expanded quintic at T = " + std::to_string( T ) +
                               " in " + std::string( __PRETTY_FUNCTION__ ) );
  }
}

// G-5: the rate factor. w = (kdot/20)^(-n) with the sentinel s < 0 tying s to the card's n, and the
// measured span 1.535 / 1.105 / 1.000 / 0.905 at the SLJ / 1 / 10 / 100 mm/s.
void testRateFactorAndSentinel()
{
  const std::vector< double > card = withMonotoneG( 1.798, -0.702, -3.172, 20.0, -1.0 );
  DufourModel                 mat( card.data(), static_cast< int >( card.size() ), elLabel );
  const double                E1 = mat.swdfmLogG( 0.6, 20.0 ); // at the reference rate, w = 1

  const std::array< double, 4 > kdot = { 1.051e-3, 2.0, 20.0, 200.0 };
  const std::array< double, 4 > w    = { 1.535, 1.105, 1.000, 0.905 };
  for ( size_t i = 0; i < kdot.size(); i++ )
    throwExceptionOnFailure( checkIfEqual( mat.swdfmLogG( 0.6, kdot[i] ) / E1, w[i], 3e-3 ),
                             "DufourModel G-5: rate factor wrong at kdot = " + std::to_string( kdot[i] ) + " in " +
                               std::string( __PRETTY_FUNCTION__ ) );

  // the sentinel must reproduce s given EXPLICITLY as n
  const std::vector< double > cardExplicit = withMonotoneG( 1.798, -0.702, -3.172, 20.0, props[12] );
  DufourModel                 matE( cardExplicit.data(), static_cast< int >( cardExplicit.size() ), elLabel );
  throwExceptionOnFailure( checkIfEqual( mat.swdfmLogG( 0.6, 1.051e-3 ), matE.swdfmLogG( 0.6, 1.051e-3 ), 1e-12 ),
                           "DufourModel G-5: the s < 0 sentinel does not equal s = n in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // kdot = 0 (an elastic increment) must stay FINITE -- it is floored, not divided by
  throwExceptionOnFailure( std::isfinite( mat.swdfmLogG( 0.6, 0.0 ) ),
                           "DufourModel G-5: kdot = 0 produced a non-finite exponent in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// G-6: with no reference rate the form must be rate-INDEPENDENT, so that a card carrying only
// ( b0, b1, b2 ) behaves identically at every rate.
void testNoReferenceRateMeansNoRateDependence()
{
  const std::vector< double > card = withMonotoneG( 1.798, -0.702, -3.172, 0.0, 0.0 );
  DufourModel                 mat( card.data(), static_cast< int >( card.size() ), elLabel );
  for ( double kd : { 0.0, 1e-3, 1.0, 1e3 } )
    throwExceptionOnFailure( checkIfEqual( mat.swdfmLogG( 0.8, kd ), mat.swdfmLogG( 0.8, 20.0 ), 1e-14 ),
                             "DufourModel G-6: the exponent moved with kdot although no reference "
                             "rate was given, in " +
                               std::string( __PRETTY_FUNCTION__ ) );
}

// G-7: the two ambiguous cards must be REFUSED -- a rate exponent with no reference rate (inert),
// and both exponent forms at once (the quintic would silently win).
void testAmbiguousShapeCardsAreRefused()
{
  auto refuses = []( const std::vector< double >& c ) {
    try {
      DufourModel mat( c.data(), static_cast< int >( c.size() ), elLabel );
    }
    catch ( const std::invalid_argument& ) {
      return true;
    }
    return false;
  };

  throwExceptionOnFailure( refuses( withMonotoneG( 1.798, -0.702, -3.172, 0.0, 0.05 ) ),
                           "DufourModel G-7: a rate exponent with no reference rate was accepted in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  std::vector< double > both = withMonotoneG( 1.798, -0.702, -3.172, 20.0, -1.0 );
  both[both.size() - 5 - 3]  = -5.75; // resurrect a cubic coefficient alongside the quintic
  throwExceptionOnFailure( refuses( both ),
                           "DufourModel G-7: a card carrying BOTH exponent forms was accepted in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// G-8: Ds_inf, the SATURATION value of the softening variable (Nguyen eq. 79-81). Ds_inf = 1 or
// absent must reproduce the previous law bit for bit; Ds_inf < 1 must BOUND omega_s so that the
// softening variable can never fail the material on its own; Ds_inf > 1 must be refused.
void testDsInfSaturatesTheSofteningVariable()
{
  // omega_s is not directly exposed, so it is read through the UNDAMAGED-state response: at a given
  // deformation the Kirchhoff stress carries the factor (1 - omega), and with cSW = 0 the only
  // damage present is omega_s. Comparing Ds_inf = 1 against Ds_inf = 0.3 therefore isolates it.
  auto card = []( double ds ) {
    std::vector< double > c = withMonotoneG( 1.936, -0.484, -3.249, 20.0, -1.0 );
    c[20]                   = 0.0; // cSW off: omega_f cannot contribute
    c.push_back( ds );
    return c;
  };
  const std::vector< double > c1 = card( 0.0 ), cS = card( 0.3 );
  DufourModel                 m1( c1.data(), static_cast< int >( c1.size() ), elLabel );
  DufourModel                 mS( cS.data(), static_cast< int >( cS.size() ), elLabel );
  throwExceptionOnFailure( checkIfEqual( m1.dsInfEff(), 1.0, 1e-14 ),
                           "DufourModel G-8: an absent Ds_inf must fall back to 1.0 in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( checkIfEqual( mS.dsInfEff(), 0.3, 1e-14 ),
                           "DufourModel G-8: Ds_inf was not read from the card in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // stretch well past yield so omega_s is appreciable, then compare the two
  Tensor33d F( 0.0 );
  F( 0, 0 ) = 1.35;
  F( 1, 1 ) = F( 2, 2 ) = 1.0 / std::sqrt( 1.35 );
  std::vector< double > sv1, svS;
  const double          A  = 0.45; // nonlocal field: well into softening
  const auto            r1 = step( sv1, F, A, 1.0, c1 );
  const auto            rS = step( svS, F, A, 1.0, cS );
  // less degradation -> MORE stress. Bounding omega_s cannot reduce the stress.
  throwExceptionOnFailure( rS.tau( 0, 0 ) > r1.tau( 0, 0 ),
                           "DufourModel G-8: Ds_inf = 0.3 did not raise the stress relative to "
                           "Ds_inf = 1, so omega_s is not being bounded, in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // and Ds_inf > 1 must throw
  bool threw = false;
  try {
    const std::vector< double > cBad = card( 1.4 );
    DufourModel                 mBad( cBad.data(), static_cast< int >( cBad.size() ), elLabel );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "DufourModel G-8: Ds_inf > 1 was accepted instead of refused in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// G-9: the QUARTIC exponent and the EXTRAPOLATION CAP.
//   E(T) = c1 T + c2 T^2 + c3 T^3 + c4 T^4,  evaluated at min(T, TCap)
// Calibrated: (c1,c2,c3,c4) = (-0.07, 18.12, -35.58, 18.85), TCap = 1.10, cSW = 1.754.
// The reference values are the ones the Python fit produced, so a drift shows up here.
void testQuarticExponentAndCap()
{
  // card order after the Prony triplets: c2, c3, c1, b0, b1, b2, kdotRef, s, dsInf, c4, TCap
  auto card = []( double TCap ) {
    std::vector< double > c = withCubicG( 18.12, -35.58, -0.07 ); // (c2, c3, c1)
    c[20]                   = 1.754;                              // cSW of the same fit
    c.push_back( 0.0 );
    c.push_back( 0.0 );
    c.push_back( 0.0 );   // b0,b1,b2 -> quartic form off
    c.push_back( 0.0 );
    c.push_back( 0.0 );   // kdotRef, s
    c.push_back( 0.0 );   // dsInf -> 1
    c.push_back( 18.85 ); // c4
    c.push_back( TCap );
    return c;
  };
  const std::vector< double > c = card( 1.10 );
  DufourModel                 mat( c.data(), static_cast< int >( c.size() ), elLabel );

  // the four triaxialities the four specimens actually sit at
  const std::array< double, 4 > T   = { 0.0, 0.584, 0.806, 1.087 };
  const std::array< double, 4 > ref = { 1.00, 3.47, 2.83, 7.05 };
  for ( size_t i = 0; i < T.size(); i++ )
    throwExceptionOnFailure( checkIfEqual( std::exp( mat.swdfmLogG( T[i] ) ), ref[i], 2e-2 ),
                             "DufourModel G-9: quartic g(T) wrong at T = " + std::to_string( T[i] ) + " in " +
                               std::string( __PRETTY_FUNCTION__ ) );

  // NON-MONOTONE by design: it must DIP between the SLJ and Arcan 45 deg
  throwExceptionOnFailure( mat.swdfmLogG( 0.584 ) > mat.swdfmLogG( 0.806 ),
                           "DufourModel G-9: g must DIP from T = 0.584 to 0.806 -- that dip is the "
                           "whole point of the quartic, see ledger 28.2, in " +
                             std::string( __PRETTY_FUNCTION__ ) );

  // THE CAP. Uncapped, g(1.5) = 9e6 and a butt joint fails at first load.
  const double gCap = std::exp( mat.swdfmLogG( 1.10 ) );
  for ( double Tq : { 1.10, 1.3, 1.5, 3.0 } )
    throwExceptionOnFailure( checkIfEqual( std::exp( mat.swdfmLogG( Tq ) ), gCap, 1e-12 ),
                             "DufourModel G-9: g is not held constant above the cap at T = " + std::to_string( Tq ) +
                               " in " + std::string( __PRETTY_FUNCTION__ ) );
  const std::vector< double > cNo = card( 0.0 );
  // no cap given with a quartic present -> must throw (the butt-joint guard)
  bool threw = false;
  try {
    DufourModel bad( cNo.data(), static_cast< int >( cNo.size() ), elLabel );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "DufourModel G-9: a quartic exponent with NO extrapolation cap was "
                           "accepted; g(1.5) = 9e6 would fail a butt joint at first load, in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// G-11..G-14: the SOFTENING TAIL ACCELERATION (swdfmBeta).
//   omega_s = Ds_inf [ 1 - exp( -x - beta x^2 ) ],  x = kappa_bar / epsF
// Measured basis: Gf = 2.890 N/mm for this card (ledger 29.24), joint needs 0.65 (Dufour's SLJ
// inverse fit). Gf ~ ld * int(1-omega_s) dkappa_bar = ld * epsF for beta = 0.
std::vector< double > withBeta( double beta )
{
  std::vector< double > c = withCubicG( -5.75, 3.50, 3.75 );
  c.push_back( 0.0 );
  c.push_back( 0.0 );
  c.push_back( 0.0 );  // b0,b1,b2 off
  c.push_back( 20.0 );
  c.push_back( -1.0 ); // kdotRef, s sentinel
  c.push_back( 0.0 );  // dsInf -> 1
  c.push_back( 0.0 );
  c.push_back( 0.0 );  // c4, TCap off
  c.push_back( beta ); // swdfmBeta
  return c;
}

// G-11: beta = 0 must reproduce the plain exponential BIT FOR BIT.
void testBetaZeroReproducesExponential()
{
  const std::vector< double > c = withBeta( 0.0 );
  DufourModel                 a( c.data(), static_cast< int >( c.size() ), elLabel );
  const double                epsF = 1.75;
  for ( double k : { 0.0, 0.01, 0.1, 0.5, 1.0, 1.75, 3.0, 10.0 } ) {
    double om, dom;
    a.softening( k, om, dom );
    throwExceptionOnFailure( checkIfEqual( om, 1.0 - std::exp( -k / epsF ), 1e-14 ) &&
                               checkIfEqual( dom, std::exp( -k / epsF ) / epsF, 1e-14 ),
                             "DufourModel G-11: beta = 0 changed the softening law at kappa_bar = " +
                               std::to_string( k ) + " in " + std::string( __PRETTY_FUNCTION__ ) );
  }
}

// G-12: THE PROPERTY THE WHOLE CHOICE RESTS ON -- the initial slope is 1/epsF for EVERY beta, so
// epsF keeps the meaning Dufour's DIC fit gave it. Also: monotone, bounded, and NEVER reaching 1 at
// finite kappa_bar (so no second failure criterion is smuggled in).
void testBetaPreservesInitialSlope()
{
  const double epsF = 1.75;
  for ( double beta : { 0.0, 1.475, 3.859, 11.266, 100.0 } ) {
    const std::vector< double > c = withBeta( beta );
    DufourModel                 mat( c.data(), static_cast< int >( c.size() ), elLabel );
    double                      om, dom;
    mat.softening( 0.0, om, dom );
    throwExceptionOnFailure( checkIfEqual( dom, 1.0 / epsF, 1e-12 ),
                             "DufourModel G-12: the initial softening slope must be 1/epsF for EVERY "
                             "beta -- that is the entire reason this form was chosen. beta = " +
                               std::to_string( beta ) + " gave " + std::to_string( dom ) + " in " +
                               std::string( __PRETTY_FUNCTION__ ) );
    double prev = -1.0;
    for ( int i = 0; i <= 2000; ++i ) {
      const double k = 20.0 * i / 2000.0;
      mat.softening( k, om, dom );
      throwExceptionOnFailure( om >= prev - 1e-14 && om <= 1.0 && dom >= -1e-14,
                               "DufourModel G-12: omega_s must be monotone, bounded by 1 and have "
                               "non-negative slope; failed at kappa_bar = " +
                                 std::to_string( k ) + ", beta = " + std::to_string( beta ) + " in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
      prev = om;
    }
    // Strictly BELOW 1 over the strain range that actually OCCURS, so no second failure criterion
    // is smuggled in. The largest kappa_bar measured in any specimen is 0.559 (Arcan 90 deg, at
    // D = 1); 1.0 is that with margin. Beyond kappa_bar ~ 3 the exponential UNDERFLOWS and
    // 1 - exp(...) == 1.0 exactly in double -- a floating-point artefact at strains no element
    // reaches, not a finite failure strain, so it is deliberately not tested there.
    mat.softening( 1.0, om, dom );
    throwExceptionOnFailure( om < 1.0,
                             "DufourModel G-12: omega_s reached 1 at kappa_bar = 1.0, well inside the "
                             "range elements actually visit -- that would add a second failure "
                             "criterion (beta = " +
                               std::to_string( beta ) + ") in " + std::string( __PRETTY_FUNCTION__ ) );
    // and no JUMP anywhere in that range: the form must be smooth (the hard cutoff was rejected
    // partly for imposing a 78% instantaneous stress drop)
    double omPrev, dPrev;
    mat.softening( 0.0, omPrev, dPrev );
    for ( int i = 1; i <= 4000; ++i ) {
      const double k = 1.0 * i / 4000.0;
      mat.softening( k, om, dom );
      throwExceptionOnFailure( om - omPrev < 0.02,
                               "DufourModel G-12: omega_s JUMPS by " + std::to_string( om - omPrev ) +
                                 " in one 2.5e-4 strain step at kappa_bar = " + std::to_string( k ) +
                                 " (beta = " + std::to_string( beta ) + "); the form must be smooth, in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
      omPrev = om;
    }
  }
}

// G-13: the returned derivative must BE the derivative of the returned value.
void testBetaDerivativeIsExact()
{
  for ( double beta : { 0.0, 1.475, 3.859, 11.266 } ) {
    const std::vector< double > c = withBeta( beta );
    DufourModel                 mat( c.data(), static_cast< int >( c.size() ), elLabel );
    for ( double k : { 0.02, 0.1, 0.3, 0.6, 1.0, 1.75, 3.0 } ) {
      const double h = 1e-7;
      double       op, dop, om1, om2, dummy;
      mat.softening( k, op, dop );
      mat.softening( k + h, om2, dummy );
      mat.softening( k - h, om1, dummy );
      const double fd = ( om2 - om1 ) / ( 2.0 * h );
      throwExceptionOnFailure( std::abs( fd - dop ) < 1e-5 * std::max( 1.0, std::abs( dop ) ),
                               "DufourModel G-13: dOmega_s disagrees with the finite difference "
                               "(analytic " +
                                 std::to_string( dop ) + " vs " + std::to_string( fd ) +
                                 ") at kappa_bar = " + std::to_string( k ) + ", beta = " + std::to_string( beta ) +
                                 " in " + std::string( __PRETTY_FUNCTION__ ) );
    }
  }
}

// G-14: the AREA int(1-omega_s) dkappa_bar -- which Gf is proportional to -- must hit the values the
// calibration relies on, and a negative beta must be REFUSED (it would heal damage).
void testBetaShrinksTheArea()
{
  const double epsF = 1.75;
  auto         area = [&]( double beta ) {
    const std::vector< double > c = withBeta( beta );
    DufourModel                 mat( c.data(), static_cast< int >( c.size() ), elLabel );
    const double                K = 80.0;
    const int                   N = 800000;
    double                      s = 0.0, om, dom;
    for ( int i = 0; i < N; ++i ) {
      mat.softening( ( i + 0.5 ) * K / N, om, dom );
      s += ( 1.0 - om ) * K / N;
    }
    return s;
  };
  throwExceptionOnFailure( checkIfEqual( area( 0.0 ), epsF, 1e-3 ),
                           "DufourModel G-14: the beta = 0 area must be epsF = 1.75, got " +
                             std::to_string( area( 0.0 ) ) + " in " + std::string( __PRETTY_FUNCTION__ ) );
  // beta values derived from the measured Gf = 2.890 N/mm: area = 1.75 * Gf_target / 2.890
  const double a11 = area( 11.266 ), a39 = area( 3.859 );
  throwExceptionOnFailure( std::abs( a11 - 1.75 * 0.65 / 2.890 ) < 0.01,
                           "DufourModel G-14: beta = 11.266 must give the area for Gf = 0.65 N/mm "
                           "(0.3936), got " +
                             std::to_string( a11 ) + " in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( std::abs( a39 - 1.75 * 1.00 / 2.890 ) < 0.01,
                           "DufourModel G-14: beta = 3.859 must give the area for Gf = 1.00 N/mm "
                           "(0.6055), got " +
                             std::to_string( a39 ) + " in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( a11 < a39 && a39 < area( 0.0 ),
                           "DufourModel G-14: larger beta must give a SMALLER area in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  bool refused = false;
  try {
    const std::vector< double > bad = withBeta( -1.0 );
    DufourModel                 m2( bad.data(), static_cast< int >( bad.size() ), elLabel );
  }
  catch ( const std::invalid_argument& ) {
    refused = true;
  }
  throwExceptionOnFailure( refused,
                           "DufourModel G-14: a NEGATIVE beta must be refused (damage would "
                           "heal) in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// G-15..G-17: LOCALIZING GRADIENT DAMAGE (Poh & Sun 2017), interaction function g(omega_f).
//   nonlocalradius = ld * sqrt( g ),  g(0) = 1, g(1) = R.
std::vector< double > withLoc( double R, double eta )
{
  std::vector< double > c = withBeta( 0.0 ); // beta off: this must be an independent mechanism
  c.push_back( R );
  c.push_back( eta );
  return c;
}

// G-15: R = 0 (and absent) must give g == 1 for EVERY omega_f -- i.e. the conventional gradient
// model, bit for bit. This is the property that makes the mechanism inert on the Arcan PEAKS, where
// the measured fraction of cells past D = 1 is 0.0% / 0.0% / 1.1%.
void testLocalizingOffIsIdentity()
{
  for ( auto c : { withBeta( 0.0 ), withLoc( 0.0, 0.0 ), withLoc( 0.0, 5.0 ) } ) {
    DufourModel mat( c.data(), static_cast< int >( c.size() ), elLabel );
    for ( double w : { 0.0, 0.01, 0.3, 0.7, 0.999, 1.0 } )
      throwExceptionOnFailure( checkIfEqual( mat.interactionG( w ), 1.0, 1e-15 ),
                               "DufourModel G-15: with R = 0 the interaction must be EXACTLY 1 at "
                               "omega_f = " +
                                 std::to_string( w ) + " (got " + std::to_string( mat.interactionG( w ) ) + ") in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
  }
}

// G-16: the end points and monotonicity. g(0) = 1 EXACTLY is what guarantees undamaged material is
// untouched; g(1) = R EXACTLY is what bounds the collapse.
void testLocalizingEndpointsAndMonotone()
{
  for ( double R : { 0.05, 0.3, 0.5, 1.0 } )
    for ( double eta : { 0.0, 1.0, 5.0, 20.0 } ) {
      std::vector< double > c = withLoc( R, eta );
      DufourModel           mat( c.data(), static_cast< int >( c.size() ), elLabel );
      throwExceptionOnFailure( checkIfEqual( mat.interactionG( 0.0 ), 1.0, 1e-12 ),
                               "DufourModel G-16: g(0) must be exactly 1 (R = " + std::to_string( R ) +
                                 ", eta = " + std::to_string( eta ) + ") in " + std::string( __PRETTY_FUNCTION__ ) );
      throwExceptionOnFailure( checkIfEqual( mat.interactionG( 1.0 ), R, 1e-12 ),
                               "DufourModel G-16: g(1) must be exactly R = " + std::to_string( R ) + ", got " +
                                 std::to_string( mat.interactionG( 1.0 ) ) + " in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
      double prev = 2.0;
      for ( int i = 0; i <= 500; ++i ) {
        const double w = i / 500.0, g = mat.interactionG( w );
        throwExceptionOnFailure( g <= prev + 1e-14 && g >= R - 1e-12 && g <= 1.0 + 1e-12,
                                 "DufourModel G-16: g must decrease monotonically and stay in [R,1]; "
                                 "failed at omega_f = " +
                                   std::to_string( w ) + " (g = " + std::to_string( g ) +
                                   ", R = " + std::to_string( R ) + ") in " + std::string( __PRETTY_FUNCTION__ ) );
        prev = g;
      }
      // clamped outside [0,1] rather than extrapolating
      throwExceptionOnFailure( checkIfEqual( mat.interactionG( -1.0 ), 1.0, 1e-12 ) &&
                                 checkIfEqual( mat.interactionG( 2.0 ), R, 1e-12 ),
                               "DufourModel G-16: omega_f outside [0,1] must be CLAMPED in " +
                                 std::string( __PRETTY_FUNCTION__ ) );
    }
}

// G-17: omegaFOfDriver must agree with the omega_f the damage law itself uses, and the gate must be
// exactly inert below initiation (D <= 1) -- that is what keeps the Arcan peaks untouched.
void testLocalizingGateIsInertBeforeInitiation()
{
  std::vector< double > c = withLoc( 0.3, 5.0 );
  DufourModel           mat( c.data(), static_cast< int >( c.size() ), elLabel );
  // Do NOT hard-code dF: recover it FROM the function and require every D to imply the SAME one.
  // That tests the functional form 1-exp(-(D-1)/dF) without depending on the test card's value,
  // which is how the first version of this test broke (it assumed 0.3; the card carries 0.1).
  const double dF = -( 1.5 - 1.0 ) / std::log( 1.0 - mat.omegaFOfDriver( 1.5 ) );
  for ( double D : { 0.0, 0.5, 0.99, 1.0 } )
    throwExceptionOnFailure( checkIfEqual( mat.omegaFOfDriver( D ), 0.0, 1e-15 ) &&
                               checkIfEqual( mat.interactionG( mat.omegaFOfDriver( D ) ), 1.0, 1e-15 ),
                             "DufourModel G-17: below initiation (D = " + std::to_string( D ) +
                               ") omega_f must be 0 and the interaction EXACTLY 1 in " +
                               std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( dF > 0.0 && dF < 10.0,
                           "DufourModel G-17: the dF recovered from omegaFOfDriver is implausible (" +
                             std::to_string( dF ) + ") in " + std::string( __PRETTY_FUNCTION__ ) );
  for ( double D : { 1.001, 1.1, 1.5, 3.0 } ) {
    const double expect = 1.0 - std::exp( -( D - 1.0 ) / dF );
    throwExceptionOnFailure( checkIfEqual( mat.omegaFOfDriver( D ), expect, 1e-10 ),
                             "DufourModel G-17: omega_f is not 1-exp(-(D-1)/dF) with a SINGLE dF; at "
                             "D = " +
                               std::to_string( D ) + " it gives " + std::to_string( mat.omegaFOfDriver( D ) ) +
                               " but the dF implied at "
                               "D = 1.5 predicts " +
                               std::to_string( expect ) + ", in " + std::string( __PRETTY_FUNCTION__ ) );
    throwExceptionOnFailure( mat.interactionG( mat.omegaFOfDriver( D ) ) < 1.0,
                             "DufourModel G-17: past initiation the interaction must have COLLAPSED "
                             "below 1 at D = " +
                               std::to_string( D ) + " in " + std::string( __PRETTY_FUNCTION__ ) );
  }
  // R outside [0,1] and negative eta must be refused
  int refused = 0;
  for ( auto bad : { withLoc( -0.1, 5.0 ), withLoc( 1.5, 5.0 ), withLoc( 0.3, -1.0 ) } ) {
    try {
      DufourModel m2( bad.data(), static_cast< int >( bad.size() ), elLabel );
    }
    catch ( const std::invalid_argument& ) {
      ++refused;
    }
  }
  throwExceptionOnFailure( refused == 3,
                           "DufourModel G-17: R < 0, R > 1 and eta < 0 must all be refused (only " +
                             std::to_string( refused ) + " of 3 were) in " + std::string( __PRETTY_FUNCTION__ ) );
}

// G-10: the COMPRESSION BRANCH. g must be < 1 and FALLING for T < 0, whatever the tension
// coefficients do. Without it the calibrated split form gives g(-0.29) = 16.6 (ledger 28.6).
void testCompressionBranchDecays()
{
  // the calibrated SPLIT form: tension exp(0.165T +20.228T^2 -41.044T^3 +22.020T^4), cap 1.10
  std::vector< double > c = withCubicG( 20.228, -41.044, 0.165 ); // (c2, c3, c1)
  c[20]                   = 1.660;
  c.push_back( 0.0 );
  c.push_back( 0.0 );
  c.push_back( 0.0 );    // b0,b1,b2 off
  c.push_back( 0.0 );
  c.push_back( 0.0 );    // kdotRef, s
  c.push_back( 0.0 );    // dsInf -> 1
  c.push_back( 22.020 ); // c4
  c.push_back( 1.10 );   // TCap
  DufourModel mat( c.data(), static_cast< int >( c.size() ), elLabel );

  // compression: strictly below 1 and monotonically falling as T decreases
  double prev = 1.0;
  for ( double T : { -0.02, -0.05, -0.10, -0.20, -0.29, -0.35 } ) {
    const double g = std::exp( mat.swdfmLogG( T ) );
    throwExceptionOnFailure( g < 1.0 && g < prev,
                             "DufourModel G-10: g must be < 1 and falling in compression, got " + std::to_string( g ) +
                               " at T = " + std::to_string( T ) + " in " + std::string( __PRETTY_FUNCTION__ ) );
    prev = g;
  }
  // the value the ledger records as the failure mode this branch prevents
  throwExceptionOnFailure( std::exp( mat.swdfmLogG( -0.29 ) ) < 0.8,
                           "DufourModel G-10: g(-0.29) should be ~0.69, not the 16.6 the unsplit "
                           "polynomial gives, in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  // continuity at zero, and the tension INVERSION still present
  throwExceptionOnFailure( checkIfEqual( std::exp( mat.swdfmLogG( 0.0 ) ), 1.0, 1e-12 ),
                           "DufourModel G-10: g(0) must be exactly 1 in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( mat.swdfmLogG( 0.584 ) > mat.swdfmLogG( 0.806 ),
                           "DufourModel G-10: the tension inversion g(0.584) > g(0.806) was lost in " +
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
    testCubicRiceTraceyWeight,
    testCubicReducesToSingleBranch,
    testLegacyTwoBranchCardIsRefused,
    testMonotoneExponentIsMonotoneForAnyB,
    testRateFactorAndSentinel,
    testNoReferenceRateMeansNoRateDependence,
    testAmbiguousShapeCardsAreRefused,
    testDsInfSaturatesTheSofteningVariable,
    testQuarticExponentAndCap,
    testCompressionBranchDecays,
    testBetaZeroReproducesExponential,
    testBetaPreservesInitialSlope,
    testBetaDerivativeIsExact,
    testBetaShrinksTheArea,
    testLocalizingOffIsIdentity,
    testLocalizingEndpointsAndMonotone,
    testLocalizingGateIsInertBeforeInitiation,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
