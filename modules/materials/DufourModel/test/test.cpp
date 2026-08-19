#include "Marmot/DufourModel.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTesting.h"
#include <array>
#include <cstdio>
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
// G-8: the tangent in the regime that actually kills simulations -- damage ACTIVE (D > 1) and
// the RATE factor on. The three existing finite-difference tangent checks (V-1, V-2, and the
// elastic one) all run with g's rate term effectively inert, so the d(logG)/d(kdot) contribution
// added to dD/dAlphaPBar has never been verified against finite differences. The coarse-mesh SLJ
// runs at s >= 0.13 diverge with "Cannot reduce increment size", which is what an inconsistent
// tangent in exactly this regime looks like.
double damagingTangentError( double sVal, bool rateOn )
{
  // A rate-dependent card with an explicit s well above the sentinel value.
  std::vector< double > card = propsVisco;
  // swdfmExtra tail starts at 28 + 3 nMaxwell, with nMaxwell held at index 27. Derive it, never
  // count backwards from the end -- the tail length varies with how many extras a deck supplies.
  const size_t base = 28 + 3 * static_cast< size_t >( card[27] );
  if ( card.size() < base + 8 )
    card.resize( base + 8, 0.0 );
  card[18] = 0.0;                       // m = 0: drive the driver from the LOCAL alphaP. A material-point test
                                        // supplies no nonlocal field (step() passes A = 0), so with m = 1 the
                                        // weighted alphaP is identically zero and the driver can never move.
  card[20]       = 1.797;               // cSW -- propsBase ships 0.0, i.e. the driver switched OFF
  card[base + 6] = rateOn ? 20.0 : 0.0; // swdfmKdotRef; 0 switches the rate term OFF
  card[base + 7] = rateOn ? sVal : 0.0; // swdfmS

  // Drive it hard enough to push the driver past D = 1, in a few committed increments.
  const double          dT = 0.02;
  std::vector< double > sv;
  Tensor33d             F     = Spatial3D::I;
  double                lastD = 0.0;
  int                   n     = 0;
  // Drive to a FIXED DEFORMATION STATE, identical for every s, so that s is the only thing that
  // varies. Stopping on a D threshold instead would reach it at a different kappa_bar for each s
  // and confound the rate exponent with the state.
  for ( ; n < 400; n++ ) {
    F( 0, 0 ) += 0.004;
    F( 1, 1 ) -= 0.0012;
    F( 2, 2 ) -= 0.0012;
    step( sv, F, 0.0, dT, card );
    // layout: Fp occupies 0..8, then alphaP 9, omega 10, damageDriver 11, alphaPBar 12
    lastD = sv[11];
    if ( sv[12] > 0.60 )
      break;
  }
  if ( !( lastD > 1.0 ) )
    throw std::runtime_error( "DufourModel G-8: could not drive the damage driver past D = 1, so "
                              "the damaging-regime tangent check would be vacuous" );

  // One more increment from that committed, damaging state: analytic vs central differences.
  Tensor33d Fnext = F;
  Fnext( 0, 0 ) += 0.004;
  Fnext( 1, 1 ) -= 0.0012;
  Fnext( 2, 2 ) -= 0.0012;

  std::vector< double >               svCont = sv;
  DufourModel::AlgorithmicModuli< 3 > tangent;
  step( svCont, Fnext, 0.0, dT, card, &tangent );

  const Tensor3333d fd = finiteDifferenceTangent( sv, Fnext, dT, card, 1e-7 );

  double num = 0.0, den = 0.0;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ ) {
          const double d = tangent.dTau_dF( i, j, k, l ) - fd( i, j, k, l );
          num += d * d;
          den += fd( i, j, k, l ) * fd( i, j, k, l );
        }
  return std::sqrt( num / std::max( den, 1e-30 ) );
}

void testDamagingRateTangent()
{
  // The analytic tangent must match central differences with damage ACTIVE, at every rate
  // exponent. The three pre-existing finite-difference checks are elastic/viscoelastic only, so
  // the d(logG)/d(kdot) contribution to dD/dAlphaPBar was previously unverified.
  //
  // Compare at a FIXED deformation state. An earlier version of this test stopped on a D
  // threshold, which reaches that threshold at a different kappa_bar for every s (0.570 down to
  // 0.364) and lands the check near the elastic-plastic transition for large s, where the FINITE
  // DIFFERENCE is inaccurate. That produced a spurious "error grows with s" signal (up to 0.019)
  // and a false report of a tangent defect. Held at one state, the tangent is exact for every s.
  for ( double sVal : { 0.0, 0.0435, 0.075, 0.114, 0.15 } ) {
    const bool   rateOn = sVal > 0.0;
    const double err    = damagingTangentError( sVal, rateOn );
    throwExceptionOnFailure( err < 1e-3,
                             "DufourModel G-8: analytic dTau/dF with damage active" +
                               std::string( rateOn ? " and the rate factor on at s = " + std::to_string( sVal )
                                                   : " and the rate term off" ) +
                               " must match central finite differences in " + std::string( __PRETTY_FUNCTION__ ) );
  }
}

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

// G-10: the COMPRESSION BRANCH. g must be < 1 and FALLING for T < 0, whatever the tension
// coefficients do. Without it the calibrated split form gives g(-0.29) = 16.6 (ledger 28.6).
void testCompressionBranchDecays()
{
  // the calibrated SPLIT form: tension exp(0.165T +20.228T^2 -41.044T^3 +22.020T^4), cap 1.10
  std::vector< double > c = withCubicG( 20.228, -41.044, 0.165 ); // (c2, c3, c1)
  c[20]                   = 1.660;
  c.push_back( 0.0 );
  c.push_back( 0.0 );
  c.push_back( 0.0 );  // b0,b1,b2 off
  c.push_back( 0.0 );
  c.push_back( 0.0 );  // kdotRef, s
  c.push_back( 0.0 );  // dsInf -> 1
  c.push_back( 0.0 );  // c4 RETIRED -- must be zero; the compression branch is E = 1.3 T
                       // and does not depend on the tension coefficients anyway
  c.push_back( 1.10 ); // TCap
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
//   base    Fp(9) + alphaP + omega + damageDriver + alphaPBar                          = 13
//   visco   PK2Ref(9) + veDev(7 x 9) + veVol(7 x 1)                                    = 79
//   triax   the clamped triaxiality the driver used, per QP, appended LAST               =  1
void testStateVarLayout()
{
  DufourModel mat( propsVisco.data(), static_cast< int >( propsVisco.size() ), elLabel );
  throwExceptionOnFailure( checkIfEqual( static_cast< double >( mat.getNumberOfRequiredStateVars() ), 93.0, 1e-12 ),
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
    testCompressionBranchDecays,
    testDamagingRateTangent,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
