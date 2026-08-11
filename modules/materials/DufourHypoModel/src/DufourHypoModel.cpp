#include "Marmot/DufourHypoModel.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotPronySeries.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <Eigen/Dense>
#include <cmath>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  /* =============================================================================================
   * FULLY ANALYTIC implementation: a coupled 7-unknown local Newton with an analytic Jacobian,
   * and the consistent tangent from the implicit function theorem.
   *
   * WHY, on the record. Three earlier local solvers all failed IN THE STRUCTURE while passing the
   * single-element check (which uses ~500 tiny increments and therefore cannot see the problem):
   *   1. scalar Newton on dLambda, flow direction at the trial stress (semi-implicit)
   *        -> global step collapsed to maxInc/8 and kept falling; ~3.6 h per Arcan deck
   *   2. same, but the stress iterated to self-consistency (fully implicit)
   *        -> collapsed to maxInc/8 more slowly; ~1 h per deck, still monotonically degrading
   *   3. logarithmic residual, solved for ln(dLambda)
   *        -> MUCH worse (maxInc/512 .. /16384): in log space Newton jumps e^3 per iteration,
   *           which throws the inner stress fixed-point outside its contraction radius
   * The cause is not the automatic differentiation (that gives the EXACT derivative of whatever
   * algorithm is written, and cost only 2x) and not the local damage (the m=0 finite-strain runs
   * converge in the normal 52 increments). It was the improvised local solver. DufourModel is
   * robust because it solves a COUPLED system with an analytic Jacobian; this file now does the
   * same in the hypoelastic setting.
   *
   * UNKNOWNS  X = [ S (6, effective Cauchy stress, Voigt) , dLambda ]        (alphaP eliminated)
   * RESIDUALS
   *   R_S = S - S_trial + Cel : ( dLambda * m(S) )                     m = dg/dSigma
   *   R_f = Phi(S) / betaP(alphaP) - ( dKappa / (dt kappaDot0) )^n
   *         with alphaP = alphaP_n + dLambda * h(S),  h = sqrt( 2/3 m:m )
   * TANGENT (implicit function theorem, R(X, dEps) = 0):
   *   dR_S/dEps = -Cel,  dR_f/dEps = 0   ->   dX/dEps = J^-1 [ Cel ; 0 ]
   *   dSigma/dEps = (1-omega) * ( dS/dEps )   plus the omega-sensitivity through alphaP
   *
   * VOIGT CONVENTIONS
   *   stress ( s11 s22 s33 s12 s13 s23 ),  strain ( e11 e22 e33 2e12 2e13 2e23 )
   * A scalar differentiated componentwise w.r.t. stress-Voigt comes out STRAIN-LIKE (dJ2/ds12 =
   * 2 s12 is the engineering shear entry), so m is directly usable as a strain increment.
   * The tensor contraction of a strain-like Voigt vector with itself is
   *   A:A = v1^2+v2^2+v3^2 + (v4^2+v5^2+v6^2)/2                        -> the matrix Wh below.
   * ============================================================================================= */

  namespace {
    constexpr double newtonTol  = 1e-10;
    constexpr int    nMaxNewton = 40;
    constexpr double tiny       = 1e-30;

    const Vector6d delta = ( Vector6d() << 1., 1., 1., 0., 0., 0. ).finished();

    /// diag(1,1,1,1/2,1/2,1/2): turns a strain-like Voigt self-product into the tensor contraction
    inline Matrix6d Wh()
    {
      Matrix6d W = Matrix6d::Zero();
      W.diagonal() << 1., 1., 1., 0.5, 0.5, 0.5;
      return W;
    }

    /// dJ2/dSigma (strain-like: deviator with doubled shear entries), and its derivative P.
    inline Vector6d dJ2dS( const Vector6d& s )
    {
      const double p = ( s( 0 ) + s( 1 ) + s( 2 ) ) / 3.0;
      Vector6d     d;
      d << s( 0 ) - p, s( 1 ) - p, s( 2 ) - p, 2. * s( 3 ), 2. * s( 4 ), 2. * s( 5 );
      return d;
    }

    /// P = d(dJ2dS)/dSigma
    inline Matrix6d Pmat()
    {
      Matrix6d P = Matrix6d::Zero();
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
          P( i, j ) = ( i == j ? 2. / 3. : -1. / 3. );
      P( 3, 3 ) = P( 4, 4 ) = P( 5, 5 ) = 2.0;
      return P;
    }

    inline double J2of( const Vector6d& s )
    {
      const double p = ( s( 0 ) + s( 1 ) + s( 2 ) ) / 3.0;
      const double a = s( 0 ) - p, b = s( 1 ) - p, c = s( 2 ) - p;
      return 0.5 * ( a * a + b * b + c * c ) + s( 3 ) * s( 3 ) + s( 4 ) * s( 4 ) + s( 5 ) * s( 5 );
    }

    inline double mac( double x )
    {
      return x > 0.0 ? x : 0.0;
    }
    inline double hvs( double x )
    {
      return x > 0.0 ? 1.0 : 0.0;
    }
  } // namespace

  double DufourHypoModel::getDensity( const double* stateVars ) const
  {
    if ( nMaterialProperties < 20 )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": card needs 20 entries" );
    return materialProperties[19];
  }

  DufourHypoModel::DufourHypoModel( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElastic::MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialNumber ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      ft( materialProperties[2] ),
      fc( materialProperties[3] ),
      Q1( materialProperties[4] ),
      Q2( materialProperties[5] ),
      b1( materialProperties[6] ),
      b2( materialProperties[7] ),
      b3( materialProperties[8] ),
      b4( materialProperties[9] ),
      b5( materialProperties[10] ),
      eta_VP( materialProperties[11] ),
      n( materialProperties[12] ),
      nuP_plus( materialProperties[13] ),
      nuP_minus( materialProperties[14] ),
      epsF( materialProperties[15] ),
      omegaMax( materialProperties[16] ),
      ld( materialProperties[17] ),
      m( materialProperties[18] ),
      // entries 20/21 are OPTIONAL: a 20-entry card gives etaRatio = 1 -> constant eta, exactly.
      // etaRatio <= 0 ALSO means "off", same sentinel convention as nuP0. This matters: the shared
      // VE card carries 0.0 in slot 20 as padding, and reading that literally would send
      // etaInf -> 0, i.e. collapse the yield surface. Same-card compatibility with DufourModel
      // depends on this guard.
      etaRatio( nMaterialProperties > 20 && materialProperties[20] > 0.0 ? materialProperties[20] : 1.0 ),
      etaKappa( nMaterialProperties > 21 ? materialProperties[21] : 1.0 ),
      // entries 22/23 OPTIONAL: nuP0 <= 0 -> constant nuP_plus, bit-identical
      nuP0( nMaterialProperties > 22 ? materialProperties[22] : 0.0 ),
      kappaNu( nMaterialProperties > 23 ? materialProperties[23] : 1.0 ),
      nMaxwell( nMaterialProperties > 27 ? static_cast< int >( materialProperties[27] ) : 0 ),
      E( 9.0 * K * G / ( 3.0 * K + G ) ),
      nu( ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) ) ),
      eta( fc / ft ),
      alphaPlus( 9.0 * ( 1.0 - 2.0 * nuP_plus ) / ( 1.0 + nuP_plus ) / 2.0 ),
      alphaMinus( 9.0 * ( 1.0 - 2.0 * nuP_minus ) / ( 1.0 + nuP_minus ) / 2.0 ),
      eta0( fc / ft ),
      etaInf( etaRatio * fc / ft )
  {
    stateLayout.add( "alphaP", 1 );
    stateLayout.add( "omega", 1 );

    // ---------------------------------------------------------------- Prony series (optional)
    using namespace ContinuumMechanics::Viscoelasticity;
    if ( nMaxwell > 0 ) {
      if ( nMaterialProperties < 28 + 3 * nMaxwell )
        throw std::invalid_argument(
          "DufourHypoModel: incomplete (gammaG, gammaK, tau) triplets on the material card" );

      pronyProps.nPronyTerms = static_cast< size_t >( nMaxwell );
      pronyProps.pronyStiffnesses.resize( 6, 6 * nMaxwell );
      pronyProps.pronyRelaxationTimes.resize( 6, 6 * nMaxwell );

      double sumG = 0.0, sumK = 0.0;
      for ( int i = 0; i < nMaxwell; ++i ) {
        const double gG  = materialProperties[28 + 3 * i];
        const double gK  = materialProperties[29 + 3 * i];
        const double tau = materialProperties[30 + 3 * i];
        if ( !( tau > 0.0 ) )
          throw std::invalid_argument( "DufourHypoModel: Maxwell relaxation time must be > 0" );
        sumG += gG;
        sumK += gK;
        pronyProps.pronyStiffnesses.block< 6, 6 >( 0, 6 * i )     = isotropicC( gK * K, gG * G );
        pronyProps.pronyRelaxationTimes.block< 6, 6 >( 0, 6 * i ) = Matrix6d::Constant( tau );
      }
      if ( sumG >= 1.0 || sumK >= 1.0 )
        throw std::invalid_argument( "DufourHypoModel: Maxwell relative moduli must sum to < 1 (the equilibrium branch "
                                     "carries the remainder)" );
      // the equilibrium branch is the REMAINDER, which is the convention the Marmot helper wants
      pronyProps.ultimateStiffnessMatrix = isotropicC( ( 1.0 - sumK ) * K, ( 1.0 - sumG ) * G );

      // 6 x (6 nMaxwell) doubles of branch state
      stateLayout.add( "prony", 36 * nMaxwell );
    }
    stateLayout.finalize();
  }

  /** Raghava part of eq 2.26 (NO /J: the working stress is already Cauchy), with dPhi/dSigma and
   *  dPhi/dkappa.  With eta held constant (etaRatio == 1) the first two outputs are bit-identical
   *  to the previous two-output version and dPhidKappa is exactly 0.
   *
   *  dPhi/dkappa = dPhi/deta * deta/dkappa, and with A = eta - 1, B = sqrt(A^2 I1^2 + 12 eta J2):
   *      dB/deta   = ( A I1^2 + 6 J2 ) / B
   *      dPhi/deta = ( I1 + dB/deta ) / ( 2 eta )  -  Phi / eta                                  */
  void DufourHypoModel::phiAndGrad( const Vector6d& s,
                                    const double    kappa,
                                    double&         Phi,
                                    Vector6d&       dPhi,
                                    double&         dPhidKappa ) const
  {
    const double I1  = s( 0 ) + s( 1 ) + s( 2 );
    const double J2  = J2of( s );
    const double eta = etaOf( kappa );
    const double A   = eta - 1.0;
    const double B   = std::sqrt( std::max( A * A * I1 * I1 + 12.0 * eta * J2, 1e-15 ) );
    Phi              = ( A * I1 + B ) / ( 2.0 * eta );
    dPhi             = ( A * delta + ( A * A * I1 * delta + 6.0 * eta * dJ2dS( s ) ) / B ) / ( 2.0 * eta );

    const double dEta = dEtaOf( kappa );
    if ( dEta == 0.0 ) {
      dPhidKappa = 0.0;
      return;
    }
    const double dBdEta = ( A * I1 * I1 + 6.0 * J2 ) / B;
    dPhidKappa          = ( ( I1 + dBdEta ) / ( 2.0 * eta ) - Phi / eta ) * dEta;
  }

  /// Non-associated potential eq 2.29 (numerator): m = dg/dSigma, M = dm/dSigma, h, dh/dSigma.
  void DufourHypoModel::flowAndGrad( const Vector6d& s,
                                     const double    aPlus,
                                     Vector6d&       mOut,
                                     Matrix6d&       MOut,
                                     double&         hOut,
                                     Vector6d&       dhOut ) const
  {
    const double          alphaPlus = aPlus; // shadows the member: may be the kappa-dependent value
    static const Matrix6d P         = Pmat();
    static const Matrix6d WH        = Wh();

    const double   I1 = s( 0 ) + s( 1 ) + s( 2 );
    const double   p  = I1 / 3.0;
    const double   J2 = J2of( s );
    const Vector6d dJ = dJ2dS( s );

    const double W = std::max( 3.0 * J2 + alphaPlus * mac( p ) * mac( p ) + alphaMinus * mac( -p ) * mac( -p ), 1e-15 );
    const double g = std::sqrt( W );

    // Q = dW/dSigma
    const double   volT = 2.0 * ( alphaPlus * mac( p ) - alphaMinus * mac( -p ) ) / 3.0;
    const Vector6d Q    = 3.0 * dJ + volT * delta;
    mOut                = Q / ( 2.0 * g );

    // d2W/dSigma2, then dm/dSigma = (d2W)/(2g) - m (x) m / g
    const Matrix6d d2W = 3.0 * P + ( 2.0 * ( alphaPlus * hvs( p ) + alphaMinus * hvs( -p ) ) / 9.0 ) *
                                     ( delta * delta.transpose() );
    MOut = d2W / ( 2.0 * g ) - mOut * mOut.transpose() / g;

    const double C = mOut.dot( WH * mOut ); // tensor contraction m:m
    hOut           = std::sqrt( std::max( 2.0 * C / 3.0, tiny ) );
    dhOut          = ( 2.0 / ( 3.0 * hOut ) ) * ( MOut.transpose() * ( WH * mOut ) );
  }

  void DufourHypoModel::computeStress( state3D&        state,
                                       Matrix6d&       dStress_dStrain,
                                       const Vector6d& dStrain,
                                       const timeInfo& timeInfo ) const
  {
    using namespace ContinuumMechanics::Elasticity;
    using namespace ContinuumMechanics::Viscoelasticity;
    const Matrix6d Cel = Isotropic::stiffnessTensor( E, nu );

    double& alphaP = stateLayout.getAs< double& >( state.stateVars, "alphaP" );
    double& omega  = stateLayout.getAs< double& >( state.stateVars, "omega" );

    const double dt        = std::max( timeInfo.dT, tiny );
    const double oneMinusD = std::max( 1.0 - omega, 1e-8 );

    // ---------------------------------------------------------------- viscoelastic predictor
    // With nMaxwell = 0 this is exactly Cel and Cel*dStrain, so a 20-entry card is bit-identical.
    // With branches present, Cve is the ALGORITHMIC viscoelastic tangent and replaces Cel in the
    // local Newton and in the consistent tangent -- using Cel there would leave the Jacobian
    // inconsistent with the residual and destroy quadratic convergence.
    Matrix6d                       Cve  = Cel;
    Vector6d                       dSve = Cel * dStrain;
    PronySeries::mapStateVarMatrix pronyState( nMaxwell > 0 ? stateLayout.getPtr( state.stateVars, "prony" ) : nullptr,
                                               6,
                                               6 * ( nMaxwell > 0 ? nMaxwell : 0 ) );
    if ( nMaxwell > 0 ) {
      dSve.setZero();
      PronySeries::evaluatePronySeries( pronyProps, dSve, Cve, pronyState, dStrain, dt, false );
    }

    // effective stress: f and F are degree-1 homogeneous, so pulling 1/(1-D) out of eqs 2.26/2.29
    // and applying it to the stress is identical to his printed form (see thesis eq. 2.21).
    const Vector6d Strial = state.stress / oneMinusD + dSve;

    const double alphaPn = alphaP;

    double   Phi, dPhiDummy;
    Vector6d dPhi;
    phiAndGrad( Strial, alphaPn, Phi, dPhi, dPhiDummy );

    const double bp0 = betaP( alphaPn );

    if ( Phi - bp0 <= 0.0 ) { // elastic
      const Vector6d Snew = Strial;
      state.stress        = oneMinusD * Snew;
      dStress_dStrain     = oneMinusD * Cve;
      if ( nMaxwell > 0 )
        PronySeries::updateStateVars( pronyProps, pronyState, dStrain, dt );
      state.elasticEnergyDensity += ( state.stress - 0.5 * oneMinusD * dSve ).dot( dStrain );
      return;
    }

    // ---------------------------------------------------------------- coupled local Newton
    const double aPlusLag = alphaPlusOf( alphaPn );

    Vector6d S       = Strial;
    double   dLambda = std::max( ( Phi - bp0 ) / ( 3.0 * G ), 1e-16 );

    Vector6d mv, dh, dPhiS;
    Matrix6d Mm;
    double   h, PhiS, dPhiSdKappa;

    Matrix< double, 7, 7 > J;
    Matrix< double, 7, 1 > R, dX;
    int                    it = 0;
    for ( ;; ) {
      // nuP_plus(kappa) evaluated at the PREVIOUS increment's converged kappa (lagged),
      // so the local Newton sees a constant and the analytic Jacobian is unchanged.
      flowAndGrad( S, aPlusLag, mv, Mm, h, dh );

      // kappa must be known BEFORE Phi now: the surface SHAPE depends on it
      const double aP = alphaPn + dLambda * h;
      phiAndGrad( S, aP, PhiS, dPhiS, dPhiSdKappa );

      const double bp   = std::max( betaP( aP ), 1e-12 );
      const double bpPr = Q1 * std::exp( -b1 * aP ) * ( 1.0 - b1 * aP ) + Q2 * b2 * std::exp( -b2 * aP ) +
                          3.0 * b3 * aP * aP + 2.0 * b4 * aP + b5;
      const double q      = std::max( dLambda * h / ( dt * std::max( eta_VP, tiny ) ), 1e-300 );
      const double qn     = std::pow( q, n );
      const double dqn_dq = n * qn / q;

      R.head< 6 >() = S - Strial + Cve * ( dLambda * mv );
      R( 6 )        = PhiS / bp - qn;

      if ( R.norm() < newtonTol * std::max( 1.0, Strial.norm() ) )
        break;
      if ( it++ == nMaxNewton )
        throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": local Newton failed" );

      J.setZero();
      J.block< 6, 6 >( 0, 0 ) = Matrix6d::Identity() + dLambda * ( Cve * Mm );
      J.block< 6, 1 >( 0, 6 ) = Cve * mv;
      // d/dkappa ( Phi/betaP ) = dPhi/dkappa / bp  -  Phi bp' / bp^2 .  The sign convention below
      // is "-c1", so the SHAPE term enters c1 with a minus.  This is the ONLY new Jacobian term:
      // everything else in the 7x7 is untouched by eta(kappa).
      const double c1         = PhiS / ( bp * bp ) * bpPr - dPhiSdKappa / bp; // through alphaP
      const double c2         = dqn_dq / ( dt * std::max( eta_VP, tiny ) );   // through q
      J.block< 1, 6 >( 6, 0 ) = ( dPhiS / bp - ( c1 + c2 ) * dLambda * dh ).transpose();
      J( 6, 6 )               = -( c1 + c2 ) * h;

      dX = J.fullPivLu().solve( -R );
      // keep dLambda on the admissible half-line without an ad-hoc halving rule
      double relax = 1.0;
      if ( dLambda + dX( 6 ) <= 0.0 )
        relax = 0.9 * dLambda / std::max( -dX( 6 ), tiny );
      S += relax * dX.head< 6 >();
      dLambda += relax * dX( 6 );
    }

    // ---------------------------------------------------------------- consistent tangent
    // R( X, dEps ) = 0 with dR_S/dEps = -Cve, dR_f/dEps = 0  ->  dX/dEps = J^-1 [ Cve ; 0 ]
    Matrix< double, 7, 6 > rhs          = Matrix< double, 7, 6 >::Zero();
    rhs.block< 6, 6 >( 0, 0 )           = Cve;
    const Matrix< double, 7, 6 > dXdEps = J.fullPivLu().solve( rhs );
    const Matrix6d               dSdEps = dXdEps.block< 6, 6 >( 0, 0 );

    // update the state
    const double dKappa   = dLambda * h;
    alphaP                = alphaPn + dKappa;
    const double omegaNew = omegaOf( alphaP );
    const double oMDnew   = std::max( 1.0 - omegaNew, 1e-8 );

    // domega/dEps through alphaP = alphaPn + dLambda h(S)
    double domega_dalphaP = 0.0;
    if ( 1.0 - std::exp( -alphaP / epsF ) < omegaMax )
      domega_dalphaP = std::exp( -alphaP / epsF ) / epsF;
    const Vector6d dalphaP_dEps = dXdEps.block< 1, 6 >( 6, 0 ).transpose() * h + ( dSdEps.transpose() * dh ) * dLambda;

    state.stress    = oMDnew * S;
    dStress_dStrain = oMDnew * dSdEps - S * ( domega_dalphaP * dalphaP_dEps ).transpose();

    omega = omegaNew;

    const Vector6d dEpsP = dLambda * mv;
    const Vector6d dEpsE = dStrain - dEpsP;
    // the Maxwell elements sit in the ELASTIC spring, so they see dEpsE, not the total increment
    if ( nMaxwell > 0 )
      PronySeries::updateStateVars( pronyProps, pronyState, dEpsE, dt );
    const Vector6d dS = Cve * dEpsE;
    state.elasticEnergyDensity += dEpsE.dot( state.stress - 0.5 * oMDnew * dS );
    state.dissipation += dEpsP.dot( dS );
  }
} // namespace Marmot::Materials
