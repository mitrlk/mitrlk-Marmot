#include "Marmot/DufourModel.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotStressMeasures.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  DufourModel::DufourModel( const double* materialProperties, int nMaterialProperties, int materialLabel )
    : MarmotMaterialGradientEnhancedFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
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
      density( nMaterialProperties > 19 ? materialProperties[19] : 0.0 ),
      // optional SWDFM entries 21-24; cSW absent or 0 -> unscaled (original) damage evolution.
      // bSW defaults to a large value (compression term switched off) and kSW to 0 (no Lode
      // dependence) so that only cSW has to be supplied to activate the driver.
      cSW( nMaterialProperties > 20 ? materialProperties[20] : 0.0 ),
      bSW( nMaterialProperties > 21 ? materialProperties[21] : 1e6 ),
      kSW( nMaterialProperties > 22 ? materialProperties[22] : 0.0 ),
      dF( nMaterialProperties > 23 ? materialProperties[23] : 1.0 ),
      volDriver( nMaterialProperties > 24 ? materialProperties[24] : 0.0 ),
      Xt( nMaterialProperties > 25 ? materialProperties[25] : 0.0 ),
      Xc( nMaterialProperties > 26 ? materialProperties[26] : 0.0 ),
      // cubic-exponent Rice-Tracey weight: the three entries AFTER the Prony triplets, in the
      // order (T^2, T^3, T). All absent -> 0 -> the published single-branch exp( 1.3 T ).
      swdfmC2( swdfmExtra( materialProperties, nMaterialProperties, 0 ) ),
      swdfmC3( swdfmExtra( materialProperties, nMaterialProperties, 1 ) ),
      swdfmC1( swdfmExtra( materialProperties, nMaterialProperties, 2 ) ),
      // monotone-by-construction exponent + rate-dependent T-sensitivity; all absent -> the cubic
      swdfmB0( swdfmExtra( materialProperties, nMaterialProperties, 3 ) ),
      swdfmB1( swdfmExtra( materialProperties, nMaterialProperties, 4 ) ),
      swdfmB2( swdfmExtra( materialProperties, nMaterialProperties, 5 ) ),
      swdfmKdotRef( swdfmExtra( materialProperties, nMaterialProperties, 6 ) ),
      swdfmS( swdfmExtra( materialProperties, nMaterialProperties, 7 ) ),
      // saturation value of the SOFTENING variable (Nguyen's Ds_inf); 0/absent -> 1.0
      dsInf( swdfmExtra( materialProperties, nMaterialProperties, 8 ) ),
      // quartic term of the exponent, and the triaxiality above which g is held constant
      swdfmC4( swdfmExtra( materialProperties, nMaterialProperties, 9 ) ),
      swdfmTCap( swdfmExtra( materialProperties, nMaterialProperties, 10 ) ),
      // quadratic acceleration of the softening tail; 0/absent -> the plain exponential
      swdfmBeta( swdfmExtra( materialProperties, nMaterialProperties, 11 ) ),
      // localizing gradient damage: floor R of the interaction function, and its steepness eta
      swdfmLocR( swdfmExtra( materialProperties, nMaterialProperties, 12 ) ),
      swdfmLocEta( swdfmExtra( materialProperties, nMaterialProperties, 13 ) ),
      // optional generalized-Maxwell entries from 27 on; absent or nMaxwell = 0 reproduces the
      // purely hyperelastic-viscoplastic model exactly.
      nMaxwell( nMaterialProperties > 27 ? static_cast< int >( materialProperties[27] ) : 0 ),
      maxwellDev( makeMaxwellProperties( materialProperties, nMaterialProperties, 0 ) ),
      maxwellVol( makeMaxwellProperties( materialProperties, nMaterialProperties, 1 ) )
  {
    if ( maxwellDev.sumGamma >= 1.0 || maxwellVol.sumGamma >= 1.0 )
      throw std::invalid_argument( "DufourModel: the sum of the Maxwell relative moduli must stay below 1, "
                                   "otherwise the equilibrium branch has non-positive stiffness" );

    // LOCALIZING-INTERACTION GATE. R is a FRACTION of ld^2, so it must lie in (0,1]; R > 1 would
    // GROW the interaction with damage, and R < 0 is meaningless. R exactly 0 collapses the operator
    // completely and destroys regularisation, which is the standard critique of localizing models --
    // refuse it and require a positive floor.
    if ( swdfmLocR < 0.0 || swdfmLocR > 1.0 )
      throw std::invalid_argument( "DufourModel: the localizing interaction floor R must lie in "
                                   "[0,1] ( 0 or absent = constant interactions, i.e. the "
                                   "conventional gradient model ); R > 1 would make interactions "
                                   "GROW with damage" );
    if ( swdfmLocEta < 0.0 )
      throw std::invalid_argument( "DufourModel: the localizing interaction steepness eta must be "
                                   "positive ( 0 or absent -> default 5 )" );

    // SOFTENING-ACCELERATION GATE. beta only enters through x = kappa_bar/epsF, so it is
    // dimensionless and a negative value would make omega_s NON-MONOTONE (damage healing) and
    // eventually drive it below zero. Refuse it.
    if ( swdfmBeta < 0.0 )
      throw std::invalid_argument( "DufourModel: swdfmBeta must be >= 0 ( 0 or absent = the plain "
                                   "exponential softening tail ); a negative value makes omega_s "
                                   "non-monotone, i.e. damage would HEAL" );

    // LEGACY-DECK GATE. The two entries after the Prony triplets used to be the two-branch
    // ( swdfmT0, swdfmExp2 ) pair; they are now ( swdfmC2, swdfmC3 ) of the cubic exponent. A deck
    // written for the old law carries e.g. ( 0.88, 4.10 ) and NO third entry, which the new law
    // would read as a perfectly valid but completely different g(T) -- a silent misinterpretation
    // of exactly the kind that has cost this project results before. Refuse it instead: any deck
    // that shapes g at all must state the linear coefficient explicitly.
    if ( ( swdfmC2 != 0.0 || swdfmC3 != 0.0 ) && swdfmC1 == 0.0 )
      throw std::invalid_argument(
        "DufourModel: the card shapes the Rice-Tracey weight (entries [28+3 nMaxwell] and "
        "[29+3 nMaxwell] are nonzero) but gives no linear coefficient at [30+3 nMaxwell]. This is "
        "the signature of a deck written for the SUPERSEDED two-branch ( swdfmT0, swdfmExp2 ) law, "
        "which would be silently re-read as the cubic ( swdfmC2, swdfmC3 ). Regenerate the deck "
        "with the cubic exponent ( c1, c2, c3 ) = ( 3.75, -5.75, 3.50 ), i.e. entries "
        "( c2, c3, c1 ) = ( -5.75, 3.50, 3.75 ) -- see DufourModel.h and HANDOFF S8.4." );

    // A rate exponent without a reference rate is silently inert -- w would never be applied -- so
    // it is far more likely to be a card written wrong than a deliberate choice.
    if ( swdfmS != 0.0 && swdfmKdotRef <= 0.0 )
      throw std::invalid_argument( "DufourModel: a SWDFM rate exponent was given at "
                                   "[35 + 3 nMaxwell] but the reference rate at [34 + 3 nMaxwell] "
                                   "is zero or absent, which switches the rate term OFF. Set the "
                                   "reference rate (calibrated: 20 /s) or clear the exponent." );

    // Both exponent forms at once is ambiguous: swdfmMonotoneForm() silently wins and the cubic
    // entries are ignored, which reads as a card that does something it does not.
    if ( swdfmMonotoneForm() && ( swdfmC1 != 0.0 || swdfmC2 != 0.0 || swdfmC3 != 0.0 ) )
      throw std::invalid_argument( "DufourModel: the card carries BOTH the cubic exponent "
                                   "( c1, c2, c3 ) and the monotone quintic ( b0, b1, b2 ). They are "
                                   "alternatives, not a sum -- clear one of the two." );

    // A quartic exponent WITHOUT a cap is the butt-joint trap: g(1.5) = 9e6 on the calibrated
    // coefficients, so any geometry reaching a triaxiality above the calibration range fails at
    // first load. Refuse it rather than let it happen silently.
    if ( swdfmC4 != 0.0 && swdfmTCap <= 0.0 )
      throw std::invalid_argument( "DufourModel: a quartic exponent coefficient was given at "
                                   "[37 + 3 nMaxwell] but no extrapolation cap at [38 + 3 nMaxwell]. "
                                   "g is calibrated only up to T = 1.087 and the quartic explodes "
                                   "beyond it (g = 493 at T = 1.3, 9e6 at 1.5). Set the cap "
                                   "(calibrated: 1.10)." );

    // Ds_inf above 1 would let the softening variable exceed full degradation on its own, which is
    // exactly what the saturation is there to prevent.
    if ( dsInf > 1.0 )
      throw std::invalid_argument( "DufourModel: Ds_inf (entry [36 + 3 nMaxwell]) must lie in (0, 1]. "
                                   "It is the SATURATION value of the softening variable; above 1 the "
                                   "bulk softening could fail the material by itself, which Nguyen's "
                                   "formulation forbids. Use 0 or omit for the unsaturated legacy "
                                   "behaviour." );
  }

  void DufourModel::computeStress( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    tangents,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement )
  {

    auto&           Fp = stateVars->Fp;
    const Tensor33d FpOld( Fp );
    double&         alphaP       = stateVars->alphaP;
    const double    alphaPOld    = alphaP;
    double&         omega        = stateVars->omega;
    const double    omegaOld     = omega;
    double&         driver       = stateVars->damageDriver;
    const double    driverOld    = driver;
    double&         alphaPBar    = stateVars->alphaPBar;
    const double    alphaPBarOld = alphaPBar;
    double&         alphaD       = stateVars->alphaD;
    const double    alphaDOld    = alphaD;
    double&         chiF         = stateVars->chiF;
    const double    chiFOld      = chiF;

    // LOCALIZING GRADIENT DAMAGE (Poh & Sun 2017): the interaction length COLLAPSES in material
    // that has already started to fail, so a forming crack stops transferring energy into its
    // neighbours. Gated on the STORED driver (lagged one increment) so l is constant within the
    // increment and the element tangent stays exact. R = 0 / absent -> g == 1 -> unchanged.
    response.nonlocalradius = ld * std::sqrt( interactionG( omegaFOfDriver( driverOld ) ) );
    double alphaP_nonlocal  = deformation.A;

    // Published to the Maxwell update, which is reached through computeMandelStress from inside
    // the return-map Newton iteration and so cannot be handed the TimeIncrement directly.
    dTCurrent = timeIncrement.dT;

    using namespace Marmot;
    using namespace Fastor;
    using namespace Eigen;
    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    Tensor33d FeTrial = deformation.F % Fastor::inverse( FpOld );
    double    betaP, dBetaP_dAlphaP;
    std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaPOld );

    Tensor33d dFp;
    dFp.eye();
    Tensor33d Fe = FeTrial;

    const double J = determinant( deformation.F );

    if ( isYielding( FeTrial, betaP, omegaOld, J ) ) {

      size_t counter = 0;

      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      VectorXd X( 11 );
      X.segment( 0, 9 ) = mV9d( FeTrial.data() );
      X( 9 )            = alphaPOld;
      Tensor33d   mandelTrial;
      Tensor3333d dMandelTrial_dFe;
      std::tie( mandelTrial, dMandelTrial_dFe ) = computeMandelStress( FeTrial );
      const double hTrial                       = std::get< 5 >( yieldFunctionFromStress( mandelTrial, betaP, J ) );
      X( 10 )                                   = timeIncrement.dT * eta_VP / std::max( hTrial, 1e-8 );
      VectorXd        dX                        = VectorXd::Zero( 11 );
      VectorXd        R                         = VectorXd::Zero( 11 );
      Eigen::MatrixXd dR_dX( 11, 11 );

      std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial, alphaPOld, timeIncrement.dT, J );

      while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {

        if ( counter > 15 ) {
          /*  std::cout << "R:" << R.transpose() << '\n'
                      << "R_norm" << R.norm() << '\n'
                      << "X:" << X.transpose() << '\n'
                      << "dX norm" << dX.norm() << '\n'
                      << "dR_dX" << '\n'
                      << dR_dX << std::endl;*/
          throw std::runtime_error( "inner newton not converged" );
        }

        dX = -dR_dX.colPivHouseholderQr().solve( R );
        X += dX;
        std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial, alphaPOld, timeIncrement.dT, J );
        counter += 1;
      }
      // update plastic deformation increment
      Fe              = X.segment( 0, 9 ).data();
      dFp             = Fastor::inverse( Fe ) % FeTrial;
      alphaP          = X( 9 );
      response.L      = volDriver != 0.0 ? alphaD : alphaP;
      Tensor33d FpNew = dFp % Fp;
      memcpy( Fp.data(), FpNew.data(), 9 * sizeof( double ) );

      using namespace ContinuumMechanics;
      double      psi_, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal;
      Tensor33d   dOmega_dTau;
      Tensor33d   Ce, dPsi_dCe, tau_eff;
      Tensor3333d dCe_dFe, d2Psi_dCedCe, dTau_dPK2_eff, dTau_dFe_partial_eff;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      // compute Kirchhoff stress FIRST: the damage evolution rate is scaled by the stress
      // triaxiality, so the stress state has to be known before omega can be evaluated.
      // Viscoelastic relaxation is applied here as well; commit = false, the branch states are
      // advanced once at the very end of this routine.
      auto [PK2_eff, dPK2_dCe_eff] = applyViscoelasticity( Tensor33d( 2.0 * dPsi_dCe ),
                                                           Tensor3333d( 2.0 * d2Psi_dCedCe ),
                                                           false );

      std::tie( tau_eff,
                dTau_dPK2_eff,
                dTau_dFe_partial_eff ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2_eff, Fe );

      // compute damage variable
      std::tie( omega,
                dOmega_dAlphaP_local,
                dOmega_dAlphaP_nonlocal,
                driver,
                alphaPBar,
                chiF,
                dOmega_dTau ) = computeOmega( volDriver != 0.0 ? alphaD : alphaP,
                                              alphaP_nonlocal,
                                              tau_eff,
                                              driverOld,
                                              alphaPBarOld,
                                              chiFOld );

      response.tau                  = tau_eff * ( 1.0 - omega );
      response.rho                  = density;
      response.elasticEnergyDensity = psi_;

      // compute tangent operator
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

      MatrixXd dYdDeformation              = MatrixXd::Zero( 11, 11 );
      dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d( einsum< IK, JL, to_IJKL >( Spatial3D::I,
                                                                                           transpose( Fastor::inverse(
                                                                                             FpOld ) ) ) )
                                                     .data() )
                                               .transpose();

      const int       idxF   = 10;
      const Tensor33d Finv_T = transpose( Fastor::inverse( deformation.F ) );

      double betaP_c, dbeta_c;
      std::tie( betaP_c, dbeta_c ) = computeBetaP( alphaP );
      Tensor33d   mandel_c;
      Tensor3333d dmandel_c;
      std::tie( mandel_c, dmandel_c ) = computeMandelStress( Fe );
      const double f_c                = std::get< 0 >( yieldFunctionFromStress( mandel_c, betaP_c, J ) );
      const double betaP_cap_c        = std::max( betaP_c, 1e-12 );
      const double ratio_conv         = Math::macauly( f_c + betaP_cap_c ) / betaP_cap_c;

      dYdDeformation.block< 1, 9 >( idxF, 0 ) = ratio_conv * mV9d( Finv_T.data() ).transpose();

      MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

      Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

      Tensor3333d dPK2_dFe_eff = einsum< ijKL, KLMN >( dPK2_dCe_eff, dCe_dFe );
      Tensor3333d dPK2_dF_eff  = einsum< ijKL, KLMN >( dPK2_dFe_eff, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      Tensor3333d dTau_dF_eff = einsum< IJKL, KLMN >( dTau_dPK2_eff, dPK2_dF_eff ) +
                                einsum< ijKL, KLMN >( dTau_dFe_partial_eff, dFe_dF );

      Tensor33d dAlphaP_local_dF = Tensor33d( Vector9d( dXdDeformation.block< 1, 9 >( 9, 0 ).transpose() ).data() );

      // dilatant plastic volume: ln Jp = ln det F - ln det Fe  (since Jp = J / Je), so
      // d(ln Jp)/dF = F^-T - Fe^-T : dFe/dF, exact from the same sensitivity solve.
      const double lnJpOld = std::log( std::max( Fastor::determinant( FpOld ), 1e-12 ) );
      const double lnJpNew = std::log( std::max( Fastor::determinant( deformation.F ), 1e-12 ) ) -
                             std::log( std::max( Fastor::determinant( Fe ), 1e-12 ) );
      const double dAlphaD = std::max( lnJpNew - lnJpOld, 0.0 ); // clamp: damage is irreversible
      alphaD               = alphaDOld + dAlphaD;
      Tensor33d dAlphaD_dF = Tensor33d( 0.0 );
      if ( dAlphaD > 0.0 ) {
        const Tensor33d FeinvT = Fastor::transpose( Fastor::inverse( Fe ) );
        // dAlphaD_dF(K,L) = Finv_T(K,L) - FeinvT(i,j) * dFe_dF(i,j,K,L)   (explicit, unambiguous)
        for ( int K = 0; K < 3; K++ ) {
          for ( int L = 0; L < 3; L++ ) {
            double acc = Finv_T( K, L );
            for ( int i = 0; i < 3; i++ )
              for ( int j = 0; j < 3; j++ )
                acc -= FeinvT( i, j ) * dFe_dF( i, j, K, L );
            dAlphaD_dF( K, L ) = acc;
          }
        }
      }

      const Tensor33d& dLocal_dF = volDriver != 0.0 ? dAlphaD_dF : dAlphaP_local_dF;

      // stress-driven part of omega: dTau/dF gains  -tau_eff (x) ( dOmega/dTau : dTau_eff/dF )
      Tensor33d dOmegaStress_dF( 0.0 );
      for ( int K = 0; K < 3; K++ )
        for ( int L = 0; L < 3; L++ ) {
          double acc = 0.0;
          for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 3; j++ )
              acc += dOmega_dTau( i, j ) * dTau_dF_eff( i, j, K, L );
          dOmegaStress_dF( K, L ) = acc;
        }

      tangents.dTau_dF = ( 1 - omega ) * dTau_dF_eff - dOmega_dAlphaP_local * Fastor::outer( tau_eff, dLocal_dF ) -
                         Fastor::outer( tau_eff, dOmegaStress_dF );
      tangents.dTau_dA = -tau_eff * dOmega_dAlphaP_nonlocal;
      tangents.dL_dF   = dLocal_dF;
    }
    else {
      using namespace Marmot::ContinuumMechanics;
      double      psi_, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal;
      Tensor33d   dOmega_dTau;
      Tensor33d   Ce, dPsi_dCe, tau_eff;
      Tensor3333d dCe_dFe, d2Psi_dCedCe, dTau_dPK2_eff, dTau_dFe_partial_eff;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      // compute Kirchhoff stress FIRST (needed by the triaxiality-scaled damage evolution)
      auto [PK2_eff, dPK2_dCe_eff] = applyViscoelasticity( Tensor33d( 2.0 * dPsi_dCe ),
                                                           Tensor3333d( 2.0 * d2Psi_dCedCe ),
                                                           false );

      std::tie( tau_eff,
                dTau_dPK2_eff,
                dTau_dFe_partial_eff ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2_eff, Fe );

      std::tie( omega,
                dOmega_dAlphaP_local,
                dOmega_dAlphaP_nonlocal,
                driver,
                alphaPBar,
                chiF,
                dOmega_dTau ) = computeOmega( volDriver != 0.0 ? alphaDOld : alphaPOld,
                                              alphaP_nonlocal,
                                              tau_eff,
                                              driverOld,
                                              alphaPBarOld,
                                              chiFOld );

      response.tau                  = tau_eff * ( 1.0 - omega );
      response.rho                  = density;
      response.elasticEnergyDensity = psi_;
      response.L                    = volDriver != 0.0 ? alphaD : alphaP;

      // compute tangent operator
      Tensor3333d dPK2_dFe    = einsum< ijKL, KLMN >( dPK2_dCe_eff, dCe_dFe );
      Tensor3333d dFe_dF      = einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Fastor::inverse( FpOld ) ) );
      Tensor3333d dPK2_dF_eff = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      Tensor3333d dTau_dF_eff = einsum< IJKL, KLMN >( dTau_dPK2_eff, dPK2_dF_eff ) +
                                einsum< ijKL, KLMN >( dTau_dFe_partial_eff, dFe_dF );

      Tensor33d dAlphaP_local_dF = Tensor33d( 0.0 );

      Tensor33d dOmegaStress_dF( 0.0 );
      for ( int K = 0; K < 3; K++ )
        for ( int L = 0; L < 3; L++ ) {
          double acc = 0.0;
          for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 3; j++ )
              acc += dOmega_dTau( i, j ) * dTau_dF_eff( i, j, K, L );
          dOmegaStress_dF( K, L ) = acc;
        }

      tangents.dTau_dF = ( 1 - omega ) * dTau_dF_eff - Fastor::outer( tau_eff, dOmegaStress_dF );
      tangents.dTau_dA = -tau_eff * dOmega_dAlphaP_nonlocal;
      tangents.dL_dF   = dAlphaP_local_dF;
    }

    // ---- advance the Maxwell branch states EXACTLY ONCE per increment, with the converged
    // elastic deformation. Everything above ran with commit = false, so that the return-map
    // Newton iteration (and the post-convergence re-evaluations of the Mandel stress) could
    // evaluate the relaxed stress as often as needed without corrupting the history.
    if ( nMaxwell > 0 ) {
      using namespace Marmot::ContinuumMechanics;
      Tensor33d   CeFinal;
      Tensor3333d dCeFinal_dFe;
      std::tie( CeFinal, dCeFinal_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      double      psiFinal;
      Tensor33d   dPsiFinal_dCe;
      Tensor3333d d2PsiFinal_dCedCe;
      std::tie( psiFinal,
                dPsiFinal_dCe,
                d2PsiFinal_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( CeFinal, K, G );

      applyViscoelasticity( Tensor33d( 2.0 * dPsiFinal_dCe ), Tensor3333d( 2.0 * d2PsiFinal_dCedCe ), true );
    }
  }

  StateView DufourModel::getStateView( const std::string& stateName )
  {
    return stateVars->getStateView( stateName );
  }

  void DufourModel::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": Not sufficient "
                                                   "stateVars!" );

    this->stateVars = std::make_unique< DufourModelStateVarManager >( stateVars_ );
  }

  void DufourModel::initializeYourself()
  {
    stateVars->Fp.eye();
    // explicit: the accumulated driver and the previous weighted alphaP must start at zero
    stateVars->damageDriver = 0.0;
    stateVars->alphaPBar    = 0.0;
    stateVars->alphaD       = 0.0;
    stateVars->chiF         = 0.0;
    // viscoelasticity: unstressed reference and quiescent Maxwell branches
    stateVars->PK2Ref.zeros();
    for ( int i = 0; i < nMaxwellMax * 9; i++ )
      stateVars->veDev[i] = 0.0;
    for ( int i = 0; i < nMaxwellMax; i++ )
      stateVars->veVol[i] = 0.0;
  }
} // namespace Marmot::Materials
