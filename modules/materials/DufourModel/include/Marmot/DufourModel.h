/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Alexander Dummer alexander.dummer@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#pragma once
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainPlasticity.h"
#include "Marmot/MarmotMaterialGEFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include <string>

namespace Marmot::Materials {

  using namespace Fastor;
  using namespace FastorStandardTensors;
  using namespace FastorIndices;

  class DufourModel : public MarmotMaterialGEFiniteStrain {
  public:
    using MarmotMaterialGEFiniteStrain::MarmotMaterialGEFiniteStrain;

    // elastic constants
    const double K, G;

    // plasticity parameters -- possibly eta as a separate parameter, instead of fc, for pressure sensitivity
    const double ft, fc, Q1, Q2, b1, b2, b3, b4, b5;

    // strain rate
    const double eta_VP, n;

    // viscoplastic flow
    const double nuP_plus, nuP_minus;

    // damage
    const double epsF, omegaMax, ld, m;

    // mass properties;
    const double density;

    DufourModel( const double* materialProperties, int nMaterialProperties, int materialLabel );

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement );

    int getNumberOfRequiredStateVars() { return DufourModelStateVarManager::layout.nRequiredStateVars; }

    double getDensity() { return density; }

    class DufourModelStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "Fp", .length = 9 },
        { .name = "alphaP", .length = 1 },
        { .name = "omega", .length = 1 },
      } );

      Fastor::TensorMap< double, 3, 3 > Fp;
      double&                           alphaP;
      double&                           omega;

      DufourModelStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          Fp( &find( "Fp" ) ),
          alphaP( find( "alphaP" ) ),
          omega( find( "omega" ) ){};
    };
    std::unique_ptr< DufourModelStateVarManager > stateVars;

    void assignStateVars( double* stateVars, int nStateVars );

    StateView getStateView( const std::string& result );

    void initializeYourself();

    // ------------------------------------------------------------
    // Damage functions
    // ------------------------------------------------------------

    std::tuple< double, double, double > computeOmega( const double alphaP_local, const double alphaP_nonlocal )
    {
      double alphaP_weighted = alphaP_nonlocal * m + alphaP_local * ( 1 - m );

      double dAlphaP_weighted_dAlphaP_local    = 1 - m;
      double dAlphaP_weighted_dAlphaP_nonlocal = m;

      if ( alphaP_weighted < 0.0 ) {
        return { 0.0, 0.0, 0.0 };
      }

      const double omega = 1.0 - exp( -alphaP_weighted / epsF );

      double dOmega_dAlphaP_weigthed = 1.0 / epsF * exp( -alphaP_weighted / epsF );
      double dOmega_dAlphaP_local    = dOmega_dAlphaP_weigthed * dAlphaP_weighted_dAlphaP_local;
      double dOmega_dAlphaP_nonlocal = dOmega_dAlphaP_weigthed * dAlphaP_weighted_dAlphaP_nonlocal;

      if ( omega > omegaMax ) {
        return { omegaMax, 0.0, 0.0 };
      }
      else {
        return { omega, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal };
      }
    }

    // ------------------------------------------------------------
    // Viscoplasticity functions
    // ------------------------------------------------------------

    std::tuple< double, Tensor33d, double, Tensor33d, Tensor33d > yieldFunction( const Tensor33d& Fe,
                                                                                 const double     betaP )
    {

      Tensor33d   mandelStress;
      Tensor3333d dMandel_dFe;
      std::tie( mandelStress, dMandel_dFe ) = computeMandelStress( Fe );

      double      f, h, df_dBetaP;
      Tensor33d   df_dMandel, dh_dMandel, df_dFe, dg_dMandel, dh_dFe;
      Tensor3333d d2g_dMandel_dMandel;

      std::tie( f,
                df_dMandel,
                df_dBetaP,
                dg_dMandel,
                d2g_dMandel_dMandel,
                h,
                dh_dMandel ) = yieldFunctionFromStress( mandelStress, betaP );
      dh_dFe                 = einsum< mn, mnij, to_ij >( dh_dMandel, dMandel_dFe );
      df_dFe                 = einsum< mn, mnij, to_ij >( df_dMandel, dMandel_dFe );

      return { f, df_dFe, df_dBetaP, dg_dMandel, dh_dFe };
    }

    std::tuple< double > yieldFunctionNominal( const Tensor33d& Fe, const double betaP, const double omega )
    {

      Tensor33d   mandelStressN;
      Tensor3333d dMandel_dFe;
      std::tie( mandelStressN, dMandel_dFe ) = computeNominalMandelStress( Fe, omega );

      const double eta = fc / ft;
      Tensor33d    dev = deviatoric( mandelStressN );
      const double J2  = 0.5 * Fastor::inner( dev, dev );
      const double I1  = trace( mandelStressN );

      const double A = eta - 1.0;
      const double B = sqrt( std::max( A * A * I1 * I1 + 12.0 * eta * J2, 1e-15 ) );

      const double f = ( A * I1 + B ) / ( 2.0 * eta ) - betaP;

      return { f };
    }

    std::tuple< double, Tensor33d, double, Tensor33d, Tensor3333d, double, Tensor33d > yieldFunctionFromStress(
      const Tensor33d& mandelStress,
      const double     betaP )
    {
      const double eta = fc / ft;
      Tensor33d    dev = deviatoric( mandelStress );
      const double J2  = 0.5 * Fastor::inner( dev, dev );
      const double I1  = trace( mandelStress );
      const double p   = I1 / 3.0;

      const double thetaPlus  = 9 * ( ( 1 - 2 * nuP_plus ) / ( 1 + nuP_plus ) ) / 2.0;
      const double thetaMinus = 9 * ( ( 1 - 2 * nuP_minus ) / ( 1 + nuP_minus ) ) / 2.0;

      const double A = eta - 1.0;
      const double B = sqrt( std::max( A * A * I1 * I1 + 12.0 * eta * J2, 1e-15 ) );

      const double f = ( A * I1 + B ) / ( 2.0 * eta ) - betaP;

      const double phiI1      = A * ( 1.0 + A * I1 / B ) / ( 2.0 * eta );
      const double phiJ2      = 3.0 / B;
      Tensor33d    df_dMandel = phiI1 * Spatial3D::I + phiJ2 * dev;

      // Unused for the non-associative flow rule, but left here for completeness
      // const double phiI1I1 = A * A * ( 1.0 / B - A * A * I1 * I1 / ( B * B * B ) ) / ( 2.0 * eta );
      // const double phiJ2J2 = -18.0 * eta / ( B * B * B );
      // const double phiI1J2 = -3.0 * A * A * I1 / ( B * B * B );

      /* Tensor3333d d2f_dMandel_dMandel = phiI1I1 * Spatial3D::I4 + phiI1J2 * Fastor::outer( Spatial3D::I, dev ) +
                                        phiI1J2 * Fastor::outer( dev, Spatial3D::I ) +
                                        phiJ2J2 * Fastor::outer( dev, dev ) +
                                        phiJ2 * ( Spatial3D::ISymm - 1.0 / 3.0 * Spatial3D::I4 ); */
      double df_dBetaP = -1.0;

      const double g = sqrt( std::max( 3 * J2 + thetaPlus * Math::macauly( p ) * Math::macauly( p ) +
                                         thetaMinus * Math::macauly( -p ) * Math::macauly( -p ),
                                       1e-15 ) );

      Tensor33d Q = 3.0 * dev +
                    2.0 * ( thetaPlus * Math::macauly( p ) - thetaMinus * Math::macauly( -p ) ) * Spatial3D::I / 3.0;

      Tensor3333d dQ_dMandel = 3.0 * ( Spatial3D::ISymm - Spatial3D::IHyd / 3.0 ) +
                               2.0 *
                                 ( thetaPlus * Math::heavisideExclude0( p ) +
                                   thetaMinus * Math::heavisideExclude0( -p ) ) *
                                 Spatial3D::IHyd / 9.0;

      Tensor33d   dg_dMandel          = Q / ( 2.0 * g );
      Tensor3333d d2g_dMandel_dMandel = dQ_dMandel / ( 2.0 * g ) - Fastor::outer( Q, Q ) / ( 4.0 * g * g * g );

      double    h          = sqrt( 2.0 * Fastor::inner( dg_dMandel, dg_dMandel ) / 3.0 );
      Tensor33d dh_dMandel = 2.0 * einsum< ijkl, kl, to_ij >( d2g_dMandel_dMandel, dg_dMandel ) / ( 3.0 * h );

      return { f, df_dMandel, df_dBetaP, dg_dMandel, d2g_dMandel_dMandel, h, dh_dMandel };
    }

    bool isYielding( const Tensor33d& Fe, const double betaP, const double omega )
    {
      double    f, df_dBetaP;
      Tensor33d df_dFe, dg_dMandel, dh_dFe;
      std::tie( f, df_dFe, df_dBetaP, dg_dMandel, dh_dFe ) = yieldFunction( Fe, betaP );
      if ( f > 0.0 )
        return true;
      else
        return false;
    }

    std::tuple< Tensor33d, Tensor3333d > computeMandelStress( const Tensor33d& Fe )
    {
      using namespace Marmot::ContinuumMechanics;
      Tensor33d   Ce;
      Tensor3333d dCe_dFe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      double      psi_;
      Tensor33d   dPsi_dCe;
      Tensor3333d d2Psi_dCedCe, dMandel_dCe;

      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      Tensor33d       PK2                      = 2.0 * dPsi_dCe;
      const Tensor33d mandel                   = Ce % PK2;
      dMandel_dCe                              = einsum< Ii, iJKL, to_IJKL >( Ce, 2. * d2Psi_dCedCe ) +
                    einsum< IK, iL, iJ, to_IJKL >( Spatial3D::I, Spatial3D::I, PK2 );
      Tensor3333d dMandel_dFe = einsum< IJKL, KLMN >( dMandel_dCe, dCe_dFe );
      return { mandel, dMandel_dFe };
    }

    std::tuple< Tensor33d, Tensor3333d > computeNominalMandelStress( const Tensor33d& Fe, const double& omega )
    {
      using namespace Marmot::ContinuumMechanics;
      Tensor33d   Ce;
      Tensor3333d dCe_dFe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      double      psi_;
      Tensor33d   dPsi_dCe;
      Tensor3333d d2Psi_dCedCe, dMandelN_dCe;

      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );

      Tensor33d       PK2     = 2.0 * dPsi_dCe;
      const Tensor33d mandelN = ( Ce % PK2 ) * ( 1.0 - omega );

      dMandelN_dCe = ( einsum< Ii, iJKL, to_IJKL >( Ce, 2. * d2Psi_dCedCe ) +
                       einsum< IK, iL, iJ, to_IJKL >( Spatial3D::I, Spatial3D::I, PK2 ) ) *
                     ( 1.0 - omega );

      Tensor3333d dMandelN_dFe = einsum< IJKL, KLMN >( dMandelN_dCe, dCe_dFe );

      return { mandelN, dMandelN_dFe };
    }

    std::tuple< double, double > computeBetaP( const double alphaP )
    {
      const double beta = ft + Q1 * alphaP * exp( -b1 * alphaP ) + Q2 * ( 1.0 - exp( -b2 * alphaP ) ) +
                          b3 * alphaP * alphaP * alphaP + b4 * alphaP * alphaP + b5 * alphaP;
      const double dBetaP_dAlphaP = Q1 * exp( -b1 * alphaP ) * ( 1.0 - b1 * alphaP ) + Q2 * b2 * exp( -b2 * alphaP ) +
                                    3.0 * b3 * alphaP * alphaP + 2.0 * b4 * alphaP + b5;
      return { beta, dBetaP_dAlphaP };
    }

    std::tuple< Eigen::VectorXd, Eigen::MatrixXd > computeResidualVectorAndTangent( const Eigen::VectorXd& X,
                                                                                    const Tensor33d&       FeTrial,
                                                                                    const double           alphaPTrial,
                                                                                    const double           dt )
    {

      const int idxA = 9;
      const int idxF = 10;
      using namespace Eigen;
      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;
      VectorXd R( 11 );
      MatrixXd dR_dX( 11, 11 );
      dR_dX.setZero();
      // initialize residual
      R.segment< 9 >( 0 ) = -mV9d( FeTrial.data() );
      R( 9 )              = -alphaPTrial;

      Tensor33d Fe( X.segment( 0, 9 ).data() );

      double       h;
      const double dLambda = X( 10 );
      const double alphaP  = X( 9 );

      double betaP, dBetaP_dAlphaP;
      std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaP );

      // compute mandel stress
      Tensor33d   mandelStress;
      Tensor3333d dMandel_dFe, d2g_dMandel_dMandel;
      std::tie( mandelStress, dMandel_dFe ) = computeMandelStress( Fe );

      double    f, df_dBetaP;
      Tensor33d df_dMandel, dg_dMandel, dh_dMandel, df_dFe, dh_dFe;
      std::tie( f,
                df_dMandel,
                df_dBetaP,
                dg_dMandel,
                d2g_dMandel_dMandel,
                h,
                dh_dMandel ) = yieldFunctionFromStress( mandelStress, betaP );

      Tensor33d   dGp = dLambda * dg_dMandel;
      Tensor33d   dFp;
      Tensor3333d ddFp_ddGp;
      std::tie( dFp, ddFp_ddGp ) = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::FirstOrderDerived::
        exponentialMap( dGp );

      Tensor3333d ddGp_dFe      = dLambda * einsum< ijmn, mnkL >( d2g_dMandel_dMandel, dMandel_dFe );
      Tensor33d   ddFp_ddLambda = einsum< IJKL, KL >( ddFp_ddGp, dg_dMandel );

      Tensor3333d ddFp_dFe = einsum< iImn, mnkL >( ddFp_ddGp, ddGp_dFe );

      Tensor3333d dFeTrial_dFe = einsum< iI, IJKL >( Fe, ddFp_dFe ) +
                                 einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Spatial3D::I % dFp ) );

      Tensor33d dFe_ddLambda = einsum< Ii, iJ >( Fe, ddFp_ddLambda );

      std::tie( f, df_dFe, df_dBetaP, dg_dMandel, dh_dFe ) = yieldFunction( Fe, betaP );

      double beta_min           = 1e-12;
      double sgn_beta           = ( betaP >= 0 ) ? 1.0 : -1.0;
      double betaP_cap          = sgn_beta * std::max( std::abs( betaP ), beta_min );
      double dBetaP_dAlphaP_cap = ( std::abs( betaP ) > beta_min ) ? sgn_beta * dBetaP_dAlphaP : 0.0;

      double    r      = Math::macauly( f ) / betaP_cap;
      double    D      = std::pow( Math::macauly( r ), ( 1.0 / n ) );
      Tensor33d dD_dFe = ( 1.0 / n ) * std::pow( ( Math::macauly( r ) ), ( ( 1.0 - n ) / n ) ) *
                         Math::heavisideExclude0( r ) * df_dFe * Math::heavisideExclude0( f ) / betaP_cap;
      double dD_dalphaP = ( 1.0 / n ) * std::pow( ( Math::macauly( r ) ), ( ( 1.0 - n ) / n ) ) *
                          Math::heavisideExclude0( r ) * ( -dBetaP_dAlphaP_cap ) *
                          ( Math::heavisideExclude0( f ) / betaP_cap + Math::macauly( f ) / ( betaP_cap * betaP_cap ) );

      const double    hmin  = 1e-8;
      const double    hsafe = std::max( h, hmin );
      const Tensor33d hZero( 0.0 );
      const Tensor33d dh_dMandel_safe = ( h > hmin ) ? dh_dMandel : hZero;
      // const Tensor33d dh_dFe_safe        = einsum< mn, mnij, to_ij >( dh_dMandel_safe, dMandel_dFe );

      // std::cout << "D: " << D << std::endl;
      // std::cout << "f: " << f << std::endl;
      // std::cout << "betaP: " << betaP << std::endl;
      // Residual
      R.segment< 9 >( 0 ) += mV9d( Tensor33d( einsum< iJ, JK >( Fe, dFp ) ).data() );
      R( idxA ) += ( alphaP - dLambda * hsafe );
      R( idxF ) = dLambda - dt * eta_VP * D / hsafe;

      // Jacobian
      // dR_dFe
      dR_dX.block< 9, 9 >( 0, 0 )    = mM9d( dFeTrial_dFe.data() ).transpose();
      dR_dX.block< 9, 1 >( 0, idxF ) = mV9d( dFe_ddLambda.data() );

      // dR_dalphaP
      Tensor33d dRa_dFe              = -dLambda * einsum< mn, mnij, to_ij >( dh_dMandel_safe, dMandel_dFe );
      dR_dX.block< 1, 9 >( idxA, 0 ) = mV9d( dRa_dFe.data() ).transpose();
      dR_dX( idxA, idxA )            = 1.0;
      dR_dX( idxA, idxF )            = -hsafe;

      // dR_dLambda
      Tensor33d dRl_dFe = -eta_VP * dt *
                          ( dD_dFe * hsafe - D * einsum< mn, mnij, to_ij >( dh_dMandel_safe, dMandel_dFe ) ) /
                          ( hsafe * hsafe );
      dR_dX.block< 1, 9 >( idxF, 0 ) = mV9d( dRl_dFe.data() ).transpose();
      dR_dX( idxF, idxA )            = -eta_VP * dt * dD_dalphaP / hsafe;
      dR_dX( idxF, idxF )            = 1;

      return { R, dR_dX };
    }
  };

} // namespace Marmot::Materials
