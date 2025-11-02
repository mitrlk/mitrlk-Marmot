#include "Marmot/DufourModel.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialGEFiniteStrain.h"
#include "Marmot/MarmotStressMeasures.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  DufourModel::DufourModel( const double* materialProperties, int nMaterialProperties, int materialLabel )
    : MarmotMaterialGEFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
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
      density( nMaterialProperties > 19 ? materialProperties[19] : 0.0 )
  {
  }

  void DufourModel::computeStress( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    tangents,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement )
  {

    auto&           Fp = stateVars->Fp;
    const Tensor33d FpOld( Fp );
    double&         alphaP    = stateVars->alphaP;
    const double    alphaPOld = alphaP;
    double&         omega     = stateVars->omega;
    const double    omegaOld  = omega;

    response.nonLocalRadius = ld;
    double alphaP_nonlocal  = deformation.A;
    response.L              = alphaP;

    using namespace Marmot;
    using namespace Fastor;
    using namespace Eigen;
    using namespace autodiff;
    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    Tensor33d FeTrial = deformation.F % Fastor::inverse( FpOld );
    double    betaP, dBetaP_dAlphaP;
    std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaPOld );

    Tensor33d dFp;
    dFp.eye();
    Tensor33d Fe = FeTrial;

    if ( isYielding( FeTrial, betaP, omegaOld ) ) {

      size_t counter = 0;

      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      VectorXd X( 11 );
      X.segment( 0, 9 )  = mV9d( FeTrial.data() );
      X( 9 )             = alphaPOld;
      X( 10 )            = 0.0;
      VectorXd        dX = VectorXd::Zero( 11 );
      VectorXd        R  = VectorXd::Zero( 11 );
      Eigen::MatrixXd dR_dX( 11, 11 );

      std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial, alphaPOld, timeIncrement.dT );

      while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {

        if ( counter > 15 ) {
          std::cout << "R:" << R.transpose() << '\n'
                    << "R_norm" << R.norm() << '\n'
                    << "X:" << X.transpose() << '\n'
                    << "dX norm" << dX.norm() << '\n'
                    << "dR_dX" << '\n'
                    << dR_dX << std::endl;
          throw std::runtime_error( "inner newton not converged" );
        }

        dX = -dR_dX.colPivHouseholderQr().solve( R );
        X += dX;
        std::tie( R, dR_dX ) = computeResidualVectorAndTangent( X, FeTrial, alphaPOld, timeIncrement.dT );
        counter += 1;
      }
      // update plastic deformation increment
      Fe              = X.segment( 0, 9 ).data();
      dFp             = Fastor::inverse( Fe ) % FeTrial;
      alphaP          = X( 9 );
      Tensor33d FpNew = dFp % Fp;
      memcpy( Fp.data(), FpNew.data(), 9 * sizeof( double ) );

      using namespace ContinuumMechanics;
      double      psi_, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      // compute damage variable
      std::tie( omega, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal ) = computeOmega( alphaP, alphaP_nonlocal );

      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe * ( 1.0 - omega );
      Tensor3333d dTau_dPK2eff, dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2eff,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;
      dTau_dPK2                     = dTau_dPK2eff * ( 1.0 - omega );

      // compute tangent operator
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

      MatrixXd dYdDeformation              = MatrixXd::Zero( 11, 11 );
      dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d( einsum< IK, JL, to_IJKL >( Spatial3D::I,
                                                                                           transpose( Fastor::inverse(
                                                                                             FpOld ) ) ) )
                                                     .data() )
                                               .transpose();
      MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

      Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe ) * ( 1.0 - omega );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
                         einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );                                 // adjust
      tangents.dTau_dA = -einsum< IJKL >( dTau_dPK2, 2. * dPsi_dCe ) * ( dOmega_dAlphaP_nonlocal ) / epsF; // check
      tangents.dL_dF   = -einsum< IJKL >( dTau_dPK2, 2. * dPsi_dCe ) * ( dAlphaP_weighted_dAlphaP_nonlocal ) /
                       epsF;                                                                               // check
    }
    else {
      using namespace Marmot::ContinuumMechanics;
      double      psi_, dAlphaP_weighted_dAlphaP_local, dAlphaP_weighted_dAlphaP_nonlocal;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      std::tie( omega,
                dAlphaP_weighted_dAlphaP_local,
                dAlphaP_weighted_dAlphaP_nonlocal ) = computeOmega( alphaPOld, alphaP_nonlocal );

      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe * ( 1.0 - omega );
      Tensor3333d dTau_dPK2eff, dTau_dPK2, dTau_dFe_partial;
      std::tie( response.tau,
                dTau_dPK2eff,
                dTau_dFe_partial )  = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2, Fe );
      response.rho                  = 1.0;
      response.elasticEnergyDensity = psi_;
      dTau_dPK2                     = dTau_dPK2eff * ( 1.0 - omega );

      // compute tangent operator
      Tensor3333d dPK2_dFe = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe ) * ( 1.0 - omega );
      Tensor3333d dFe_dF   = einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Fastor::inverse( FpOld ) ) );
      Tensor3333d dPK2_dF  = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) + einsum< ijKL, KLMN >( dTau_dFe_partial, dFe_dF );
      tangents.dTau_dA = -einsum< IJKL >( dTau_dPK2, 2. * dPsi_dCe ) * ( dAlphaP_weighted_dAlphaP_nonlocal ) /
                         epsF;                                                                                  // check
      tangents.dL_dF = -einsum< IJKL >( dTau_dPK2, 2. * dPsi_dCe ) * ( dAlphaP_weighted_dAlphaP_local ) / epsF; // check
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
  }
} // namespace Marmot::Materials
