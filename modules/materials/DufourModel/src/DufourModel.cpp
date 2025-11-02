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

    response.nonlocalradius = ld;
    double alphaP_nonlocal  = deformation.A;

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
      response.L      = alphaP;
      Tensor33d FpNew = dFp % Fp;
      memcpy( Fp.data(), FpNew.data(), 9 * sizeof( double ) );

      using namespace ContinuumMechanics;
      double      psi_, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal;
      Tensor33d   Ce, dPsi_dCe, tau_eff;
      Tensor3333d dCe_dFe, d2Psi_dCedCe, dTau_dPK2_eff, dTau_dFe_partial_eff;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      // compute damage variable
      std::tie( omega, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal ) = computeOmega( alphaP, alphaP_nonlocal );

      // compute Kirchhoff stress
      Tensor33d PK2_eff = 2.0 * dPsi_dCe;

      std::tie( tau_eff,
                dTau_dPK2_eff,
                dTau_dFe_partial_eff ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2_eff, Fe );

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
      MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

      Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

      Tensor3333d dPK2_dFe_eff = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dPK2_dF_eff  = einsum< ijKL, KLMN >( dPK2_dFe_eff, dFe_dF );

      /* tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dPK2, dPK2_dF ) +
       * dTau_dF_partial; */
      Tensor3333d dTau_dF_eff = einsum< IJKL, KLMN >( dTau_dPK2_eff, dPK2_dF_eff ) +
                                einsum< ijKL, KLMN >( dTau_dFe_partial_eff, dFe_dF );

      Tensor33d dAlphaP_local_dF = Tensor33d( Vector9d( dXdDeformation.block< 1, 9 >( 9, 0 ).transpose() ).data() );

      tangents.dTau_dF = ( 1 - omega ) * dTau_dF_eff -
                         dOmega_dAlphaP_local * Fastor::outer( tau_eff, dAlphaP_local_dF );
      tangents.dTau_dA = -tau_eff * dOmega_dAlphaP_nonlocal;
      tangents.dL_dF   = dAlphaP_local_dF;
    }
    else {
      using namespace Marmot::ContinuumMechanics;
      double      psi_, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal;
      Tensor33d   Ce, dPsi_dCe, tau_eff;
      Tensor3333d dCe_dFe, d2Psi_dCedCe, dTau_dPK2_eff, dTau_dFe_partial_eff;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy
      // Green deformation
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      std::tie( omega, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal ) = computeOmega( alphaPOld, alphaP_nonlocal );

      // compute Kirchhoff stress
      Tensor33d PK2_eff = 2. * dPsi_dCe;

      std::tie( tau_eff,
                dTau_dPK2_eff,
                dTau_dFe_partial_eff ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2_eff, Fe );

      response.tau                  = tau_eff * ( 1.0 - omega );
      response.rho                  = density;
      response.elasticEnergyDensity = psi_;
      response.L                    = alphaP;

      // compute tangent operator
      Tensor3333d dPK2_dFe    = einsum< ijKL, KLMN >( 2. * d2Psi_dCedCe, dCe_dFe );
      Tensor3333d dFe_dF      = einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Fastor::inverse( FpOld ) ) );
      Tensor3333d dPK2_dF_eff = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      Tensor3333d dTau_dF_eff = einsum< IJKL, KLMN >( dTau_dPK2_eff, dPK2_dF_eff ) +
                                einsum< ijKL, KLMN >( dTau_dFe_partial_eff, dFe_dF );

      Tensor33d dAlphaP_local_dF = Tensor33d( 0.0 );

      tangents.dTau_dF = ( 1 - omega ) * dTau_dF_eff;
      tangents.dTau_dA = -tau_eff * dOmega_dAlphaP_nonlocal;
      tangents.dL_dF   = dAlphaP_local_dF;
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
