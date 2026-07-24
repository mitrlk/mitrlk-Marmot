/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
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
#include "Fastor/Fastor.h"
#include "Marmot/MarmotUtils.h"
#include <Fastor/tensor/Tensor.h>
#include <string>
#include <tuple>

/**
 * @class MarmotMaterialGradientEnhancedFiniteStrain
 * @brief Standalone base class for gradient-enhanced (nonlocal) finite-strain
 *        materials with a single nonlocal field.
 *
 * Analytic-tangent counterpart of upstream's AD-based
 * MarmotMaterialGradientEnhancedFiniteStrainAD. Self-contained: since v26
 * removed the common `MarmotMaterial` root in favour of independent per-kind
 * bases, this class inlines the material-property / state-variable boilerplate
 * itself (mirroring MarmotMaterialFiniteStrain) and adds the gradient-enhanced
 * finite-strain constitutive interface (nonlocal field A, coupled tangents).
 */
class MarmotMaterialGradientEnhancedFiniteStrain {

protected:
  const double* materialProperties;  ///< Pointer to array of material property values.
  const int     nMaterialProperties; ///< Number of material properties.

  double* stateVars = nullptr;       ///< Pointer to array of state variables.
  int     nStateVars = 0;            ///< Number of assigned state variables.

public:
  const int materialNumber; ///< Identifier for material type/implementation.

  MarmotMaterialGradientEnhancedFiniteStrain( const double* materialProperties_,
                                                int           nMaterialProperties_,
                                                int           materialNumber_ )
    : materialProperties( materialProperties_ ),
      nMaterialProperties( nMaterialProperties_ ),
      materialNumber( materialNumber_ )
  {
  }

  virtual ~MarmotMaterialGradientEnhancedFiniteStrain() = default;

  // ----------------------------------------------------------- state variables
  virtual int getNumberOfRequiredStateVars() = 0;

  virtual void assignStateVars( double* stateVars_, int nStateVars_ )
  {
    this->stateVars  = stateVars_;
    this->nStateVars = nStateVars_;
  }

  double* getAssignedStateVars() { return stateVars; }
  int     getNumberOfAssignedStateVars() { return nStateVars; }

  virtual StateView getStateView( const std::string& stateName ) = 0;

  virtual void initializeYourself()
  {
    for ( int i = 0; i < this->getNumberOfAssignedStateVars(); i++ )
      this->stateVars[i] = 0;
  }

  virtual double getDensity() { return 1.0; }

  // -------------------------------------------------- gradient-enhanced response
  template < int nDim >
  struct ConstitutiveResponse {
    Fastor::Tensor< double, nDim, nDim > tau;                  ///< Kirchhoff stress
    double                               rho;                  ///< mass density
    double                               elasticEnergyDensity; ///< elastic energy per unit volume
    double                               nonlocalradius;       ///< nonlocal interaction radius
    double                               L;                    ///< Local field
  };

  template < int nDim >
  struct AlgorithmicModuli {
    Fastor::Tensor< double, nDim, nDim, nDim, nDim > dTau_dF; ///< tangent operator w.r.t. deformation gradient
    Fastor::Tensor< double, nDim, nDim >             dTau_dA; ///< tangent operator w.r.t. nonlocal field
    Fastor::Tensor< double, nDim, nDim > dL_dF; ///< tangent operator of local field w.r.t. deformation gradient
  };

  template < int nDim >
  struct Deformation {
    Fastor::Tensor< double, nDim, nDim > F; ///< deformation gradient
    double                               A; ///< nonlocal field
  };

  struct TimeIncrement {
    const double time; ///< time at the beginning of the increment
    const double dT;   ///< size of the time increment
  };

  virtual void computeStress( ConstitutiveResponse< 3 >& response,
                              AlgorithmicModuli< 3 >&    tangents,
                              const Deformation< 3 >&,
                              const TimeIncrement& ) = 0;

  virtual void computeStress( ConstitutiveResponse< 3 >&                  response,
                              AlgorithmicModuli< 3 >&                     tangents,
                              const Deformation< 3 >&                     deformation,
                              const TimeIncrement&                        timeIncrement,
                              const std::tuple< double, double, double >& eigenDeformation );

  virtual void computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    algorithmicModuli,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement );

  virtual void computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                   AlgorithmicModuli< 3 >&                     algorithmicModuli,
                                   const Deformation< 3 >&                     deformation,
                                   const TimeIncrement&                        timeIncrement,
                                   const std::tuple< double, double, double >& eigenDeformation );

  virtual void computePlaneStress( ConstitutiveResponse< 2 >& response,
                                   AlgorithmicModuli< 2 >&    algorithmicModuli,
                                   const Deformation< 2 >&    deformation,
                                   const TimeIncrement&       timeIncrement );

  std::tuple< double, double, double > findEigenDeformationForEigenStress(
    const std::tuple< double, double, double >& initialGuess,
    const std::tuple< double, double, double >& eigenStress );
};
