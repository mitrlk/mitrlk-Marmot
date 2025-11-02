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
#include "Marmot/MarmotMaterialGEFiniteStrain.h"
#include "Marmot/MarmotJournal.h"

void MarmotMaterialGEFiniteStrain::computeStress( ConstitutiveResponse< 3 >&                  response,
                                                  AlgorithmicModuli< 3 >&                     tangents,
                                                  const Deformation< 3 >&                     deformation,
                                                  const TimeIncrement&                        timeIncrement,
                                                  const std::tuple< double, double, double >& eigenDeformation )
{
  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                            << "Not yet implemented for gradient-enhanced finite strain materials." );
}

void MarmotMaterialGEFiniteStrain::computePlaneStrain( ConstitutiveResponse< 3 >& response,
                                                       AlgorithmicModuli< 3 >&    algorithmicModuli,
                                                       const Deformation< 3 >&    deformation,
                                                       const TimeIncrement&       timeIncrement )
{
  return computeStress( response, algorithmicModuli, deformation, timeIncrement );
}

void MarmotMaterialGEFiniteStrain::computePlaneStrain( ConstitutiveResponse< 3 >&                  response,
                                                       AlgorithmicModuli< 3 >&                     algorithmicModuli,
                                                       const Deformation< 3 >&                     deformation,
                                                       const TimeIncrement&                        timeIncrement,
                                                       const std::tuple< double, double, double >& eigenDeformation )
{
  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                            << "Not yet implemented for gradient-enhanced finite strain materials." );
}

void MarmotMaterialGEFiniteStrain::computePlaneStress( ConstitutiveResponse< 2 >& response,
                                                       AlgorithmicModuli< 2 >&    algorithmicModuli,
                                                       const Deformation< 2 >&    deformation,
                                                       const TimeIncrement&       timeIncrement )
{
  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                            << "Not yet implemented for gradient-enhanced finite strain materials." );
}

std::tuple< double, double, double > MarmotMaterialGEFiniteStrain::findEigenDeformationForEigenStress(
  const std::tuple< double, double, double >& initialGuess,
  const std::tuple< double, double, double >& eigenStress )
{
  throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                            << "Not yet implemented for gradient-enhanced finite strain materials." );
}
