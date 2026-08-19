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
#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include <array>
#include <cstring>
#include <string>
#include <vector>

namespace Marmot::Materials {

  using namespace Fastor;
  using namespace FastorStandardTensors;
  using namespace FastorIndices;

  class DufourModel : public MarmotMaterialGradientEnhancedFiniteStrain {
  public:
    using MarmotMaterialGradientEnhancedFiniteStrain::MarmotMaterialGradientEnhancedFiniteStrain;

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

    // Stress-weighted ductile fracture (SWDFM) initiation driver -- OPTIONAL card entries 21-24.
    //
    // Abrari Vajari, Neuner, Kammardi Arunachala, Ziccarelli, Deierlein, Linder,
    // CMAME 400 (2022) 115467, eq. (65)-(66), after Rice & Tracey and Smith et al.:
    //
    //   D = cSW * int [ exp( 1.3 T ) - 1/bSW * exp( -1.3 T ) ] * exp( kSW ( |zeta| - 1 ) ) dAlphaP
    //
    // with T the stress triaxiality and zeta = cos( 3 theta ) the Lode angle PARAMETER.
    // Crack initiation at D = 1; dF then governs how fast omega grows past initiation.
    //
    // cSW = 0 reproduces the unscaled original model EXACTLY (see computeOmega).
    const double cSW, bSW, dF;

    // ------------------------------------------------------------------------------------------
    // CUBIC-EXPONENT Rice-Tracey weight in the SWDFM driver -- OPTIONAL, the three card entries
    // that follow the Prony triplets, i.e. [28 + 3 nMaxwell] .. [30 + 3 nMaxwell]:
    //
    //   [28 + 3 nMaxwell]  swdfmC2     coefficient of T^2
    //   [29 + 3 nMaxwell]  swdfmC3     coefficient of T^3
    //   [30 + 3 nMaxwell]  swdfmC1     coefficient of T   (0 or absent -> swdfmExponent = 1.3)
    //
    //   g(T) = exp( c1 T + c2 T^2 + c3 T^3 ) - exp( -( c1 T + c2 T^2 + c3 T^3 ) ) / bSW
    //
    // All three zero (or absent) reproduces the published single-branch law exp( 1.3 T ) bit for
    // bit, so every deck that does not carry the tail is unaffected.
    //
    // CALIBRATED VALUE (HANDOFF S8.4, `calibration/refit_g.py`, no FE required):
    //
    //   g(T) = exp( 3.75 T - 5.75 T^2 + 3.50 T^3 ),  cSW = 1.797
    //
    // i.e. c1 = 3.75, c2 = -5.75, c3 = 3.50. Monotone increasing over T = 0..1.3 (verified
    // numerically; min slope 1.5e-3), so it is admissible as a void-growth weight.
    //
    // WHY THIS REPLACED THE TWO-BRANCH KINK. The required constant cSW was measured on all nine
    // Arcan runs plus the single-lap joint (ten runs). g must go as 1/kappa_bar_f: 2.05 at
    // T = 0.07, 4.65 at 0.81, 12.39 at 1.08, i.e. a local exponent of 1.12 below T ~ 0.8 and 3.51
    // above it. The superseded two-branch law reproduced that with a KINK at swdfmT0 = 0.88 --
    // and the SLJ's damage-critical element sits at T ~ 0.59, so 0 % of its driver accumulated
    // above the threshold and the whole term did nothing for the one specimen that was not in the
    // fit (HANDOFF S7.2). One smooth rising exponent does the same job everywhere and needs no
    // threshold. Measured spread of the required cSW over all ten runs: 2.42x -> 1.85x, against a
    // within-angle rate-scatter floor of 1.18-1.37x. The SLJ then requires 2.53, the same as
    // Arcan 0 deg, and is no longer the outlier.
    //
    // DOUBLE-COUNTING GATE (ledger 21.5 gate 3): this makes the DAMAGE law the owner of the
    // cavitation mechanism. If a tension-side cap is ever added to the yield surface as well,
    // one of the two must be removed -- nuP_plus already owns the plastic flow DIRECTION.
    const double swdfmC2, swdfmC3, swdfmC1;

    // ------------------------------------------------------------------------------------------
    // MONOTONE-BY-CONSTRUCTION exponent + RATE-DEPENDENT T-SENSITIVITY -- OPTIONAL, the five card
    // entries after the cubic ones, i.e. [31 + 3 nMaxwell] .. [35 + 3 nMaxwell]:
    //
    //   [31 + 3 nMaxwell]  swdfmB0 |  E'(T) = ( b0 + b1 T + b2 T^2 )^2, so E is non-decreasing
    //   [32 + 3 nMaxwell]  swdfmB1 |  for ANY parameter values -- monotonicity is structural, not
    //   [33 + 3 nMaxwell]  swdfmB2 |  a fitted accident. All three zero -> use the cubic above.
    //   [34 + 3 nMaxwell]  swdfmKdotRef   reference rate [1/s]; 0 (or absent) -> NO rate term
    //   [35 + 3 nMaxwell]  swdfmS         rate exponent; NEGATIVE -> tie it to n (see below)
    //
    //   g(T, kdot) = exp( E(T) * w ),    w = ( kdot / swdfmKdotRef )^( -s )
    //   E(T) = b0^2 T + b0 b1 T^2 + (b1^2 + 2 b0 b2)/3 T^3 + (b1 b2)/2 T^4 + (b2^2)/5 T^5
    //
    // with kdot = d alphaPBar / dt, the rate at which the driver itself accumulates. w > 0 always,
    // so g stays non-decreasing in T at EVERY rate.
    //
    // CALIBRATED (HANDOFF S9.10, ten runs: 9 Arcan + the single-lap joint):
    //   ( b0, b1, b2 ) = ( 1.798, -0.702, -3.172 ),  swdfmKdotRef = 20,  swdfmS < 0,  cSW = 2.006
    // Spread of the required constant: 2.42x (published) -> 1.85x (cubic) -> 1.43x, against a floor
    // of 1.18-1.37x set by the within-angle scatter across rates, which is experimental.
    //
    // WHY swdfmS < 0 MEANS "USE n". The freely fitted rate exponent is 0.0502; the model's own
    // VISCOPLASTIC exponent n, identified by Dufour from three loading rates, is 0.0435. Fixing s
    // to n costs 0.01x of spread (1.43x against 1.42x). So the rate dependence of FAILURE is the
    // rate dependence the FLOW RULE already has, and it costs NO new parameter. Encoding it as a
    // sentinel rather than as a repeated literal keeps the two tied if n is ever recalibrated.
    //
    // RATE FLOOR. On an elastic increment d alphaPBar = 0, so kdot = 0 and w would be infinite;
    // the increment contributes nothing to D, but 0 * inf is NaN, so kdot is floored at
    // swdfmKdotMin. The floor also bounds w: at 1e-8 /s against a reference of 20 it gives
    // w <= 2.5, and the slowest real run (the SLJ at 1 mm/min) sits at 1e-3 /s, five decades above.
    //
    // MEASURED SPAN, so it is clear this is not a lever on one specimen: w = 1.535 / 1.105 / 1.000
    // / 0.905 at the SLJ / 1 / 10 / 100 mm/s -- a factor 1.7 across FOUR decades. A rejected
    // candidate with a log^2 rate factor scored better but applied x3.53 to the SLJ alone and
    // x1.07 to all nine Arcan runs (ledger 25.10); that is why the form here acts on the
    // T-SENSITIVITY and not as a common factor.
    //
    // DO NOT EXTRAPOLATE IN T. g = 1.00 / 1.77 / 2.48 / 2.72 / 2.80 / 4.30 / 22.0 at
    // T = 0 / 0.2 / ... / 1.2 while the data reaches only T = 1.09 (Arcan 0 deg). The quintic rises
    // very steeply past that; cap it before using this card on a new geometry.
    const double swdfmB0, swdfmB1, swdfmB2, swdfmKdotRef, swdfmS;

    // ------------------------------------------------------------------------------------------
    // QUARTIC term and the EXTRAPOLATION CAP -- OPTIONAL, entries [37 + 3 nMaxwell] and
    // [38 + 3 nMaxwell]:
    //
    //   [37 + 3 nMaxwell]  swdfmC4    coefficient of T^4 in the cubic-form exponent
    //   [38 + 3 nMaxwell]  swdfmTCap  triaxiality above which g is HELD CONSTANT; 0 -> no cap
    //
    //   E(T) = c1 T + c2 T^2 + c3 T^3 + c4 T^4 ,   evaluated at min(T, swdfmTCap)
    //
    // WHY THE QUARTIC. The required weight is NON-MONOTONE (ledger 28.2): at the cSW that Arcan
    // 90 deg pins (1.88, because 90 deg sits at T ~ 0 where g == 1) the four specimens demand
    // g(0) = 1.01, g(0.584) = 3.38 (SLJ), g(0.806) = 2.62 (Arcan 45 deg), g(1.087) = 5.66 (0 deg).
    // g must PEAK near T ~ 0.6 and DIP by T ~ 0.85. A cubic cannot hold that and still rise steeply
    // afterwards; the calibrated quartic
    //
    //   E(T) = -0.07 T + 18.12 T^2 - 35.58 T^3 + 18.85 T^4 ,   cSW = 1.754
    //
    // turns at T = 0.559 and 0.854 -- essentially exactly the SLJ's hot element and Arcan 45 deg --
    // and collapses the four SPECIMEN medians to 1.11x, below the 1.18-1.37x rate-scatter floor.
    // The physical reading is cavitation of the rubber phase: it needs hydrostatic tension to switch
    // on but is suppressed again at high constraint, so it has a WINDOW, superposed on Rice-Tracey
    // void growth at high T. Monotone forms are structurally incapable and were the reason nothing
    // worked; the monotonicity rule of ledger 25.6 was itself the error.
    //
    // WHY THE CAP IS NOT OPTIONAL IN PRACTICE. g is constrained by data at FOUR triaxialities only
    // (0, 0.584, 0.806, 1.087). Beyond the last one the quartic explodes:
    //   g(1.087) = 7.0   g(1.15) = 15.5   g(1.20) = 37.1   g(1.30) = 493   g(1.50) = 9.0e6
    // A thin adhesive layer between stiff substrates in pure tension -- a BUTT JOINT -- approaches
    // hydrostatic tension and sits well above T = 1.1. Uncapped, this card would give such a
    // geometry a damage rate thousands of times too high and fail it at first load. With
    // swdfmTCap = 1.10 the model instead HOLDS the last value it has evidence for, which is a
    // statement about the calibration range rather than an invented trend.
    const double swdfmTCap;

    /** omega_f from the driver. ONE definition, used both by computeOmega and by the lagged
     *  interaction gate, so the two cannot drift apart. */
    double omegaFOfDriver( const double D ) const { return D > 1.0 ? 1.0 - std::exp( -( D - 1.0 ) / dF ) : 0.0; }

    /** Softening variable and its derivative w.r.t. the (weighted) driving strain.
     *
     * ONE function so the law, its derivative and the tests cannot drift apart.
     *   omega_s = Ds_inf ( 1 - E ),  E = exp( -x - beta x^2 ),  x = kappa_bar / epsF
     *   d omega_s / d kappa_bar = Ds_inf E ( 1 + 2 beta x ) / epsF
     * At x = 0 that is Ds_inf / epsF for every beta -- the property the whole choice rests on.
     */
    void softening( const double kbar, double& omega_s, double& dOmega_s ) const
    {
      const double x = kbar / epsF;
      const double E = std::exp( -x );
      omega_s        = 1.0 - E;
      dOmega_s       = E / epsF;
    }

    inline const static double swdfmKdotMin = 1e-8;

    /// Read one of the optional shape entries, which sit AFTER the Prony triplets.
    static double swdfmExtra( const double* materialProperties, int nMaterialProperties, int which )
    {
      const int n = nMaterialProperties > 27 ? static_cast< int >( materialProperties[27] ) : 0;
      const int i = 28 + 3 * ( n > 0 ? n : 0 ) + which;
      return nMaterialProperties > i ? materialProperties[i] : 0.0;
    }

    /// True when the card carries the monotone-by-construction quintic instead of the cubic.
    bool swdfmMonotoneForm() const { return swdfmB0 != 0.0 || swdfmB1 != 0.0 || swdfmB2 != 0.0; }

    /** Exponent of the Rice-Tracey weight, including the rate factor.
     *
     * Kept as ONE function so that g, its unit tests and its documentation cannot drift apart, and
     * so that the fallbacks (c1 = swdfmExponent, s = n) each live in exactly one place.
     *
     * @param T     stress triaxiality, already clamped by the caller
     * @param kdot  d alphaPBar / dt of THIS increment; floored internally
     */
    double swdfmLogG( const double T, const double kdot = 0.0, double* kdotDLogGDKdot = nullptr ) const
    {
      if ( kdotDLogGDKdot )
        *kdotDLogGDKdot = 0.0;
      // COMPRESSION BRANCH. A single polynomial cannot behave on both sides of T = 0: for T < 0 the
      // even powers stay positive while the odd ones flip, so the calibrated tension coefficients give
      // g(-0.29) = 16.6 -- damage 17x FASTER in compression than in pure shear. That destroyed Arcan
      // 90 deg (its field spans T = -0.292 .. +0.002, median -0.003) and the SLJ (median cell at
      // T = -0.044). Ledger 28.6. So compression gets its own decaying branch, using the published
      // Rice-Tracey exponent and NO new parameter:
      //
      //     T <= 0 :  E = swdfmExponent * T        ( = 1.3 T, so g < 1 and falling )
      //
      // Continuous at T = 0 with g(0) = 1 either way, which is what lets Arcan 90 deg -- the only
      // specimen sitting at T ~ 0, with its driver uniform to max/median = 1.0 -- set the drive scale.
      // Applies to BOTH exponent forms; a purely monotone card is unaffected because its own value at
      // T < 0 is already below 1.
      if ( T <= 0.0 )
        return swdfmExponent * T;

      // CAP: above the calibrated range the weight is held at its last evidenced value.
      const double Tc = swdfmTCap > 0.0 ? std::min( T, swdfmTCap ) : T;

      double E;
      if ( swdfmMonotoneForm() ) {
        // E(T) = INT_0^T ( b0 + b1 s + b2 s^2 )^2 ds, expanded
        const double b0 = swdfmB0, b1 = swdfmB1, b2 = swdfmB2;
        E = ( ( ( ( b2 * b2 / 5.0 ) * Tc + b1 * b2 / 2.0 ) * Tc + ( b1 * b1 + 2.0 * b0 * b2 ) / 3.0 ) * Tc + b0 * b1 ) *
              Tc * Tc +
            b0 * b0 * Tc;
      }
      else {
        const double c1 = swdfmC1 != 0.0 ? swdfmC1 : swdfmExponent;
        E               = ( ( swdfmC3 * Tc + swdfmC2 ) * Tc + c1 ) * Tc;
      }

      if ( swdfmKdotRef > 0.0 ) {
        const double s = swdfmS < 0.0 ? n : swdfmS; // sentinel: tie the rate exponent to n
        E *= std::pow( std::max( kdot, swdfmKdotMin ) / swdfmKdotRef, -s );
        // w = (kdot/kdotRef)^(-s)  =>  kdot dE/dkdot = -s E exactly. Reported so that the
        // tangent can carry the rate sensitivity instead of omitting it: at the calibrated
        // s = n = 0.0435 the omission was mild (w spans 1.7x over four decades) but at
        // s >= 0.15 it makes the return map diverge (SLJ kC2s150 died at step 12).
        // Above the floor the derivative of max(kdot, kdotMin) is 1; below it kdot is
        // clamped, the increment is elastic and contributes nothing, so 0 is correct.
        if ( kdotDLogGDKdot && kdot > swdfmKdotMin )
          *kdotDLogGDKdot = -s * E;
      }
      return E;
    }

    // ------------------------------------------------------------------------------------------
    // Generalized-Maxwell (Prony) VISCOELASTICITY -- optional card entries 28 onwards.
    //
    // After Nguyen, Lani, Pardoen, Morelle & Noels, Int. J. Solids Struct. 96 (2016), Sec. 3.1.3:
    // a hyperelastic equilibrium spring in parallel with nMaxwell Maxwell branches. The branch
    // update is the recursive (history-free) exponential scheme of Simo (1987), reused verbatim
    // from Marmot's own ContinuumMechanics::FiniteStrain::Viscoelasticity helper, which is also
    // what CompressibleFiniteStrainLinearViscoelasticity uses.
    //
    // The relaxation acts on the SECOND PIOLA-KIRCHHOFF stress of the elastic (Pence-Gou B)
    // potential, i.e. UPSTREAM of the return map: the yield function therefore sees the relaxed
    // stress, exactly as in the existing hyperelastic-viscoplastic split.
    //
    // Card layout, per branch i = 0 .. nMaxwell-1:
    //     [27]          nMaxwell            number of branches (0..nMaxwellMax)
    //     [28 + 3 i]    gammaG_i            DEVIATORIC relative modulus of branch i
    //     [29 + 3 i]    gammaK_i            VOLUMETRIC relative modulus of branch i
    //     [30 + 3 i]    tau_i               relaxation time of branch i  [s]
    //
    // NOTE on the parameter convention: gammaG_i / gammaK_i are RELATIVE (dimensionless) moduli,
    // gammaG_i = G_i / G_0, not absolute stiffnesses. The equilibrium branch carries the
    // remainder ( 1 - sum_i gamma_i ), which is the convention the Marmot helper requires and
    // which makes the reduction exact: nMaxwell = 0, or all gamma_i = 0, reproduces the current
    // purely hyperelastic-viscoplastic model bit for bit.
    //
    // Deviatoric and volumetric relaxation are kept SEPARATE (own gamma, shared tau) because a
    // structural adhesive does not creep equally in shear and in bulk. Setting gammaK_i = 0
    // gives shear-only creep; setting gammaK_i = gammaG_i gives proportional relaxation of the
    // whole PK2 tensor, which is what the in-tree VE material does.
    //
    // COMPILE-TIME ceiling on the number of branches; it only sizes the static state-var layout.
    // The number actually used is the CARD entry nMaxwell (<= this). Set to 7 because Dufour
    // PUBLISHES a 7-term series for SikaPower-498 (IJAA 2016, Tab. 1) -- these are measured
    // parameters, not quantities identified from our three loading rates, so the usual
    // "n rates cannot identify more than ~n branches" argument does not apply and there is no
    // reason to condense the series and pay the approximation error.
    inline const static int nMaxwellMax = 7;

    const int nMaxwell;

    /// Deviatoric branch set: weights gammaG_i, relaxation times tau_i.
    const ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties maxwellDev;
    /// Volumetric branch set: weights gammaK_i, the same relaxation times tau_i.
    const ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties maxwellVol;

    /// Time increment of the step currently being solved. Set once at the top of computeStress so
    /// that the Maxwell update is available to computeMandelStress inside the return-map Newton
    /// loop, which has no access to the TimeIncrement.
    double dTCurrent = 0.0;

    /** Build a MaxwellProperties set from the (gammaG, gammaK, tau) triplets on the card.
     *
     * @param which 0 -> deviatoric weights (gammaG), 1 -> volumetric weights (gammaK).
     *
     * The Marmot helper expects an interleaved (gamma, tau) pair vector, so the triplets are
     * repacked here; the validation (tau > 0, nMaxwell >= 0) is left to createMaxwellProperties.
     */
    static ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties makeMaxwellProperties(
      const double* materialProperties,
      int           nMaterialProperties,
      int           which )
    {
      const int n = nMaterialProperties > 27 ? static_cast< int >( materialProperties[27] ) : 0;
      if ( n <= 0 )
        return ContinuumMechanics::FiniteStrain::Viscoelasticity::createMaxwellProperties( 0, nullptr );
      if ( n > nMaxwellMax )
        throw std::invalid_argument( "DufourModel: too many Maxwell branches for the state-var layout" );
      if ( nMaterialProperties < 28 + 3 * n )
        throw std::invalid_argument( "DufourModel: incomplete (gammaG, gammaK, tau) triplets on the material card" );

      std::vector< double > pairs( 2 * n );
      for ( int i = 0; i < n; i++ ) {
        pairs[2 * i]     = materialProperties[28 + 3 * i + which]; // gammaG_i or gammaK_i
        pairs[2 * i + 1] = materialProperties[30 + 3 * i];         // tau_i (shared)
      }
      return ContinuumMechanics::FiniteStrain::Viscoelasticity::createMaxwellProperties( n, pairs.data() );
    }

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
        // accumulated damage driver  D = int d(alphaP_weighted) / epsF_eff( eta ),
        // and the previous converged alphaP_weighted needed to form its increment.
        // Only used when c_eta != 0; both stay consistent with omega otherwise.
        { .name = "damageDriver", .length = 1 },
        { .name = "alphaPBar", .length = 1 },
        // ---- viscoelasticity -------------------------------------------------------------
        // PK2Ref  : the hyperelastic (unrelaxed) PK2 stress of the last CONVERGED increment.
        //           The Maxwell recurrence is driven by its increment, so it has to persist.
        // veDev   : deviatoric branch stresses Q_i, 9 doubles each.
        // veVol   : volumetric branch stresses q_i. One SCALAR each: the volumetric branch
        //           stress is spherical by construction, so only its (0,0) entry is stored.
        // Allocated for nMaxwellMax branches unconditionally so that the layout stays static;
        // unused branches simply remain zero. 79 doubles when nMaxwellMax = 7.
        { .name = "PK2Ref", .length = 9 },
        { .name = "veDev", .length = nMaxwellMax * 9 },
        { .name = "veVol", .length = nMaxwellMax },
      } );

      Fastor::TensorMap< double, 3, 3 > Fp;
      double&                           alphaP;
      double&                           omega;
      double&                           damageDriver;
      double&                           alphaPBar;
      Fastor::TensorMap< double, 3, 3 > PK2Ref;
      double*                           veDev;
      double*                           veVol;

      DufourModelStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          Fp( &find( "Fp" ) ),
          alphaP( find( "alphaP" ) ),
          omega( find( "omega" ) ),
          damageDriver( find( "damageDriver" ) ),
          alphaPBar( find( "alphaPBar" ) ),
          PK2Ref( &find( "PK2Ref" ) ),
          veDev( &find( "veDev" ) ),
          veVol( &find( "veVol" ) ){};
    };
    std::unique_ptr< DufourModelStateVarManager > stateVars;

    void assignStateVars( double* stateVars, int nStateVars );

    StateView getStateView( const std::string& result );

    void initializeYourself();

    // ------------------------------------------------------------
    // Damage functions
    // ------------------------------------------------------------

    /** Stress triaxiality eta = p/q of a Kirchhoff stress.
     *
     * eta is a RATIO of invariants, so it is identical for the Kirchhoff, Cauchy and effective
     * (undamaged) stress - the factors J and (1-omega) cancel. It may therefore be evaluated on
     * whichever measure is at hand.
     */
    static double triaxiality( const Tensor33d& tau )
    {
      const double p = ( tau( 0, 0 ) + tau( 1, 1 ) + tau( 2, 2 ) ) / 3.0;

      double J2 = 0.0;
      for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 3; j++ ) {
          const double s_ij = tau( i, j ) - ( i == j ? p : 0.0 );
          J2 += s_ij * s_ij;
        }
      }
      J2 *= 0.5;

      const double q = std::sqrt( std::max( 3.0 * J2, 0.0 ) );

      // near-hydrostatic or unstressed: eta is meaningless, fall back to 0 (shear-like)
      if ( q < 1e-12 * std::max( 1.0, std::abs( p ) ) ) {
        return 0.0;
      }
      return p / q;
    }

    /** Lode angle PARAMETER zeta = cos( 3 theta ) = 3 sqrt(3) J3 / ( 2 J2^(3/2) ).
     *
     * zeta = +1 axisymmetric tension, 0 pure shear, -1 axisymmetric compression. This is the
     * normalisation used in CMAME 400 (2022) 115467 eq. (64); note it differs from
     * thetaBar = 1 - 6 theta / pi (same endpoints, monotonically related, not equal).
     * No arccos is needed, so it is cheaper and free of branch issues.
     */

    /** Paraboloidal (Melro) failure measure of a stress tensor: phiBar = 1 on the surface. */

    /// Rice-Tracey / Smith exponent in the SWDFM driver. A micromechanical constant, NOT fitted.
    inline const static double swdfmExponent = 1.3;

    /// eta beyond these bounds is clamped so that epsF_eff stays bounded away from 0 and infinity.
    inline const static double etaMin = -0.5;
    inline const static double etaMax = 1.2;

    /** Damage variable omega and its derivatives wrt the local / nonlocal alphaP.
     *
     * Unscaled (c_eta == 0), IDENTICAL to the original formulation:
     *
     *     omega = 1 - exp( -alphaP_weighted / epsF )
     *
     * Triaxiality-scaled (c_eta != 0), the damage rate is scaled by the stress state and the
     * driver is accumulated INCREMENTALLY, because eta evolves during loading:
     *
     *     epsF_eff( eta ) = epsF * exp( -c_eta * clamp( eta ) )
     *     D  = D_old + max( 0, alphaP_weighted - alphaPBar_old ) / epsF_eff
     *     omega = 1 - exp( -D )
     *
     * The functional form of omega is unchanged; only the rate at which its driver accumulates
     * depends on the stress state. Accumulation makes damage irreversible.
     *
     * Returns { omega, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal, D_new, alphaPBar_new }.
     */
    std::tuple< double, double, double, double, double > computeOmega( const double     alphaP_local,
                                                                       const double     alphaP_nonlocal,
                                                                       const Tensor33d& tau_eff,
                                                                       const double     D_old,
                                                                       const double     alphaPBar_old )
    {
      const double alphaP_weighted = alphaP_nonlocal * m + alphaP_local * ( 1 - m );

      const double dAlphaP_weighted_dAlphaP_local    = 1 - m;
      const double dAlphaP_weighted_dAlphaP_nonlocal = m;

      if ( alphaP_weighted < 0.0 ) {
        return { 0.0, 0.0, 0.0, D_old, alphaPBar_old };
      }

      const double alphaPBar_new = alphaP_weighted;

      // ---- (1) SOFTENING variable: the original law, ALWAYS active and already calibrated.
      //          This is what reproduces the pre-peak response; removing it before initiation
      //          would lose the calibrated peak forces.
      //          Ds_inf BOUNDS it (see the declaration): Nguyen's Ds saturates and CANNOT fail the
      //          material; only the failure variable may. Ds_inf = 1 is the previous behaviour.
      //          The tail may be TRUNCATED (epsFCut > 0): see the declaration of epsFCut. Both the
      //          value and the derivative come from softening() so they cannot disagree.
      double omega_s, dOmega_s;
      softening( alphaP_weighted, omega_s, dOmega_s );

      // ---- (2) FAILURE variable: zero until the SWDFM driver reaches D = 1, then steep.
      //          Two-variable structure after Nguyen, Lani, Pardoen, Morelle & Noels,
      //          Int. J. Solids Struct. 96 (2016), which separates gradual softening from the
      //          final failure stage.
      double omega_f  = 0.0;
      double dOmega_f = 0.0;
      double D        = D_old;

      if ( cSW != 0.0 ) {
        const double T = std::min( std::max( triaxiality( tau_eff ), etaMin ), etaMax );

        // When the driver is already the dilatant plastic volume, the Rice-Tracey factor would
        // DOUBLE-COUNT the pressure sensitivity (exp(1.3T) dAlphaP is itself a void-growth
        // proxy), so it is switched off and the stress-state dependence comes from alphaD alone.
        // The driver increment is needed BEFORE g, because g may depend on the RATE at which the
        // driver accumulates (see swdfmLogG). dTCurrent is set at the top of computeStress.
        const double dAlphaPBar = std::max( alphaP_weighted - alphaPBar_old, 0.0 );
        const double kdot       = dTCurrent > 0.0 ? dAlphaPBar / dTCurrent : 0.0;

        // Rice-Tracey weight with the exponent E(T) * w(kdot) (see the declaration above); the
        // compression term keeps its mirrored form exp(-E)/bSW.
        double       kdotDLogGDKdot = 0.0;
        const double logG           = swdfmLogG( T, kdot, &kdotDLogGDKdot );

        double g = exp( logG ) - exp( -logG ) / bSW;
        g        = std::max( g, 0.0 ); // damage is irreversible under monotonic loading

        D = D_old + cSW * g * dAlphaPBar;

        if ( D > 1.0 ) {
          omega_f = 1.0 - exp( -( D - 1.0 ) / dF );
          // dD/dAlphaPBar carries BOTH the direct term and the rate sensitivity of g:
          //   dD/dk = cSW ( g + dk * dg/dkdot / dt ),  dk/(kdot dt) == 1
          //         = cSW ( g - s logG * lode * ( e^logG + e^-logG / bSW ) )
          // The kdot cancels exactly, so no extra division and no small-kdot blow-up.
          double dD_dAlphaPBar = cSW * g;
          if ( kdotDLogGDKdot != 0.0 )
            dD_dAlphaPBar += cSW * kdotDLogGDKdot * ( exp( logG ) + exp( -logG ) / bSW );
          dOmega_f = exp( -( D - 1.0 ) / dF ) / dF * dD_dAlphaPBar;
        }
      }

      const double omega                   = 1.0 - ( 1.0 - omega_s ) * ( 1.0 - omega_f );
      const double dOmega_dAlphaP_weigthed = ( 1.0 - omega_f ) * dOmega_s + ( 1.0 - omega_s ) * dOmega_f;

      // T, zeta AND the rate factor w are held fixed in the tangent (the dOmega/dT * dT/dTau and
      // dOmega/dw * dw/dAlphaPBar contributions are deliberately omitted), so the structure of the
      // tangent is unchanged. Omitting w is the same approximation already made for T and zeta and
      // is mild for the same reason: w varies by a factor 1.7 over four decades of rate, so its
      // derivative is small next to dOmega/dAlphaPBar itself.
      const double dOmega_dAlphaP_local    = dOmega_dAlphaP_weigthed * dAlphaP_weighted_dAlphaP_local;
      const double dOmega_dAlphaP_nonlocal = dOmega_dAlphaP_weigthed * dAlphaP_weighted_dAlphaP_nonlocal;

      if ( omega > omegaMax ) {
        return { omegaMax, 0.0, 0.0, D, alphaPBar_new };
      }
      else {
        return { omega, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal, D, alphaPBar_new };
      }
    }

    // ------------------------------------------------------------
    // Viscoelasticity (generalized Maxwell / Prony)
    // ------------------------------------------------------------

    /** Relax the instantaneous hyperelastic PK2 stress through the generalized Maxwell chain.
     *
     * The deviatoric and volumetric parts are relaxed independently, each by its own set of
     * relative moduli but sharing the relaxation times, so that shear creep and bulk creep can
     * differ. Both calls go through Marmot's own
     * ContinuumMechanics::FiniteStrain::Viscoelasticity::evaluateGeneralizedMaxwellModel, i.e.
     * the recurrence Q_i^{n+1} = exp(-dt/tau_i) Q_i^n + gamma_i (1-exp(-dt/tau_i))/(dt/tau_i) dS
     * and the matching closed-form tangent scaling are NOT reimplemented here.
     *
     * The volumetric branch stresses stay spherical for all time (the driver is spherical and
     * the recurrence is linear), so only their scalar magnitude is persisted; the 3x3 form is
     * rebuilt on entry and collapsed again on exit.
     *
     * @param PK2_0      instantaneous (unrelaxed) PK2 stress of the Pence-Gou potential
     * @param dPK2_0_dCe its derivative wrt the elastic right Cauchy-Green tensor
     * @param commit     false -> the branch states are advanced on a SCRATCH copy, leaving the
     *                            converged history untouched, so that this may be called as
     *                            often as the return-map Newton iteration needs;
     *                   true  -> the advanced branch states and the new reference stress are
     *                            written back to stateVars. Call exactly ONCE per increment,
     *                            with the converged elastic deformation.
     * @return { PK2 including relaxation, its derivative wrt Ce }
     */
    std::tuple< Tensor33d, Tensor3333d > applyViscoelasticity( const Tensor33d&   PK2_0,
                                                               const Tensor3333d& dPK2_0_dCe,
                                                               const bool         commit )
    {
      // no Maxwell branches -> bit-for-bit the original hyperelastic-viscoplastic model
      if ( nMaxwell == 0 )
        return { PK2_0, dPK2_0_dCe };

      using namespace ContinuumMechanics::FiniteStrain::Viscoelasticity;

      // ---- split the instantaneous stress and its tangent into volumetric / deviatoric parts.
      //      d(PK2vol)_ij/dCe_kl = 1/3 delta_ij d(tr PK2)/dCe_kl; the contraction below yields
      //      d(tr PK2)/dCe because 2 d2Psi/dCdC is major symmetric.
      const double    p0       = trace( PK2_0 ) / 3.0;
      const Tensor33d PK2Vol_0 = p0 * Spatial3D::I;
      const Tensor33d PK2Dev_0 = PK2_0 - PK2Vol_0;

      const Tensor3333d dPK2Vol_dCe = ( 1.0 / 3.0 ) *
                                      Fastor::outer( Spatial3D::I,
                                                     Tensor33d( einsum< ijkl, kl >( dPK2_0_dCe, Spatial3D::I ) ) );
      const Tensor3333d dPK2Dev_dCe = dPK2_0_dCe - dPK2Vol_dCe;

      // ---- the recurrence is driven by the increment since the last CONVERGED increment
      const Tensor33d PK2RefOld( stateVars->PK2Ref );
      const double    pRefOld = trace( PK2RefOld ) / 3.0;
      const Tensor33d dPK2Vol = ( p0 - pRefOld ) * Spatial3D::I;
      const Tensor33d dPK2Dev = PK2Dev_0 - ( PK2RefOld - pRefOld * Spatial3D::I );

      // ---- branch states. The helper writes through its pointer, so an uncommitted evaluation
      //      is given a scratch buffer instead of the persistent state.
      std::array< double, nMaxwellMax * 9 > devState{};
      std::array< double, nMaxwellMax * 9 > volState{};
      std::memcpy( devState.data(), stateVars->veDev, nMaxwell * 9 * sizeof( double ) );
      for ( int i = 0; i < nMaxwell; i++ )
        for ( int d = 0; d < 3; d++ )
          volState[i * 9 + d * 4] = stateVars->veVol[i]; // rebuild q_i * I

      Tensor33d   PK2Dev = PK2Dev_0, PK2Vol = PK2Vol_0;
      Tensor3333d tanDev = dPK2Dev_dCe, tanVol = dPK2Vol_dCe;

      evaluateGeneralizedMaxwellModel( PK2Dev, tanDev, dPK2Dev, dTCurrent, maxwellDev, devState.data() );
      evaluateGeneralizedMaxwellModel( PK2Vol, tanVol, dPK2Vol, dTCurrent, maxwellVol, volState.data() );

      if ( commit ) {
        std::memcpy( stateVars->veDev, devState.data(), nMaxwell * 9 * sizeof( double ) );
        for ( int i = 0; i < nMaxwell; i++ )
          stateVars->veVol[i] = volState[i * 9]; // spherical -> keep the magnitude only
        std::memcpy( stateVars->PK2Ref.data(), PK2_0.data(), 9 * sizeof( double ) );
      }

      return { Tensor33d( PK2Dev + PK2Vol ), Tensor3333d( tanDev + tanVol ) };
    }

    // ------------------------------------------------------------
    // Viscoplasticity functions
    // ------------------------------------------------------------

    std::tuple< double, Tensor33d, double, Tensor33d, Tensor33d > yieldFunction( const Tensor33d& Fe,
                                                                                 const double     betaP,
                                                                                 const double     J )
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
                dh_dMandel ) = yieldFunctionFromStress( mandelStress, betaP, J );
      dh_dFe                 = einsum< mn, mnij, to_ij >( dh_dMandel, dMandel_dFe );
      df_dFe                 = einsum< mn, mnij, to_ij >( df_dMandel, dMandel_dFe );

      return { f, df_dFe, df_dBetaP, dg_dMandel, dh_dFe };
    }

    std::tuple< double > yieldFunctionNominal( const Tensor33d& Fe,
                                               const double     betaP,
                                               const double     omega,
                                               const double     J )
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

      const double f = ( A * I1 + B ) / ( 2.0 * eta ) / J - betaP;

      return { f };
    }

    std::tuple< double, Tensor33d, double, Tensor33d, Tensor3333d, double, Tensor33d > yieldFunctionFromStress(
      const Tensor33d& mandelStress,
      const double     betaP,
      const double     J )
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

      const double f = ( A * I1 + B ) / ( 2.0 * eta ) / J - betaP;

      const double phiI1      = A * ( 1.0 + A * I1 / B ) / ( 2.0 * eta );
      const double phiJ2      = 3.0 / B;
      Tensor33d    df_dMandel = ( phiI1 * Spatial3D::I + phiJ2 * dev ) / J;

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

    bool isYielding( const Tensor33d& Fe, const double betaP, const double omega, const double J )
    {
      double    f, df_dBetaP;
      Tensor33d df_dFe, dg_dMandel, dh_dFe;
      std::tie( f, df_dFe, df_dBetaP, dg_dMandel, dh_dFe ) = yieldFunction( Fe, betaP, J );
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
      // Viscoelastic relaxation happens HERE, upstream of the return map, so that the yield
      // function sees the relaxed stress. commit = false: this is called inside the Newton loop.
      auto [PK2,
            dPK2_dCe] = applyViscoelasticity( Tensor33d( 2.0 * dPsi_dCe ), Tensor3333d( 2.0 * d2Psi_dCedCe ), false );

      const Tensor33d mandel = Ce % PK2;
      dMandel_dCe            = einsum< Ii, iJKL, to_IJKL >( Ce, dPK2_dCe ) +
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

      auto [PK2,
            dPK2_dCe] = applyViscoelasticity( Tensor33d( 2.0 * dPsi_dCe ), Tensor3333d( 2.0 * d2Psi_dCedCe ), false );

      const Tensor33d mandelN = ( Ce % PK2 ) * ( 1.0 - omega );

      dMandelN_dCe = ( einsum< Ii, iJKL, to_IJKL >( Ce, dPK2_dCe ) +
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
                                                                                    const double           dt,
                                                                                    const double           J )
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
                dh_dMandel ) = yieldFunctionFromStress( mandelStress, betaP, J );

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

      std::tie( f, df_dFe, df_dBetaP, dg_dMandel, dh_dFe ) = yieldFunction( Fe, betaP, J );

      double beta_min           = 1e-12;
      double sgn_beta           = ( betaP >= 0 ) ? 1.0 : -1.0;
      double betaP_cap          = sgn_beta * std::max( std::abs( betaP ), beta_min );
      double dBetaP_dAlphaP_cap = ( std::abs( betaP ) > beta_min ) ? sgn_beta * dBetaP_dAlphaP : 0.0;

      double ratio = Math::macauly( f + betaP_cap ) / betaP_cap;
      // double    D          = std::pow( ratio, 1.0 / n );
      // Tensor33d dD_dFe     = ( 1.0 / n ) * std::pow( ratio, ( 1.0 - n ) / n ) * df_dFe / betaP_cap;
      // double    dD_dalphaP = -( D / ( n * betaP_cap ) ) * dBetaP_dAlphaP_cap;

      const double    hmin  = 1e-8;
      const double    hsafe = std::max( h, hmin );
      const Tensor33d hZero( 0.0 );
      const Tensor33d dh_dMandel_safe = ( h > hmin ) ? dh_dMandel : hZero;
      const Tensor33d dh_dFe_safe     = einsum< mn, mnij, to_ij >( dh_dMandel_safe, dMandel_dFe );

      double q  = dLambda * hsafe / ( dt * eta_VP );
      double qn = std::pow( std::max( q, 1e-30 ), n );

      // std::cout << "D: " << D << std::endl;
      // std::cout << "f: " << f << std::endl;
      // std::cout << "betaP: " << betaP << std::endl;
      // Residual
      R.segment< 9 >( 0 ) += mV9d( Tensor33d( einsum< iJ, JK >( Fe, dFp ) ).data() );
      R( idxA ) += ( alphaP - dLambda * hsafe );
      R( idxF ) = ratio - qn;

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
      Tensor33d dRl_dFe = df_dFe / betaP_cap - ( n * qn / hsafe ) * dh_dFe_safe;

      dR_dX.block< 1, 9 >( idxF, 0 ) = mV9d( dRl_dFe.data() ).transpose();
      dR_dX( idxF, idxA )            = -dBetaP_dAlphaP_cap * ratio / betaP_cap;
      dR_dX( idxF, idxF )            = -n * qn / std::max( dLambda, 1e-30 );

      return { R, dR_dX };
    }
  };

} // namespace Marmot::Materials
