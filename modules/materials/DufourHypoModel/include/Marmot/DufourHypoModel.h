/* ---------------------------------------------------------------------
 *  DufourHypoModel
 *
 *  A 1:1 HYPOELASTIC port of Dufour's mesoscopic adhesive model, written to
 *  answer one question: does his published Arcan agreement come from the
 *  FORMULATION rather than from the parameters?
 *
 *  Our in-tree DufourModel is a finite-strain, gradient-enhanced, Mandel-stress
 *  formulation. Dufour's own model (thesis eq. 2.9 / IJAA eq. 1) is HYPOELASTIC:
 *  an additive split of the rate of deformation with a logarithmic corotational
 *  rate, whose working stress IS the Cauchy stress. Consequences, and the four
 *  and only four differences from DufourModel:
 *
 *    1. Cauchy stress integrated in rate form, not Mandel from Fe.
 *    2. NO 1/J push-forward in the yield function. In DufourModel the `/J` is a
 *       REQUIRED Mandel->Cauchy conversion (f is degree-1 homogeneous); here the
 *       working stress is already Cauchy so no factor is needed, exactly as his
 *       eq. 2.26 / paper eq. 14 contain no Jacobian anywhere.
 *    3. Small-strain kinematics.
 *    4. LOCAL damage, D = D(kappa), not the implicit-gradient nonlocal kappa_bar.
 *       He runs local damage; that is also why he can use 2 elements through the
 *       0.3 mm bondline while our ld = 0.12 forces 5.
 *
 *  Everything else -- the yield surface, the hardening law, the non-associated
 *  flow potential, the rate law and the damage law -- is character-for-character
 *  the same as DufourModel, deliberately, so that a single-element comparison
 *  between the two isolates the formulation and nothing else. At small strain the
 *  two MUST agree to O(strain^2); that is the verification of this file.
 *
 *  Thesis equations implemented:
 *    (2.26)  f = [ (eta-1) I1 + sqrt( (eta-1)^2 I1^2 + 12 eta J2 ) ] / (2 eta)
 *                - ( ft + R(kappa) ) * ( kappaDot / kappaDot0 )^n
 *    (2.28)  R(kappa) = Q1 k exp(-b1 k) + Q2 (1 - exp(-b2 k)) + b3 k^3 + b4 k^2 + b5 k
 *    (2.29)  F = sqrt( 3 J2 + alphaPlus <p>^2 + alphaMinus <-p>^2 )      (non-associated)
 *    (2.31)  alphaPlus  = (9/2) (1 - 2 nuP_plus ) / (1 + nuP_plus )
 *    (2.32)  alphaMinus = (9/2) (1 - 2 nuP_minus) / (1 + nuP_minus)
 *    (2.42)  D = 1 - exp( -kappa / kappaC ),  capped at omegaMax
 *
 *  Damage acts through the EFFECTIVE stress: the return map runs on
 *  sigmaEff = sigma / (1 - D) and the reported stress is (1 - D) sigmaEff. Both f
 *  and F are degree-1 homogeneous in stress, so this is identical to carrying the
 *  1/(1-D) inside them as printed in eqs 2.26 and 2.29.
 *
 *  CARD -- identical 20-entry layout to DufourModel, so every existing .json card
 *  works unchanged:
 *    [0] K   [1] G   [2] ft  [3] fc  [4] Q1  [5] Q2  [6] b1  [7] b2  [8] b3
 *    [9] b4  [10] b5 [11] eta_VP (= kappaDot0)  [12] n
 *    [13] nuP_plus  [14] nuP_minus  [15] epsF (= kappaC)  [16] omegaMax
 *    [17] ld  [18] m  [19] density
 *  ld and m are READ BUT UNUSED -- this model is local by construction. They are
 *  kept so the cards stay interchangeable.
 * ---------------------------------------------------------------------
 */

#pragma once
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotPronySeries.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Marmot;

namespace Marmot::Materials {

  class DufourHypoModel : public MarmotMaterialHypoElastic {
  public:
    using MarmotMaterialHypoElastic::MarmotMaterialHypoElastic;

    // elastic
    const double K, G;
    // yield surface + hardening
    const double ft, fc, Q1, Q2, b1, b2, b3, b4, b5;
    // rate law
    const double eta_VP, n;
    // non-associated flow (plastic Poisson ratios)
    const double nuP_plus, nuP_minus;
    // damage
    const double epsF, omegaMax;
    // read for card compatibility with DufourModel, UNUSED here (local model)
    const double ld, m;

    /* OPTIONAL card entries 20/21 -- EVOLVING PRESSURE SENSITIVITY, eta = fc/ft -> eta(kappa).
     *
     *   eta(kappa) = etaInf + ( eta0 - etaInf ) exp( -kappa / etaKappa ),   etaInf = etaRatio*eta0
     *
     * etaRatio = 1 (the default when the card has only 20 entries) gives eta == eta0 EXACTLY, so
     * every existing 20-entry deck is bit-for-bit unchanged.  The construction is the standard one
     * from pressure-sensitive plasticity: Grassl & Jirasek, IJSS 43 (2006) 7166-7196, sec. 3.1.3 --
     * "It controls the size AND SHAPE of the yield surface".  There the shape variable is qh(kappa_p);
     * here it is eta(kappa).  Dufour fixes eta = 1.8 by two tangent constructions at INITIAL yield
     * and then holds it to kappa = 0.63; this lets it relax. */
    const double etaRatio, etaKappa;

    /* OPTIONAL card entries 22/23 -- EVOLVING PLASTIC POISSON RATIO, nuP_plus -> nuP_plus(kappa).
     *
     *   nuP_plus(kappa) = nuP_plus + ( nuP0 - nuP_plus ) exp( -kappa / kappaNu )
     *
     * nuP0 = 0 (the default when the card is shorter) gives nuP_plus == the card value EXACTLY, so
     * every existing deck is unchanged. This is NOT a new mechanism: Dufour's own sect. 2.4.2.3
     * says the MEASURED nu_p "decreases ... to a value of 0.3", and figure 2.10a shows the decay.
     * He took the ASYMPTOTE as a constant and discarded the small-strain end as DIC noise. Frozen
     * from that figure (digitize_fig210.py / fit_nup_kappa.py): nuP0 = 0.4130, kappaNu = 0.01750.
     *
     * IMPLEMENTATION: evaluated at the CONVERGED kappa OF THE PREVIOUS INCREMENT (lagged), so the
     * local Newton sees it as a constant and the analytic Jacobian is untouched. Error is O(dKappa)
     * and must be checked by halving maxInc. */
    const double nuP0, kappaNu;

    /* OPTIONAL card entry 28 onwards -- GENERALIZED-MAXWELL (PRONY) VISCOELASTICITY.
     *
     * SAME CARD SLOTS AS DufourModel, deliberately, so ONE card drives both models:
     *     [27]          nMaxwell        number of branches (0 = off)
     *     [28 + 3 i]    gammaG_i        deviatoric relative modulus of branch i
     *     [29 + 3 i]    gammaK_i        volumetric relative modulus of branch i
     *     [30 + 3 i]    tau_i           relaxation time of branch i [s]
     *
     * gamma_i are RELATIVE moduli (G_i/G_0), and the equilibrium branch carries the remainder
     * ( 1 - sum_i gamma_i ). nMaxwell = 0, or a card shorter than 28 entries, reproduces the
     * previous purely elastoplastic model BIT FOR BIT.
     *
     * WHERE IT ACTS: on the ELASTIC PREDICTOR, i.e. upstream of the return map, so the yield
     * function sees the relaxed stress. This mirrors DufourModel, where relaxation acts on the
     * PK2 stress of the elastic potential before the viscoplastic update. The algorithmic
     * viscoelastic tangent C_ve replaces Cel EVERYWHERE in the local Newton and in the consistent
     * tangent, so the quadratic convergence of the 7x7 system is preserved.
     *
     * The branch state is advanced with the ELASTIC part of the strain increment,
     * dEps - dLambda * m, not the total -- the Maxwell elements sit in the elastic spring.
     *
     * NOTE this model is SMALL STRAIN, so it uses Marmot's small-strain
     * ContinuumMechanics::Viscoelasticity::PronySeries helper, NOT the finite-strain one that
     * DufourModel uses. The two are not expected to agree beyond small strain -- that difference
     * is precisely what this twin exists to measure. */
    const int                                                            nMaxwell;
    Marmot::ContinuumMechanics::Viscoelasticity::PronySeries::Properties pronyProps;

    /// isotropic Voigt stiffness from a (K, G) pair. Built directly rather than through (E, nu)
    /// because a branch may carry gammaK = 0, which makes E = 0 and nu = -1 and is not
    /// representable in (E, nu) form.
    static Marmot::Matrix6d isotropicC( const double Kb, const double Gb )
    {
      const double     lam = Kb - 2.0 * Gb / 3.0;
      Marmot::Matrix6d C   = Marmot::Matrix6d::Zero();
      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j )
          C( i, j ) = lam;
        C( i, i ) += 2.0 * Gb;
        C( 3 + i, 3 + i ) = Gb;
      }
      return C;
    }

    // derived
    const double E, nu, eta, alphaPlus, alphaMinus;
    const double eta0, etaInf;

    /// alphaPlus evaluated at a given kappa. nuP0 <= 0 -> the constant card value, exactly.
    inline double alphaPlusOf( const double kappa ) const
    {
      if ( nuP0 <= 0.0 )
        return alphaPlus;
      const double nuP = nuP_plus + ( nuP0 - nuP_plus ) * std::exp( -kappa / kappaNu );
      return 9.0 * ( 1.0 - 2.0 * nuP ) / ( 1.0 + nuP ) / 2.0;
    }

    /// eta(kappa) and its derivative. etaRatio == 1 -> ( eta0, 0 ) exactly.
    inline double etaOf( const double kappa ) const
    {
      return etaInf + ( eta0 - etaInf ) * std::exp( -kappa / etaKappa );
    }
    inline double dEtaOf( const double kappa ) const
    {
      return -( eta0 - etaInf ) / etaKappa * std::exp( -kappa / etaKappa );
    }

    DufourHypoModel( const double* materialProperties, int nMaterialProperties, int materialNumber );

    double getDensity( const double* stateVars ) const override;

    /** Fully analytic: coupled 7-unknown local Newton with an analytic Jacobian, and the
     *  consistent tangent from the implicit function theorem. See the .cpp header comment for
     *  why the three earlier AD-based local solvers were abandoned. */
    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStress_dStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override;

  protected:
    /** Raghava part of eq 2.26 (no /J), its stress gradient, AND dPhi/dkappa.
     *  dPhi/dkappa is nonzero only through eta(kappa); it is identically 0 when etaRatio == 1,
     *  which is what makes the extension free for every existing deck. */
    void phiAndGrad( const Marmot::Vector6d& s,
                     const double            kappa,
                     double&                 Phi,
                     Marmot::Vector6d&       dPhi,
                     double&                 dPhidKappa ) const;

    /// Non-associated potential eq 2.29: m = dg/dSigma, M = dm/dSigma, h = sqrt(2/3 m:m), dh/dSigma.
    /** Non-associated potential eq 2.29.  `aPlus` is alphaPlus, passed in so it can be evaluated
     *  at a LAGGED kappa when nuP_plus(kappa) is active; pass `alphaPlus` for the constant case. */
    void flowAndGrad( const Marmot::Vector6d& s,
                      const double            aPlus,
                      Marmot::Vector6d&       mOut,
                      Marmot::Matrix6d&       MOut,
                      double&                 hOut,
                      Marmot::Vector6d&       dhOut ) const;

    /// Isotropic hardening, thesis eq. 2.28.
    template < typename T >
    T R( const T kappa ) const
    {
      return Q1 * kappa * exp( -b1 * kappa ) + Q2 * ( 1.0 - exp( -b2 * kappa ) ) + b3 * kappa * kappa * kappa +
             b4 * kappa * kappa + b5 * kappa;
    }

    /// Current strength beta_p = ft + R(kappa).
    template < typename T >
    T betaP( const T kappa ) const
    {
      return ft + R( kappa );
    }

    /// Damage, thesis eq. 2.42, capped at omegaMax.
    double omegaOf( const double kappa ) const
    {
      const double D = 1.0 - std::exp( -kappa / epsF );
      return std::min( D, omegaMax );
    }
  };
} // namespace Marmot::Materials
