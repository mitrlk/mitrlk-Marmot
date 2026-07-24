#include "Marmot/DufourModel.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool DufourModelRegistered = MarmotMaterialGradientEnhancedFiniteStrainFactory::
      registerMaterial< class DufourModel >( "DUFOURMODEL" );

  } // namespace Registration
} // namespace Marmot::Materials
