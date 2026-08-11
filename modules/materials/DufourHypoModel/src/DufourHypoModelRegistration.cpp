#include "Marmot/DufourHypoModel.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      DufourHypoModelRegistered = MarmotMaterialHypoElasticFactory::registerMaterial< class DufourHypoModel >(
        "DUFOURHYPOMODEL" );

  } // namespace Registration
} // namespace Marmot::Materials
