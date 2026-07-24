#include "Marmot/Marmot.h"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterial.h"
#include <cassert>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>

namespace MarmotLibrary {

  // MaterialFactory

  // Construct-on-first-use idiom: the maps live inside accessor functions so they are
  // guaranteed to be constructed before any registration static initializer uses them,
  // regardless of translation-unit initialization order (which is unspecified across TUs).
  std::unordered_map< std::string, int >& MarmotMaterialFactory::materialNameToCodeAssociation()
  {
    static std::unordered_map< std::string, int > theMap;
    return theMap;
  }
  std::unordered_map< int, MarmotMaterialFactory::materialFactoryFunction >&
    MarmotMaterialFactory::materialFactoryFunctionByCode()
  {
    static std::unordered_map< int, materialFactoryFunction > theMap;
    return theMap;
  }

  bool MarmotMaterialFactory::registerMaterial( int                     materialCode,
                                                const std::string&      materialName,
                                                materialFactoryFunction factoryFunction )
  {
    assert( materialNameToCodeAssociation().find( materialName ) == materialNameToCodeAssociation().end() );
    assert( materialFactoryFunctionByCode().find( materialCode ) == materialFactoryFunctionByCode().end() );

    materialNameToCodeAssociation()[materialName] = materialCode;
    materialFactoryFunctionByCode()[materialCode] = factoryFunction;

    return true;
  }

  int MarmotMaterialFactory::getMaterialCodeFromName( const std::string& materialName )
  {
    try {
      return materialNameToCodeAssociation().at( materialName );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( MakeString() << "Invalid material " << materialName << " requested!" );
    }
  }

  MarmotMaterial* MarmotMaterialFactory::createMaterial( int           materialCode,
                                                         const double* materialProperties,
                                                         int           nMaterialProperties,
                                                         int           materialNumber )
  {
    try {
      return materialFactoryFunctionByCode().at(
        materialCode )( materialProperties, nMaterialProperties, materialNumber );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( MakeString() << "Invalid material " << materialCode << " requested!" );
    }
  }

  // ElementFactory

  std::unordered_map< std::string, int >& MarmotElementFactory::elementNameToCodeAssociation()
  {
    static std::unordered_map< std::string, int > theMap;
    return theMap;
  }
  std::unordered_map< int, MarmotElementFactory::elementFactoryFunction >&
    MarmotElementFactory::elementFactoryFunctionByCode()
  {
    static std::unordered_map< int, elementFactoryFunction > theMap;
    return theMap;
  }

  bool MarmotElementFactory::registerElement( const std::string&     elementName,
                                              int                    elementCode,
                                              elementFactoryFunction factoryFunction )
  {
    assert( elementNameToCodeAssociation().find( elementName ) == elementNameToCodeAssociation().end() );
    assert( elementFactoryFunctionByCode().find( elementCode ) == elementFactoryFunctionByCode().end() );

    elementNameToCodeAssociation()[elementName] = elementCode;
    elementFactoryFunctionByCode()[elementCode] = factoryFunction;

    return true;
  }

  int MarmotElementFactory::getElementCodeFromName( const std::string& elementName )
  {
    try {
      return elementNameToCodeAssociation().at( elementName );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( MakeString() << "Invalid element " << elementName << " requested!" );
    }
  }

  MarmotElement* MarmotElementFactory::createElement( int elementCode, int elementNumber )
  {
    try {
      return elementFactoryFunctionByCode().at( elementCode )( elementNumber );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( MakeString() << "Invalid element " << elementCode << " requested!" );
    }
  }

} // namespace MarmotLibrary
