# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for GEDisplacementFiniteStrainULElement
add_marmot_test("TestGEDisplacementFiniteStrainULElement" "${CURR_TEST_SOURCE_DIR}/TestGEDisplacementFiniteStrainULElement.cpp")
