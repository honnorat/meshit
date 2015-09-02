CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

SET(MESHGEN_PATH ..)
INCLUDE_DIRECTORIES(${MESHGEN_PATH})

MACRO(meshit_add_test examplename)
  SET(targetname test_${examplename})
  SET(filename ${examplename}.cpp)
  ADD_EXECUTABLE(${targetname} ${filename})
  TARGET_LINK_LIBRARIES(${targetname} meshgen)
ENDMACRO(meshit_add_test)

FILE(GLOB tests "*.cpp")
FOREACH( test ${tests} )
	GET_FILENAME_COMPONENT(test_name ${test} NAME_WE)
	MESSAGE("found test: " ${test_name})
	meshit_add_test( ${test_name} )
ENDFOREACH()
