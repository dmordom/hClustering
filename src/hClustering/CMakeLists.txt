PROJECT( hClustering_1.0 )

cmake_minimum_required( VERSION 2.6 )

# guard against in-source builds
IF( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
  MESSAGE( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. (you may need to remove CMakeCache.txt" )
ENDIF()

# The build types
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE Release
         CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE
       )
ENDIF( NOT CMAKE_BUILD_TYPE )

# guard against typos in build-type strings
STRING( TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
IF( NOT cmake_build_type_tolower STREQUAL "debug" AND
    NOT cmake_build_type_tolower STREQUAL "release" AND
    NOT cmake_build_type_tolower STREQUAL "relwithdebinfo" AND
    NOT cmake_build_type_tolower STREQUAL "")
  MESSAGE( SEND_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo  and \"\" (case-insensitive).")
ENDIF()

# find and link lipsia library
FIND_LIBRARY( VISTA NAMES libviaio.a PATHS /../../lib/ )


MESSAGE( STATUS "Libraries: ${VISTA}" )
LINK_LIBRARIES( ${VISTA} )


# find and link boost packages
FIND_PACKAGE( Boost 1.39.0 REQUIRED program_options thread filesystem system regex )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
LINK_LIBRARIES ( ${Boost_LIBRARIES} )


FIND_PACKAGE( OpenMP )
SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

INCLUDE_DIRECTORIES( /usr/include ../ ../common ../../include )

SET( COMMON_SRCS
	../common/cnbTreeBuilder.cpp
	../common/compactTractChar.cpp
	../common/compactTract.cpp
	../common/distBlock.cpp
	../common/graphTreeBuilder.cpp
	../common/partitionMatcher.cpp
	../common/protoNode.cpp
	../common/randCnbTreeBuilder.cpp
	../common/treeComparer.cpp
	../common/treeManager.cpp
	../common/vistaManager.cpp
	../common/WFileParser.cpp
	../common/WHcoord.cpp
	../common/WHnode.cpp
	../common/WHtree.cpp
	../common/WHtreePartition.cpp
	../common/WHtreeProcesser.cpp
	../common/WStringUtils.cpp
)

SET( MAIN_SRCS
	basetracts.cpp
	buildctree.cpp
	buildgraphtree.cpp
	compact2full.cpp
	comparetrees.cpp
	cpcc.cpp
	fliptracts.cpp
	fliptree.cpp
	full2compact.cpp
	fulltract.cpp
	matchpartition.cpp
	pairdist.cpp
	partitiontree.cpp
	prunetree.cpp
	randtracts.cpp
	treetracts.cpp
)

foreach( binarysourcefile ${MAIN_SRCS} )
    string( REPLACE ".cpp" "" binary ${binarysourcefile} )
    add_executable( ${binary} ${binarysourcefile} ${COMMON_SRCS} )
endforeach( binarysourcefile ${MAIN_SRCS} )





