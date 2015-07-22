PROJECT( distBlocks )

cmake_minimum_required( VERSION 2.6 )

FIND_LIBRARY( VISTA NAMES libviaio.a PATHS /../../lib/ )

MESSAGE( STATUS "Libraries: ${VISTA}" )

FIND_PACKAGE( OpenMP )
SET( CMAKE_C_FLAGS
    "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}"
)
SET( CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}"
)

INCLUDE_DIRECTORIES(
    /usr/include
)

LINK_LIBRARIES(
    ${VISTA}
)

SET( SRCS
    main.cpp

    # classes
    coordinate.cpp
    tnode.cpp

    # functions
    output.cpp
    get_roi.cpp
    get_names.cpp
    do_Vista.cpp
    tract_dist.cpp
    do_dist_blocks.cpp

)

ADD_EXECUTABLE( distBlocks
    ${SRCS}
)

