CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(src_GENERAL
    general/bitarray.cpp
    general/flags.cpp
    general/hashtabl.cpp
    general/optmem.cpp
    general/sort.cpp
    general/symbolta.cpp
    general/table.cpp)

SET(src_GEOM2D
    geom2d/genmesh2d.cpp
    geom2d/geom2dmesh.cpp
    geom2d/geometry2d.cpp)

SET(src_GPRIM
    gprim/adtree.cpp
    gprim/geom2d.cpp
    gprim/geom3d.cpp
    gprim/geomfuncs.cpp
    gprim/geomtest3d.cpp
    gprim/spline.cpp
    gprim/splinegeometry.cpp
    gprim/transform3d.cpp)

SET(src_INTERFACE
    interface/writegmsh2.cpp
    interface/writeuser.cpp)

SET(src_LINALG
    linalg/bfgs.cpp
    linalg/densemat.cpp
    linalg/linopt.cpp
    linalg/linsearch.cpp
    linalg/polynomial.cpp)

SET(src_MESHING
    meshing/adfront2.cpp
    meshing/bisect.cpp
    meshing/global.cpp
    meshing/improve2.cpp
    meshing/improve2gen.cpp
    meshing/localh.cpp
    meshing/meshclass.cpp
    meshing/meshing2.cpp
    meshing/meshtool.cpp
    meshing/meshtype.cpp
    meshing/msghandler.cpp
    meshing/netrule2.cpp
    meshing/parser2.cpp
    meshing/quadrls.cpp
    meshing/refine.cpp
    meshing/ruler2.cpp
    meshing/smoothing2.cpp
    meshing/topology.cpp
    meshing/triarls.cpp)

#######################
# Build library
#######################

ADD_LIBRARY(meshit STATIC
            ${src_GENERAL}
            ${src_GEOM2D}
            ${src_GPRIM}
            ${src_INTERFACE}
            ${src_LINALG}
            ${src_MESHING})

#######################
# Export header files
#######################

LIST(APPEND cmd_export_headers COMMAND ${CMAKE_COMMAND} -E make_directory ${INCLUDE_OUTPUT_PATH}/meshit/)

FOREACH(hdr_dir . general geom2d gprim interface linalg meshing)
    MESSAGE("Exporting headers for ${hdr_dir}...")
    FILE(GLOB headers_list "${CMAKE_CURRENT_SOURCE_DIR}/${hdr_dir}/*.hpp")
    LIST(APPEND cmd_export_headers COMMAND ${CMAKE_COMMAND} -E make_directory ${INCLUDE_OUTPUT_PATH}/meshit/${hdr_dir})
    FOREACH(hdr_file ${headers_list})
        LIST(APPEND cmd_export_headers COMMAND ${CMAKE_COMMAND} -E copy ${hdr_file} ${INCLUDE_OUTPUT_PATH}/meshit/${hdr_dir})
    ENDFOREACH(hdr_file)
ENDFOREACH(hdr_dir)
ADD_CUSTOM_COMMAND(TARGET meshit PRE_LINK ${cmd_export_headers})

#######################
# Export variables
#######################

SET(MESHIT_LINK_LIBRARIES meshit PARENT_SCOPE)
