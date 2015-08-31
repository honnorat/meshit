CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

PROJECT( meshit )

SET(src_GENERAL
	src/general/bitarray.cpp
	src/general/flags.cpp
	src/general/hashtabl.cpp
	src/general/optmem.cpp
	src/general/symbolta.cpp
	src/general/table.cpp
)
SET(src_GEOM2D
	src/geom2d/genmesh2d.cpp
	src/geom2d/geom2dmesh.cpp
	src/geom2d/geometry2d.cpp
)
SET(src_GPRIM
	src/gprim/adtree.cpp
	src/gprim/geom2d.cpp
	src/gprim/geom3d.cpp
	src/gprim/geomfuncs.cpp
	src/gprim/geomtest3d.cpp
	src/gprim/spline.cpp
	src/gprim/splinegeometry.cpp
)
SET(src_IO
	src/io/writeuser.cpp
	src/io/writegmsh2.cpp
)
SET(src_LINALG
	src/linalg/bfgs.cpp
	src/linalg/densemat.cpp
	src/linalg/linopt.cpp
	src/linalg/linsearch.cpp
	src/linalg/polynomial.cpp
)
SET(src_MESHING
	src/meshing/adfront2.cpp
	src/meshing/clusters.cpp
	src/meshing/global.cpp
	src/meshing/improve2.cpp
	src/meshing/improve2gen.cpp
	src/meshing/localh.cpp
	src/meshing/meshclass.cpp
	src/meshing/meshing2.cpp
	src/meshing/meshtool.cpp
	src/meshing/meshtype.cpp
	src/meshing/netrule2.cpp
	src/meshing/optimize2d.cpp
	src/meshing/quadrls.cpp
	src/meshing/refinement.cpp
	src/meshing/ruler2.cpp
	src/meshing/smoothing2.cpp
	src/meshing/topology.cpp
	src/meshing/triarls.cpp
)

ADD_LIBRARY(meshit STATIC
	${src_GENERAL}
	${src_GEOM2D}
	${src_GPRIM}
	${src_IO}
	${src_LINALG}
	${src_MESHING}
)
TARGET_LINK_LIBRARIES(meshit)
