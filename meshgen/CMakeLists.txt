CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

PROJECT( meshgen )

SET(src_CSG
#	csg/algprim.cpp
#	csg/brick.cpp
#	csg/bspline2d.cpp
#	csg/csgeom.cpp
#	csg/csgparser.cpp
#	csg/csgpkg.cpp
#	csg/curve2d.cpp
#	csg/edgeflw.cpp
#	csg/explicitcurve2d.cpp
#	csg/extrusion.cpp
#	csg/gencyl.cpp
#	csg/genmesh.cpp
#	csg/identify.cpp
#	csg/manifold.cpp
#	csg/meshsurf.cpp
#	csg/polyhedra.cpp
#	csg/revolution.cpp
#	csg/singularref.cpp
#	csg/solid.cpp
#	csg/specpoin.cpp
#	csg/spline3d.cpp
#	csg/surface.cpp
#	csg/triapprox.cpp
#	csg/vscsg.cpp
)
SET(src_GENERAL
	general/array.cpp
	general/bitarray.cpp
	general/dynamicmem.cpp
	general/flags.cpp
	general/gzstream.cpp
	general/hashtabl.cpp
	general/mystring.cpp
	general/ngexception.cpp
	general/optmem.cpp
	general/seti.cpp
	general/sort.cpp
	general/spbita2d.cpp
	general/symbolta.cpp
	general/table.cpp
)
SET(src_GEOM2D
	geom2d/genmesh2d.cpp
	geom2d/geom2dmesh.cpp
	geom2d/geom2dpkg.cpp
	geom2d/geometry2d.cpp
	geom2d/vsgeom2d.cpp
)
SET(src_GPRIM
	gprim/adtree.cpp
	gprim/geom2d.cpp
	gprim/geom3d.cpp
	gprim/geomfuncs.cpp
	gprim/geomtest3d.cpp
	gprim/spline.cpp
	gprim/splinegeometry.cpp
	gprim/transform3d.cpp
)
SET(src_INTERFACE
#	interface/nginterface.cpp
#	interface/nginterface_v2.cpp
	interface/read_fnf_mesh.cpp
	interface/readtetmesh.cpp
	interface/readuser.cpp
	interface/writeabaqus.cpp
	interface/writediffpack.cpp
	interface/writedolfin.cpp
	interface/writeelmer.cpp
	interface/writefeap.cpp
	interface/writefluent.cpp
	interface/writegmsh2.cpp
	interface/writegmsh.cpp
	interface/writejcm.cpp
	interface/writeOpenFOAM15x.cpp
	interface/writepermas.cpp
	interface/writetecplot.cpp
	interface/writetet.cpp
	interface/writetochnog.cpp
	interface/writeuser.cpp
	interface/wuchemnitz.cpp
)
SET(src_LINALG
	linalg/bfgs.cpp
	linalg/densemat.cpp
	linalg/linopt.cpp
	linalg/linsearch.cpp
	linalg/polynomial.cpp
)
SET(src_MESHING
	meshing/adfront2.cpp
	meshing/adfront3.cpp
	meshing/basegeom.cpp
	meshing/bcfunctions.cpp
	meshing/bisect.cpp
	meshing/clusters.cpp
	meshing/curvedelems.cpp
	meshing/delaunay2d.cpp
	meshing/delaunay.cpp
	meshing/geomsearch.cpp
	meshing/global.cpp
	meshing/hprefinement.cpp
	meshing/improve2.cpp
	meshing/improve2gen.cpp
	meshing/improve3.cpp
	meshing/localh.cpp
	meshing/meshclass.cpp
#	meshing/meshfunc2d.cpp
	meshing/meshfunc.cpp
	meshing/meshing2.cpp
	meshing/meshing3.cpp
	meshing/meshtool.cpp
	meshing/meshtype.cpp
	meshing/msghandler.cpp
	meshing/netrule2.cpp
	meshing/netrule3.cpp
	meshing/parallelmesh.cpp
	meshing/paralleltop.cpp
	meshing/parser2.cpp
	meshing/parser3.cpp
	meshing/prism2rls.cpp
	meshing/pyramid2rls.cpp
	meshing/pyramidrls.cpp
	meshing/quadrls.cpp
	meshing/refine.cpp
	meshing/ruler2.cpp
	meshing/ruler3.cpp
	meshing/secondorder.cpp
	meshing/smoothing2.5.cpp
	meshing/smoothing2.cpp
	meshing/smoothing3.cpp
	meshing/specials.cpp
	meshing/tetrarls.cpp
	meshing/topology.cpp
	meshing/triarls.cpp
	meshing/validate.cpp
#	meshing/zrefine.cpp
)
SET (src_OCC
#	occ/occconstruction.cpp
#	occ/occgenmesh.cpp
#	occ/occgeom.cpp
#	occ/occmeshsurf.cpp
#	occ/occpkg.cpp
#	occ/vsocc.cpp
)
SET(src_STLGEOM
#	stlgeom/meshstlsurface.cpp
#	stlgeom/stlgeomchart.cpp
#	stlgeom/stlgeom.cpp
#	stlgeom/stlgeommesh.cpp
#	stlgeom/stlline.cpp
#	stlgeom/stlpkg.cpp
#	stlgeom/stltool.cpp
#	stlgeom/stltopology.cpp
#	stlgeom/vsstl.cpp
)
SET(src_VISUALIZATION
#	visualization/importsolution.cpp
#	visualization/meshdoc.cpp
#	visualization/mvdraw.cpp
#	visualization/vsfieldlines.cpp
#	visualization/vsmesh.cpp
#	visualization/vssolution.cpp
)
ADD_LIBRARY(meshgen STATIC
	${src_CSG}
	${src_GENERAL}
	${src_GEOM2D}
	${src_GPRIM}
	${src_LINALG}
	${src_MESHING}
	${src_OCC}
	${src_STLGEOM}
	${src_VISUALIZATION}
)
TARGET_LINK_LIBRARIES(meshgen)
