#ifndef FILE_MESHTOOL
#define FILE_MESHTOOL

#include "meshclass.hpp"

namespace meshit {
///
void MeshQuality2d (const Mesh & mesh);

///
void MeshQuality3d (const Mesh & mesh,
			   Array<int> * inclass = NULL);

///
void SaveEdges (const Mesh & mesh, 
		       const char * geomfile, 
		       double h, 
		       char * filename);

///
void SaveSurfaceMesh (const Mesh & mesh,
			     double h,
			     char * filename);
/*
///
extern void Save2DMesh (
         const Mesh & mesh2d,
	 const Array<class SplineSegment*> * splines,
         ostream & outfile);
*/

class Surface;
///
void SaveVolumeMesh (
         const Array<Point3d> & points,
         const Array<Element> & elements,
         const Array<Element> & volelements,
         const Array<Surface*> & surfaces,
         char * filename);

///
void SaveVolumeMesh (const Mesh & mesh, 
		     const class CSGeometry & geometry,
		     char * filename);

///
int CheckCode ();


///
double CalcTetBadness (const Point3d & p1, const Point3d & p2,
			      const Point3d & p3, const Point3d & p4, 
			      double h,
			      const MeshingParameters & mp);
///
double CalcTetBadnessGrad (const Point3d & p1, const Point3d & p2,
				  const Point3d & p3, const Point3d & p4, 
				  double h, int pi,
				  Vec<3> & grad,
				  const MeshingParameters & mp);


/** Calculates volume of an element.
  The volume of the tetrahedron el is computed
 */
// extern double CalcVolume (const Array<Point3d> & points,
//        const Element & el);  

/** The total volume of all elements is computed.
  This function calculates the volume of the mesh */
double CalcVolume (const Array<Point3d> & points, 
	const Array<Element> & elements);

///
int CheckSurfaceMesh (const Mesh & mesh);

///
int CheckSurfaceMesh2 (const Mesh & mesh);
///
int CheckMesh3D (const Mesh & mesh);
///
void RemoveProblem (Mesh & mesh, int domainnr);

}
#endif
