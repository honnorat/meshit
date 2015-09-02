#include <mystdlib.h>


#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>


namespace netgen
{
  Array<SpecialPoint> specpoints;
  static Array<MeshPoint> spoints;

#define TCL_OK 0
#define TCL_ERROR 1



  static void FindPoints (CSGeometry & geom, Mesh & mesh)
  {
    PrintMessage (1, "Start Findpoints");

    const char * savetask = multithread.task;
    multithread.task = "Find points";

    for (int i = 0; i < geom.GetNUserPoints(); i++)
      {
	mesh.AddPoint(geom.GetUserPoint (i));
	mesh.Points().Last().Singularity (geom.GetUserPointRefFactor(i));
	mesh.AddLockedPoint (PointIndex (i+1));
      }

    SpecialPointCalculation spc;

    spc.SetIdEps(geom.GetIdEps());

    if (spoints.size() == 0)
      spc.CalcSpecialPoints (geom, spoints);
    
    PrintMessage (2, "Analyze spec points");
    spc.AnalyzeSpecialPoints (geom, spoints, specpoints);
  
    PrintMessage (5, "done");

    std::cerr << specpoints.size() << " special points:" <<std::endl;
    for (int i = 0; i < specpoints.size(); i++)
      specpoints[i].Print std::cerr;

    /*
      for (int i = 1; i <= geom.identifications.Size(); i++)
      geom.identifications.Elem(i)->IdentifySpecialPoints (specpoints);
    */
    multithread.task = savetask;
  }






  static void FindEdges (CSGeometry & geom, Mesh & mesh, const bool setmeshsize = false)
  {
    EdgeCalculation ec (geom, specpoints);
    ec.SetIdEps(geom.GetIdEps());
    ec.Calc (mparam.maxh, mesh);

    for (int i = 0; i < geom.singedges.size(); i++)
      {
	geom.singedges[i]->FindPointsOnEdge (mesh);
	if(setmeshsize)
	  geom.singedges[i]->SetMeshSize(mesh,10.*geom.BoundingBox().Diam());
      }
    for (int i = 0; i < geom.singpoints.size(); i++)
      geom.singpoints[i]->FindPoints (mesh);

    for (int i = 1; i <= mesh.GetNSeg(); i++)
      {
	//std::cerr << "segment " << mesh.LineSegment(i) <<std::endl;
	int ok = 0;
	for (int k = 1; k <= mesh.GetNFD(); k++)
	  if (mesh.GetFaceDescriptor(k).SegmentFits (mesh.LineSegment(i)))
	    {
	      ok = k;
	      //std::cerr << "fits to " << k <<std::endl;
	    }

	if (!ok)
	  {
	    ok = mesh.AddFaceDescriptor (FaceDescriptor (mesh.LineSegment(i)));
	    //std::cerr << "did not find, now " << ok <<std::endl;
	  }

	//std::cerr << "change from " << mesh.LineSegment(i).si;
	mesh.LineSegment(i).si = ok;
	//std::cerr << " to " << mesh.LineSegment(i).si <<std::endl;
      }

    if (geom.identifications.size())
      {
	PrintMessage (3, "Find Identifications");
	for (int i = 0; i < geom.identifications.size(); i++)
	  {
	    geom.identifications[i]->IdentifyPoints (mesh);
	    //std::cerr << "identification " << i << " is " 
	    //	       << *geom.identifications[i] <<std::endl;
	    
	  }
	for (int i = 0; i < geom.identifications.size(); i++)
	  geom.identifications[i]->IdentifyFaces (mesh);
      }


    // find intersecting segments
    PrintMessage (3, "Check intersecting edges");
    
    Point3d pmin, pmax;
    mesh.GetBox (pmin, pmax);
    Box3dTree segtree (pmin, pmax);
    
    for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
      {
	if (mesh[si].seginfo)
	  {
	    Box<3> hbox;
	    hbox.Set (mesh[mesh[si][0]]);
	    hbox.Add (mesh[mesh[si][1]]);
	    segtree.Insert (hbox.PMin(), hbox.PMax(), si);
	  }
      }

    Array<int> loc;
    if (!ec.point_on_edge_problem)
      for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
	{
	  if (!mesh[si].seginfo) continue;

	  Box<3> hbox;
	  hbox.Set (mesh[mesh[si][0]]);
	  hbox.Add (mesh[mesh[si][1]]);
	  hbox.Increase (1e-6);
	  segtree.GetIntersecting (hbox.PMin(), hbox.PMax(), loc);
	  	  
	  // for (SegmentIndex sj = 0; sj < si; sj++)
	  for (int j = 0; j < loc.size(); j++)
	    {
	      SegmentIndex sj = loc[j];
	      if (sj >= si) continue;
	      if (!mesh[si].seginfo || !mesh[sj].seginfo) continue;
	      if (mesh[mesh[si][0]].GetLayer() != mesh[mesh[sj][1]].GetLayer()) continue;
	      
	      Point<3> pi1 = mesh[mesh[si][0]];
	      Point<3> pi2 = mesh[mesh[si][1]];
	      Point<3> pj1 = mesh[mesh[sj][0]];
	      Point<3> pj2 = mesh[mesh[sj][1]];
	      Vec<3> vi = pi2 - pi1;
	      Vec<3> vj = pj2 - pj1;
	      
	      if (sqr (vi * vj) > (1.-1e-6) * Abs2 (vi) * Abs2 (vj)) continue;
	      
	      // pi1 + vi t = pj1 + vj s
	      Mat<3,2> mat;
	      Vec<3> rhs;
	      Vec<2> sol;
	      
	      for (int jj = 0; jj < 3; jj++)
		{ 
		  mat(jj,0) = vi(jj); 
		  mat(jj,1) = -vj(jj); 
		  rhs(jj) = pj1(jj)-pi1(jj); 
		}
	      
	      mat.Solve (rhs, sol);

	      //std::cerr << "mat " << mat <<std::endl << "rhs " << rhs <<std::endl << "sol " << sol <<std::endl;
	      
	      if (sol(0) > 1e-6 && sol(0) < 1-1e-6 &&
		  sol(1) > 1e-6 && sol(1) < 1-1e-6 &&
		  Abs (rhs - mat*sol) < 1e-6)
		{
		  Point<3> ip = pi1 + sol(0) * vi;
		  
		  //std::cerr << "ip " << ip <<std::endl;

		  Point<3> pip = ip;
		  ProjectToEdge (geom.GetSurface (mesh[si].surfnr1),
				 geom.GetSurface (mesh[si].surfnr2), pip);
		  
		  //std::cerr << "Dist (ip, pip_si) " << Dist (ip, pip) <<std::endl;
		  if (Dist (ip, pip) > 1e-6*geom.MaxSize()) continue;
		  pip = ip;
		  ProjectToEdge (geom.GetSurface (mesh[sj].surfnr1),
				 geom.GetSurface (mesh[sj].surfnr2), pip);

		  //std::cerr << "Dist (ip, pip_sj) " << Dist (ip, pip) <<std::endl;
		  if (Dist (ip, pip) > 1e-6*geom.MaxSize()) continue;
		  
		  
		  
		  std::cout << "Intersection at " << ip <<std::endl;
		  
		  geom.AddUserPoint (ip);
		  spoints.push_back (MeshPoint (ip, mesh[mesh[si][0]].GetLayer()));
		  mesh.AddPoint (ip);
		  
		  std::cerr << "found intersection at " << ip <<std::endl;
		  std::cerr << "sol = " << sol <<std::endl;
		  std::cerr << "res = " << (rhs - mat*sol) <<std::endl;
		  std::cerr << "segs = " << pi1 << " - " << pi2 <<std::endl;
		  std::cerr << "and = " << pj1 << " - " << pj2 <<std::endl <<std::endl;
		}
	    }
	}  
  }






  static void MeshSurface (CSGeometry & geom, Mesh & mesh)
  {
    const char * savetask = multithread.task;
    multithread.task = "Surface meshing";
  
    Array<Segment> segments;
    int noldp = mesh.GetNP();

    // find master faces from identified
    Array<int> masterface(mesh.GetNFD());
    for (int i = 1; i <= mesh.GetNFD(); i++)
      masterface.Elem(i) = i;
  
    Array<INDEX_2> fpairs;
    bool changed;
    do
      {
	changed = 0;
	for (int i = 0; i < geom.identifications.size(); i++)
	  {
	    geom.identifications[i]->GetIdentifiedFaces (fpairs);

	    for (int j = 0; j < fpairs.size(); j++)
	      {
		if (masterface.Get(fpairs[j].I1()) <
		    masterface.Get(fpairs[j].I2()))
		  {
		    changed = 1;
		    masterface.Elem(fpairs[j].I2()) =
		      masterface.Elem(fpairs[j].I1());
		  }
		if (masterface.Get(fpairs[j].I2()) <
		    masterface.Get(fpairs[j].I1()))
		  {
		    changed = 1;
		    masterface.Elem(fpairs[j].I1()) =
		      masterface.Elem(fpairs[j].I2());
		  }
	      }
	  }
      }
    while (changed);


    int bccnt=0;
    for (int k = 0; k < geom.GetNSurf(); k++)
      bccnt = max2 (bccnt, geom.GetSurface(k)->GetBCProperty());

    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	bool increased = false;

	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	const Surface * surf = geom.GetSurface(fd.SurfNr());

	if (fd.TLOSurface() && 
	    geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetBCProp() > 0)
	  fd.SetBCProperty (geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetBCProp());
	else if (surf -> GetBCProperty() != -1)
	  fd.SetBCProperty (surf->GetBCProperty());
	else
	  {
	    bccnt++;
	    fd.SetBCProperty (bccnt);
	    increased = true;
	  }      

	for (int l = 0; l < geom.bcmodifications.size(); l++)
	  {
	    if (geom.GetSurfaceClassRepresentant (fd.SurfNr()) == 
		geom.GetSurfaceClassRepresentant (geom.bcmodifications[l].si) &&
		(fd.DomainIn() == geom.bcmodifications[l].tlonr+1 ||
		 fd.DomainOut() == geom.bcmodifications[l].tlonr+1))
	      {
		if(geom.bcmodifications[l].bcname == NULL)
		  fd.SetBCProperty (geom.bcmodifications[l].bcnr);
		else
		  {
		    if(!increased)
		      {
			bccnt++;
			fd.SetBCProperty (bccnt);
			increased = true;
		      }
		  }
	      }
	  }
      }

    mesh.SetNBCNames( bccnt );

    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	const Surface * surf = geom.GetSurface(fd.SurfNr());
	if (fd.TLOSurface() )
	  {
	    int bcp = fd.BCProperty();
	    string nextbcname = geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetBCName();
	    if ( nextbcname != "default" )
	      mesh.SetBCName ( bcp - 1 , nextbcname );
	  }
	else // if (surf -> GetBCProperty() != -1)
	  {
	    int bcp = fd.BCProperty();
	    string nextbcname = surf->GetBCName();
	    if ( nextbcname != "default" )
	      mesh.SetBCName ( bcp - 1, nextbcname );
	  }
      }
    
    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	fd.SetBCName ( mesh.GetBCNamePtr ( fd.BCProperty() - 1 ) );
      }
    

    //!!
    
    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	//const Surface * surf = geom.GetSurface(fd.SurfNr());

	for (int l = 0; l < geom.bcmodifications.size(); l++)
	  {
	    if (geom.GetSurfaceClassRepresentant (fd.SurfNr()) == 
		geom.GetSurfaceClassRepresentant (geom.bcmodifications[l].si) &&
		(fd.DomainIn() == geom.bcmodifications[l].tlonr+1 ||
		 fd.DomainOut() == geom.bcmodifications[l].tlonr+1) &&
		geom.bcmodifications[l].bcname != NULL
		)
	      {
		int bcp = fd.BCProperty();
		mesh.SetBCName ( bcp - 1, *(geom.bcmodifications[l].bcname) );
		fd.SetBCName ( mesh.GetBCNamePtr ( bcp - 1) );
	      }
	  }
      }

    for(int k = 0; k<geom.bcmodifications.size(); k++)
      {
	delete geom.bcmodifications[k].bcname;
	geom.bcmodifications[k].bcname = NULL;
      }

    //!!


    for (int j = 0; j < geom.singfaces.size(); j++)
      {
	Array<int> surfs;
	geom.GetIndependentSurfaceIndices (geom.singfaces[j]->GetSolid(),
					   geom.BoundingBox(), surfs);
	for (int k = 1; k <= mesh.GetNFD(); k++)
	  {
	    FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
	    for (int l = 0; l < surfs.size(); l++)
	      if (surfs[l] == fd.SurfNr())
		{
		  if (geom.singfaces[j]->GetDomainNr() == fd.DomainIn())
		    fd.SetDomainInSingular (1);
		  if (geom.singfaces[j]->GetDomainNr() == fd.DomainOut())
		    fd.SetDomainOutSingular (1);
		}
	  }
      }
    

    // assemble edge hash-table
    mesh.CalcSurfacesOfNode();

    for (int k = 1; k <= mesh.GetNFD(); k++)
      {
	multithread.percent = 100.0 * k / (mesh.GetNFD()+1e-10);

	if (masterface.Get(k) != k)
	  continue;

	FaceDescriptor & fd = mesh.GetFaceDescriptor(k);

	std::cerr << "Surface " << k <<std::endl;
	std::cerr << "Face Descriptor: " << fd <<std::endl;
	PrintMessage (1, "Surface ", k, " / ", mesh.GetNFD());

	int oldnf = mesh.GetNSE();
      
	const Surface * surf =
	  geom.GetSurface((mesh.GetFaceDescriptor(k).SurfNr()));


	Meshing2Surfaces meshing(*surf, mparam, geom.BoundingBox());

        double eps = 1e-8 * geom.MaxSize();
	for (PointIndex pi = PointIndex::BASE; pi < noldp+PointIndex::BASE; pi++)
	  { 
	    // if(surf->PointOnSurface(mesh[pi]))
	    meshing.AddPoint (mesh[pi], pi, NULL,
			      (surf->PointOnSurface(mesh[pi], eps) != 0));
	  }

	segments.resize (0);

	for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
	  if (mesh[si].si == k)
	    {
	      segments.push_back (mesh[si]);
	      std::cerr << "appending segment " << mesh[si] <<std::endl;
	      //<< " from " << mesh[mesh[si][0]]
	      //	 << " to " <<mesh[mesh[si][1]]<<std::endl;
	    }

	std::cerr << "num-segments " << segments.size() <<std::endl;

	for (int i = 1; i <= geom.identifications.size(); i++)
	  {
	    geom.identifications.Get(i)->
	      BuildSurfaceElements(segments, mesh, surf);
	  }

	for (int si = 0; si < segments.size(); si++)
	  {
	    PointGeomInfo gi;
	    gi.trignum = k;
	    meshing.AddBoundaryElement (segments[si][0] + 1 - PointIndex::BASE, 
					segments[si][1] + 1 - PointIndex::BASE, 
					gi, gi);
	  }

	double maxh = mparam.maxh;
	if (fd.DomainIn() != 0)
	  {
	    const Solid * s1 = 
	      geom.GetTopLevelObject(fd.DomainIn()-1) -> GetSolid();
	    if (s1->GetMaxH() < maxh)
	      maxh = s1->GetMaxH();
	    maxh = min2(maxh, geom.GetTopLevelObject(fd.DomainIn()-1)->GetMaxH());
	  }
	if (fd.DomainOut() != 0)
	  {
	    const Solid * s1 = 
	      geom.GetTopLevelObject(fd.DomainOut()-1) -> GetSolid();
	    if (s1->GetMaxH() < maxh)
	      maxh = s1->GetMaxH();
	    maxh = min2(maxh, geom.GetTopLevelObject(fd.DomainOut()-1)->GetMaxH());
	  }
	if (fd.TLOSurface() != 0)
	  {
	    double hi = geom.GetTopLevelObject(fd.TLOSurface()-1) -> GetMaxH();
	    if (hi < maxh) maxh = hi;
	  }

	std::cerr << "domin = " << fd.DomainIn() << ", domout = " << fd.DomainOut()
		   << ", tlo-surf = " << fd.TLOSurface()
		   << " mpram.maxh = " << mparam.maxh << ", maxh = " << maxh <<std::endl;

	mparam.checkoverlap = 0;

	MESHING2_RESULT res =
	  meshing.GenerateMesh (mesh, mparam, maxh, k);

	if (res != MESHING2_OK)
	  {
	    PrintError ("Problem in Surface mesh generation");
	    throw NgException ("Problem in Surface mesh generation");
	  }

	if (multithread.terminate) return;
      
	for (SurfaceElementIndex sei = oldnf; sei < mesh.GetNSE(); sei++)
	  mesh[sei].SetIndex (k);


	//      mesh.CalcSurfacesOfNode();

	if (segments.size())   
	  { 
	    // surface was meshed, not copied

	    static int timer = NgProfiler::CreateTimer ("total surface mesh optimization");
	    NgProfiler::RegionTimer reg (timer);


	    PrintMessage (2, "Optimize Surface");
	    for (int i = 1; i <= mparam.optsteps2d; i++)
	      {
		if (multithread.terminate) return;
		
		{
		  MeshOptimize2dSurfaces meshopt(geom);
		  meshopt.SetFaceIndex (k);
		  meshopt.SetImproveEdges (0);
		  meshopt.SetMetricWeight (mparam.elsizeweight);
		  meshopt.SetWriteStatus (0);
		  
		  meshopt.Edgestd::swapping (mesh, (i > mparam.optsteps2d/2));
		}
		
		if (multithread.terminate) return;
		{
		  //		mesh.CalcSurfacesOfNode();
		
		  MeshOptimize2dSurfaces meshopt(geom);
		  meshopt.SetFaceIndex (k);
		  meshopt.SetImproveEdges (0);
		  meshopt.SetMetricWeight (mparam.elsizeweight);
		  meshopt.SetWriteStatus (0);

		  meshopt.ImproveMesh (mesh, mparam);
		}
		
		{
		  MeshOptimize2dSurfaces meshopt(geom);
		  meshopt.SetFaceIndex (k);
		  meshopt.SetImproveEdges (0);
		  meshopt.SetMetricWeight (mparam.elsizeweight);
		  meshopt.SetWriteStatus (0);

		  meshopt.CombineImprove (mesh);
		  //		mesh.CalcSurfacesOfNode();
		}
		
		if (multithread.terminate) return;
		{
		  MeshOptimize2dSurfaces meshopt(geom);
		  meshopt.SetFaceIndex (k);
		  meshopt.SetImproveEdges (0);
		  meshopt.SetMetricWeight (mparam.elsizeweight);
		  meshopt.SetWriteStatus (0);

		  meshopt.ImproveMesh (mesh, mparam);
		}
	      }
	  }


	PrintMessage (3, (mesh.GetNSE() - oldnf), " elements, ", mesh.GetNP(), " points");

	extern void Render();
	Render();
      }
    
    mesh.Compress();

    do
      {
	changed = 0;
	for (int k = 1; k <= mesh.GetNFD(); k++)
	  {
	    multithread.percent = 100.0 * k / (mesh.GetNFD()+1e-10);
	  
	    if (masterface.Get(k) == k)
	      continue;

	    FaceDescriptor & fd = mesh.GetFaceDescriptor(k);

	    std::cerr << "Surface " << k <<std::endl;
	    std::cerr << "Face Descriptor: " << fd <<std::endl;
	    PrintMessage (2, "Surface ", k);

	    int oldnf = mesh.GetNSE();
      
	    const Surface * surf =
	      geom.GetSurface((mesh.GetFaceDescriptor(k).SurfNr()));

	    /*
	      if (surf -> GetBCProperty() != -1)
	      fd.SetBCProperty (surf->GetBCProperty());
	      else
	      {
	      bccnt++;
	      fd.SetBCProperty (bccnt);
	      }
	    */
  
	    segments.resize (0);
	    for (int i = 1; i <= mesh.GetNSeg(); i++)
	      {
		Segment * seg = &mesh.LineSegment(i);
		if (seg->si == k)
		  segments.push_back (*seg);
	      }

	    for (int i = 1; i <= geom.identifications.size(); i++)
	      {
		geom.identifications.Elem(i)->GetIdentifiedFaces (fpairs);
		int found = 0;
		for (int j = 1; j <= fpairs.size(); j++)
		  if (fpairs.Get(j).I1() == k || fpairs.Get(j).I2() == k)
		    found = 1;

		if (!found)
		  continue;

		geom.identifications.Get(i)->
		  BuildSurfaceElements(segments, mesh, surf);
		if (!segments.size())
		  break;
	      }

	  
	    if (multithread.terminate) return;

	    for (SurfaceElementIndex  sei = oldnf; sei < mesh.GetNSE(); sei++)
	      mesh[sei].SetIndex (k);


	    if (!segments.size())
	      {
		masterface.Elem(k) = k;
		changed = 1; 
	      }

	    PrintMessage (3, (mesh.GetNSE() - oldnf), " elements, ", mesh.GetNP(), " points");
	  }
      
	extern void Render();
	Render();
      }
    while (changed);

    
    mesh.SplitSeparatedFaces();
    mesh.CalcSurfacesOfNode();

    multithread.task = savetask;
  }



  int CSGGenerateMesh (CSGeometry & geom, 
		       Mesh *& mesh, MeshingParameters & mparam,
		       int perfstepsstart, int perfstepsend)
  {
    if (mesh && mesh->GetNSE() &&
	!geom.GetNSolids())
      {
	if (perfstepsstart < MESHCONST_MESHVOLUME)
	  perfstepsstart = MESHCONST_MESHVOLUME;
      }

    if (perfstepsstart <= MESHCONST_ANALYSE)
      {
        if (mesh)
          mesh -> DeleteMesh();
        else
          mesh = new Mesh();

	mesh->SetGlobalH (mparam.maxh);
	mesh->SetMinimalH (mparam.minh);

	Array<double> maxhdom(geom.GetNTopLevelObjects());
	for (int i = 0; i < maxhdom.size(); i++)
	  maxhdom[i] = geom.GetTopLevelObject(i)->GetMaxH();

	mesh->SetMaxHDomain (maxhdom);

	if (mparam.uselocalh)
	  {
	    double maxsize = geom.MaxSize(); 
	    mesh->SetLocalH (Point<3>(-maxsize, -maxsize, -maxsize),
			     Point<3>(maxsize, maxsize, maxsize),
			     mparam.grading);

	    mesh -> LoadLocalMeshSize (mparam.meshsizefilename);
	  }

	spoints.resize(0);
	FindPoints (geom, *mesh);
      
	PrintMessage (5, "find points done");
      }


    if (multithread.terminate || perfstepsend <= MESHCONST_ANALYSE) 
      return TCL_OK;


    if (perfstepsstart <= MESHCONST_MESHEDGES)
      {
	FindEdges (geom, *mesh, true);
	if (multithread.terminate) return TCL_OK;
	if (mparam.uselocalh)
	  {
	    mesh->CalcLocalH(mparam.grading);
	    mesh->DeleteMesh();
	    
	    FindPoints (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;
	    FindEdges (geom, *mesh, true);
	    if (multithread.terminate) return TCL_OK;
	    
	    mesh->DeleteMesh();
	  
	    FindPoints (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;
	    FindEdges (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;
	  }
      }
  
    if (multithread.terminate || perfstepsend <= MESHCONST_MESHEDGES)
      return TCL_OK;


    if (perfstepsstart <= MESHCONST_MESHSURFACE)
      {
	MeshSurface (geom, *mesh);  
	if (multithread.terminate) return TCL_OK;
	if (mparam.uselocalh && 0)
	  {
	    mesh->CalcLocalH(mparam.grading);      
	    mesh->DeleteMesh();

	    FindPoints (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;
	    FindEdges (geom, *mesh);
	    if (multithread.terminate) return TCL_OK;

	    MeshSurface (geom, *mesh);  
	    if (multithread.terminate) return TCL_OK;
	  }

	MeshQuality2d (*mesh);
	mesh->CalcSurfacesOfNode();
      }
  
    if (multithread.terminate || perfstepsend <= MESHCONST_OPTSURFACE)
      return TCL_OK;


    if (perfstepsstart <= MESHCONST_MESHVOLUME)
      {
	multithread.task = "Volume meshing";

	MESHING3_RESULT res =
	  MeshVolume (mparam, *mesh);

	if (res != MESHING3_OK) return TCL_ERROR;
      
	if (multithread.terminate) return TCL_OK;
      
	RemoveIllegalElements (*mesh);
	if (multithread.terminate) return TCL_OK;

	MeshQuality3d (*mesh);
      
	for (int i = 0; i < geom.GetNTopLevelObjects(); i++)
	  mesh->SetMaterial (i+1, geom.GetTopLevelObject(i)->GetMaterial().c_str());

      }

    if (multithread.terminate || perfstepsend <= MESHCONST_MESHVOLUME)
      return TCL_OK;


    if (perfstepsstart <= MESHCONST_OPTVOLUME)
      {
	multithread.task = "Volume optimization";
      
	OptimizeVolume (mparam, *mesh);
	if (multithread.terminate) return TCL_OK;
      }

    return TCL_OK;
  }
}
