#include <meshgen.hpp>
#include "stlgeom.hpp"

#include "../general/ngexception.hpp"
#include "../meshing/meshing3.hpp"
#include "../meshing/global.hpp"

namespace netgen
{

static void STLFindEdges (STLGeometry & geom,
			  class Mesh & mesh)
{
  double h = mparam.maxh;

  // mark edge points:
  //int ngp = geom.GetNP();

  geom.RestrictLocalH(mesh, h);
  
  PushStatusF("Mesh Lines");

  Array<STLLine*> meshlines;
  Array<Point3d> meshpoints;

  PrintMessage(3,"Mesh Lines");

  /*
  std::cout << geom.GetNLines() << " lines" <<std::endl;
  double totnp = 0;
  for (int i = 1; i <= geom.GetNLines(); i++)
    totnp += geom.GetLine(i)->NP();
  std::cout << "avg np per line " << totnp/geom.GetNLines() <<std::endl;
  */

  for (int i = 1; i <= geom.GetNLines(); i++)
    {
      meshlines.push_back(geom.GetLine(i)->Mesh(geom.GetPoints(), meshpoints, h, mesh)); 
      SetThreadPercent(100.0 * (double)i/(double)geom.GetNLines());
    }

  geom.meshpoints.resize(0); //testing
  geom.meshlines.resize(0);  //testing
  for (int i = 1; i <= meshpoints.size(); i++)
    {
      geom.meshpoints.push_back(meshpoints.Get(i)); //testing
      mesh.AddPoint(meshpoints.Get(i));
    }
  //(++++++++++++++testing
  for (int i = 1; i <= geom.GetNLines(); i++)
    {
      geom.meshlines.push_back(meshlines.Get(i));
    }
  //++++++++++++++testing)

  PrintMessage(7,"feed with edges");

  for (int i = 1; i <= meshlines.size(); i++)
    {
      STLLine* line = meshlines.Get(i);
      std::cerr << "store line " << i <<std::endl;
      for (int j = 1; j <= line->GetNS(); j++)
	{
	  int p1, p2;
	  
	  line->GetSeg(j, p1, p2);
	  int trig1, trig2, trig1b, trig2b;

	  if (p1 == p2) 
	    std::cout << "Add Segment, p1 == p2 == " << p1 <<std::endl;

	  // Test auf geschlossener Rand mit 2 Segmenten 
	      
	  if ((j == 2) && (line->GetNS() == 2))
	    {
	      int oldp1, oldp2;
	      line->GetSeg (1, oldp1, oldp2);
	      if (oldp1 == p2 && oldp2 == p1)
		{
		  PrintMessage(7,"MESSAGE: don't use second segment");
		  continue;
		}
	    }


	  //mesh point number
	  //p1 = geom2meshnum.Get(p1); // for unmeshed lines!!!
	  //p2 = geom2meshnum.Get(p2); // for unmeshed lines!!!
	  
	  //left and right trigs
	  trig1 = line->GetLeftTrig(j);
	  trig2 = line->GetRightTrig(j);
	  trig1b = line->GetLeftTrig(j+1);
	  trig2b = line->GetRightTrig(j+1);
	  
	  std::cerr << "j = " << j << ", p1 = " << p1 << ", p2 = " << p2 <<std::endl;
	  std::cerr << "segm-trigs: "
		   << "trig1 = " << trig1
		   << ", trig1b = " << trig1b
		   << ", trig2 = " << trig2
		   << ", trig2b = " << trig2b <<std::endl;

	  if (trig1 <= 0 || trig2 <= 0 || trig1b <= 0 || trig2b <= 0)
	    {
	      std::cout << "negative trigs, "
		   << ", trig1 = " << trig1
		   << ", trig1b = " << trig1b
		   << ", trig2 = " << trig2
		   << ", trig2b = " << trig2b <<std::endl;
	    }
	  /*
	  std::cerr << "   trigs p1: " << trig1 << " - " << trig2 <<std::endl;
	  std::cerr << "   trigs p2: " << trig1b << " - " << trig2b <<std::endl;
	  std::cerr << "   charts p1: " << geom.GetChartNr(trig1) << " - " << geom.GetChartNr(trig2) <<std::endl;
	  std::cerr << "   charts p2: " << geom.GetChartNr(trig1b) << " - " << geom.GetChartNr(trig2b) <<std::endl;
	  */
	  Point3d hp, hp2;
	  Segment seg;
	  seg[0] = p1;
	  seg[1] = p2;
	  seg.si = geom.GetTriangle(trig1).GetFaceNum();
	  seg.edgenr = i;

	  seg.epgeominfo[0].edgenr = i;
	  seg.epgeominfo[0].dist = line->GetDist(j);
	  seg.epgeominfo[1].edgenr = i;
	  seg.epgeominfo[1].dist = line->GetDist(j+1);
	  /*
	  std::cerr << "seg = " 
		     << "edgenr " << seg.epgeominfo[0].edgenr
		     << " dist " << seg.epgeominfo[0].dist
		     << " edgenr " << seg.epgeominfo[1].edgenr
		     << " dist " << seg.epgeominfo[1].dist <<std::endl;
	  */
	  
	  seg.geominfo[0].trignum = trig1;
	  seg.geominfo[1].trignum = trig1b;

	  /*
	  geom.SelectChartOfTriangle (trig1);
	  hp = hp2 = mesh.Point (seg[0]);
	  seg.geominfo[0].trignum = geom.Project (hp);

	  std::cerr << "hp = " << hp2 << ", hp proj = " << hp << ", trignum = " << seg.geominfo[0].trignum <<std::endl;
	  if (Dist (hp, hp2) > 1e-5 || seg.geominfo[0].trignum == 0) 
	    {
	      std::cerr << "PROBLEM" <<std::endl;
	    }

	  geom.SelectChartOfTriangle (trig1b);
	  hp = hp2 = mesh.Point (seg[1]);
	  seg.geominfo[1].trignum = geom.Project (hp);

	  std::cerr << "hp = " << hp2 << ", hp proj = " << hp << ", trignum = " << seg.geominfo[1].trignum <<std::endl;
	  if (Dist (hp, hp2) > 1e-5 || seg.geominfo[1].trignum == 0) 
	    {
	      std::cerr << "PROBLEM" <<std::endl;
	    }
	  */


	  if (Dist (mesh.Point(seg[0]), mesh.Point(seg[1])) < 1e-10)
	    {
	      std::cerr << "ERROR: Line segment of length 0" <<std::endl;
	      std::cerr << "pi1, 2 = " << seg[0] << ", " << seg[1] <<std::endl;
	      std::cerr << "p1, 2 = " << mesh.Point(seg[0])
			 << ", " << mesh.Point(seg[1]) <<std::endl;
	      throw NgException ("Line segment of length 0");
	    }
	  
	  mesh.AddSegment (seg);


	  Segment seg2;
	  seg2[0] = p2;
	  seg2[1] = p1;
	  seg2.si = geom.GetTriangle(trig2).GetFaceNum();
	  seg2.edgenr = i;

	  seg2.epgeominfo[0].edgenr = i;
	  seg2.epgeominfo[0].dist = line->GetDist(j+1);
	  seg2.epgeominfo[1].edgenr = i;
	  seg2.epgeominfo[1].dist = line->GetDist(j);
	  /*
	  std::cerr << "seg = " 
		     << "edgenr " << seg2.epgeominfo[0].edgenr
		     << " dist " << seg2.epgeominfo[0].dist
		     << " edgenr " << seg2.epgeominfo[1].edgenr
		     << " dist " << seg2.epgeominfo[1].dist <<std::endl;
	  */
	  
	  seg2.geominfo[0].trignum = trig2b;
	  seg2.geominfo[1].trignum = trig2;
	  
	  /*
	  geom.SelectChartOfTriangle (trig2);
	  hp = hp2 = mesh.Point (seg[0]);
	  seg2.geominfo[0].trignum = geom.Project (hp);

	  std::cerr << "hp = " << hp2 << ", hp proj = " << hp << ", trignum = " << seg.geominfo[0].trignum <<std::endl;
	  if (Dist (hp, hp2) > 1e-5 || seg2.geominfo[0].trignum == 0) 
	    {
	      std::cerr << "Get GeomInfo PROBLEM" <<std::endl;
	    }


	  geom.SelectChartOfTriangle (trig2b);
	  hp = hp2 = mesh.Point (seg[1]);
	  seg2.geominfo[1].trignum = geom.Project (hp);
	  std::cerr << "hp = " << hp2 << ", hp proj = " << hp << ", trignum = " << seg.geominfo[1].trignum <<std::endl;
	  if (Dist (hp, hp2) > 1e-5 || seg2.geominfo[1].trignum == 0) 
	    {
	      std::cerr << "Get GeomInfo PROBLEM" <<std::endl;
	    }
	  */	  

	  mesh.AddSegment (seg2);
	}
    }

  PopStatus();
}




void STLSurfaceMeshing1 (STLGeometry & geom, class Mesh & mesh,
			 int retrynr);


int STLSurfaceMeshing (STLGeometry & geom, class Mesh & mesh)
{
  PrintFnStart("Do Surface Meshing");

  geom.PrepareSurfaceMeshing();

  if (mesh.GetNSeg() == 0)
    STLFindEdges (geom, mesh);

  int nopen;
  int outercnt = 20;

  for (int i = 1; i <= mesh.GetNSeg(); i++)
    {
      const Segment & seg = mesh.LineSegment (i);
      if (seg.geominfo[0].trignum <= 0 || seg.geominfo[1].trignum <= 0)
	std::cerr << "Problem with segment " << i << ": " << seg <<std::endl;
    }


  do
    {
      outercnt--;
      if (outercnt <= 0)
	return MESHING3_OUTERSTEPSEXCEEDED;
      
      if (multithread.terminate) return MESHING3_TERMINATE;

      mesh.FindOpenSegments();
      nopen = mesh.GetNOpenSegments();

      if (nopen)
	{
	  int trialcnt = 0;
	  while (nopen && trialcnt <= 5)
	    {
	      if (multithread.terminate) { return MESHING3_TERMINATE; }

	      trialcnt++;
	      STLSurfaceMeshing1 (geom, mesh, trialcnt);

	      mesh.FindOpenSegments();
	      nopen = mesh.GetNOpenSegments();

	      if (nopen)
		{
		  geom.ClearMarkedSegs();
		  for (int i = 1; i <= nopen; i++)
		    {
		      const Segment & seg = mesh.GetOpenSegment (i);
		      geom.AddMarkedSeg(mesh.Point(seg[0]),mesh.Point(seg[1]));
		    }

		  geom.InitMarkedTrigs();
		  for (int i = 1; i <= nopen; i++)
		    {
		      const Segment & seg = mesh.GetOpenSegment (i);
		      geom.SetMarkedTrig(seg.geominfo[0].trignum,1);
		      geom.SetMarkedTrig(seg.geominfo[1].trignum,1);
		    }

		  MeshOptimizeSTLSurface optmesh(geom);
		  optmesh.SetFaceIndex (0);
		  optmesh.SetImproveEdges (0);
		  optmesh.SetMetricWeight (0);
		  
		  mesh.CalcSurfacesOfNode();
		  optmesh.EdgeSwapping (mesh, 0);
		  mesh.CalcSurfacesOfNode();
		  optmesh.ImproveMesh (mesh, mparam);
		}

	      mesh.Compress();
	      mesh.FindOpenSegments();
	      nopen = mesh.GetNOpenSegments();

	      if (trialcnt <= 5 && nopen)
		{
		  mesh.RemoveOneLayerSurfaceElements();

		  if (trialcnt >= 4)
		    {
		      mesh.FindOpenSegments();
		      mesh.RemoveOneLayerSurfaceElements();

		      mesh.FindOpenSegments ();		  
		      nopen = mesh.GetNOpenSegments();
		    }
		}
	    }


	  if (multithread.terminate)
	    return MESHING3_TERMINATE;

	  if (nopen)
	    {
	      
	      PrintMessage(3,"Meshing failed, trying to refine");

	      mesh.FindOpenSegments ();
	      nopen = mesh.GetNOpenSegments();
			  
	      mesh.FindOpenSegments ();
	      mesh.RemoveOneLayerSurfaceElements();
	      mesh.FindOpenSegments ();
	      mesh.RemoveOneLayerSurfaceElements();

	      // Open edge-segments will be refined !
	      INDEX_2_HASHTABLE<int> openseght (nopen+1);
	      for (int i = 1; i <= mesh.GetNOpenSegments(); i++)
		{
		  const Segment & seg = mesh.GetOpenSegment (i);
		  INDEX_2 i2(seg[0], seg[1]);
		  i2.Sort();
		  openseght.Set (i2, 1);
		}

	      
	      mesh.FindOpenSegments ();
	      mesh.RemoveOneLayerSurfaceElements();
	      mesh.FindOpenSegments ();
	      mesh.RemoveOneLayerSurfaceElements();
	      

	      INDEX_2_HASHTABLE<int> newpht(100);

	      int nsegold = mesh.GetNSeg();
	      for (int i = 1; i <= nsegold; i++)
		{
		  Segment seg = mesh.LineSegment(i);
		  INDEX_2 i2(seg[0], seg[1]);
		  i2.Sort();
		  if (openseght.Used (i2))
		    {
		      // segment will be split
		      PrintMessage(7,"Split segment ", int(seg[0]), "-", int(seg[1]));
	      
		      Segment nseg1, nseg2;
		      EdgePointGeomInfo newgi;
		      
		      const EdgePointGeomInfo & gi1 = seg.epgeominfo[0];
		      const EdgePointGeomInfo & gi2 = seg.epgeominfo[1];
		      
		      newgi.dist = 0.5 * (gi1.dist + gi2.dist);
		      newgi.edgenr = gi1.edgenr;

		      int hi;
		      
		      Point3d newp;
		      int newpi;
		      
		      if (!newpht.Used (i2))
			{
			  newp = geom.GetLine (gi1.edgenr)->
			    GetPointInDist (geom.GetPoints(), newgi.dist, hi);
			  newpi = mesh.AddPoint (newp);
			  newpht.Set (i2, newpi);
			}
		      else
			{
			  newpi = newpht.Get (i2);
			  newp = mesh.Point (newpi);
			}

		      nseg1 = seg;
		      nseg2 = seg;
		      nseg1[1] = newpi;
		      nseg1.epgeominfo[1] = newgi;
		      
		      nseg2[0] = newpi;
		      nseg2.epgeominfo[0] = newgi;
		      
		      mesh.LineSegment(i) = nseg1;
		      mesh.AddSegment (nseg2);
		      
		      mesh.RestrictLocalH (Center (mesh.Point(nseg1[0]),
						   mesh.Point(nseg1[1])),
					   Dist (mesh.Point(nseg1[0]),
						 mesh.Point(nseg1[1])));
		      mesh.RestrictLocalH (Center (mesh.Point(nseg2[0]),
						   mesh.Point(nseg2[1])),
					   Dist (mesh.Point(nseg2[0]),
						 mesh.Point(nseg2[1])));
		    }
		}

	    }

	  nopen = -1;
	}
    
      else

	{
	  PrintMessage(5,"mesh is closed, verifying ...");

	  // no open elements, check wrong elemetns (intersecting..)



	  PrintMessage(5,"check overlapping");
	  // 	  mesh.FindOpenElements(); // would leed to locked points
	  if(mesh.CheckOverlappingBoundary())
	    return MESHING3_BADSURFACEMESH;

	  geom.InitMarkedTrigs();

	  for (int i = 1; i <= mesh.GetNSE(); i++)
	    if (mesh.SurfaceElement(i).BadElement())
	      {
		int trig = mesh.SurfaceElement(i).PNum(1);
		geom.SetMarkedTrig(trig,1);
		PrintMessage(7, "overlapping element, will be removed");
	      }
	  
	  

	  Array<Point3d> refpts;
	  Array<double> refh;

	  // was commented:

	  for (int i = 1; i <= mesh.GetNSE(); i++)
	    if (mesh.SurfaceElement(i).BadElement())
	      {
		for (int j = 1; j <= 3; j++)
		  {
		    refpts.push_back (mesh.Point (mesh.SurfaceElement(i).PNum(j)));
		    refh.push_back (mesh.GetH (refpts.Last()) / 2);
		  }
		mesh.DeleteSurfaceElement(i);
	      }
	  	  
	  // delete wrong oriented element
	  for (int i = 1; i <= mesh.GetNSE(); i++)
	    {
	      const Element2d & el = mesh.SurfaceElement(i);
	      if (!el.PNum(1))
		continue;

	      Vec3d n = Cross (Vec3d (mesh.Point(el.PNum(1)), 
				      mesh.Point(el.PNum(2))),
			       Vec3d (mesh.Point(el.PNum(1)), 
				      mesh.Point(el.PNum(3))));
	      Vec3d ng = geom.GetTriangle(el.GeomInfoPi(1).trignum).Normal();
	      if (n * ng < 0)
		{
		  refpts.push_back (mesh.Point (mesh.SurfaceElement(i).PNum(1)));
		  refh.push_back (mesh.GetH (refpts.Last()) / 2);
		  mesh.DeleteSurfaceElement(i);
		}
	    }
	  // end comments

	  for (int i = 1; i <= refpts.size(); i++)
	    mesh.RestrictLocalH (refpts.Get(i), refh.Get(i));

	  mesh.RemoveOneLayerSurfaceElements();

	  mesh.Compress();
	  
	  mesh.FindOpenSegments ();
	  nopen = mesh.GetNOpenSegments();

	  /*
	  if (!nopen)
	    {
	      // mesh is still ok

	      void STLSurfaceOptimization (STLGeometry & geom,
					   class Mesh & mesh,
					   MeshingParameters & mparam)
	      
	    }
	  */
	}
      
    }
  while (nopen);

  mesh.Compress();
  mesh.CalcSurfacesOfNode();

  return MESHING3_OK;
}






void STLSurfaceMeshing1 (STLGeometry & geom,
			 class Mesh & mesh,
			 int retrynr)
{
  double h = mparam.maxh;

  mesh.FindOpenSegments();
  
  Array<int> spiralps(0);
  spiralps.resize(0);
  for (int i = 1; i <= geom.GetNP(); i++)
    if (geom.GetSpiralPoint(i)) 
      spiralps.push_back(i);
  
  PrintMessage(7,"NO spiralpoints = ", spiralps.size());
  //int spfound;

  /*
  Array<int> meshsp(mesh.GetNP());
  meshsp = 0;
  for (int i = 1; i <= mesh.GetNP(); i++)
    for (int j = 1; j <= spiralps.Size(); j++)
      if (Dist2(geom.GetPoint(spiralps.Get(j)), mesh.Point(i)) < 1e-20) 
	meshsp.Elem(i) = spiralps.Get(j);
  Array<PointIndex> imeshsp;
  for (int i = 1; i <= meshsp.Size(); i++)
    if (meshsp.Elem(i)) imeshsp.Append(i);
  */
  Array<PointIndex> imeshsp;
  Array<int> ispiral_point;
  for (int i = 1; i <= mesh.GetNP(); i++)
    {
      for (int j = 1; j <= spiralps.size(); j++)
	if (Dist2(geom.GetPoint(spiralps.Get(j)), mesh.Point(i)) < 1e-20) 
	  {
	    imeshsp.push_back(i);
	    ispiral_point.push_back(spiralps.Get(j));
	    break;
	  }
    }
  mesh.SurfaceArea().ReCalc();

  // int oldnp = mesh.GetNP();

  Array<int,PointIndex::BASE> compress(mesh.GetNP());
  compress = 0;
  Array<PointIndex> icompress; 

  Array<int, 1> opensegsperface(mesh.GetNFD());
  opensegsperface = 0;
  for (int i = 1; i <= mesh.GetNOpenSegments(); i++)
    opensegsperface[mesh.GetOpenSegment(i).si]++;
  
  TABLE<int, 1> opensegments(mesh.GetNFD());
  for (int i = 1; i <= mesh.GetNOpenSegments(); i++)
    {
      const Segment & seg = mesh.GetOpenSegment (i);
      if (seg.si < 1 || seg.si > mesh.GetNFD())
	std::cerr << "segment index " << seg.si << " out of range [1, " << mesh.GetNFD() << "]" <<std::endl;
      opensegments.Add (seg.si, i);
    }
  

  for (int fnr = 1; fnr <= mesh.GetNFD(); fnr++)
    {
      if (!opensegsperface[fnr]) continue;
      if (multithread.terminate) return;

      PrintMessage(5,"Meshing surface ", fnr, "/", mesh.GetNFD());
      MeshingSTLSurface meshing (geom, mparam);
      
      // compress = 0;
      icompress.resize(0); 
      int cntused = 0;

      for (int i = 0; i < imeshsp.size(); i++)
	{
	  compress[imeshsp[i]] = ++cntused;
	  icompress.push_back(imeshsp[i]);
	}

      FlatArray<int> segs = opensegments[fnr];
      for (int hi = 0; hi < segs.size(); hi++)
	{
	  int i = segs[hi];
	  const Segment & seg = mesh.GetOpenSegment (i);
	  for (int j = 0; j < 2; j++)
	    if (compress[seg[j]] == 0)
	      {
		compress[seg[j]] = ++cntused;
		icompress.push_back(seg[j]);
	      }
	}

      for (int hi = 0; hi < icompress.size(); hi++)
	{
	  PointIndex pi = icompress[hi];
	  
	  /*
	  // int sppointnum = meshsp.Get(i);
	  int sppointnum = 0;
	  if (hi < ispiral_point.Size())
	    sppointnum = ispiral_point[hi];

	  if (sppointnum)
	    {
	  */
	  if (hi < ispiral_point.size())
	    {
	      int sppointnum = ispiral_point[hi];

	      MultiPointGeomInfo mgi;
	      
	      int ntrigs = geom.NOTrigsPerPoint(sppointnum);
	      for (int j = 0; j < ntrigs; j++)
		{
		  PointGeomInfo gi;
		  gi.trignum = geom.TrigPerPoint(sppointnum, j+1);
		  mgi.AddPointGeomInfo (gi);
		}
	      
	      // Einfuegen von ConePoint: Point bekommt alle
	      // Dreiecke (werden dann intern kopiert)
	      // Ein Segment zum ConePoint muss vorhanden sein !!!
	      
	      // meshing.AddPoint (mesh.Point(i), i, &mgi);
	      meshing.AddPoint (mesh[pi], pi, &mgi);
	    }
	  else
	    meshing.AddPoint (mesh[pi], pi);
	}

      // FlatArray<int> segs = opensegments[fnr];
      for (int hi = 0; hi < segs.size(); hi++)
	{
	  int i = segs[hi];
	  const Segment & seg = mesh.GetOpenSegment (i);
	  meshing.AddBoundaryElement (compress[seg[0]], compress[seg[1]], 
				      seg.geominfo[0], seg.geominfo[1]);
	}



      PrintMessage(3,"start meshing, trialcnt = ", retrynr);
      
      meshing.GenerateMesh (mesh, mparam, h, fnr);  
      
      for (int i = 0; i < icompress.size(); i++)
	compress[icompress[i]] = 0;
      
      
      extern void Render();
      Render();
    }     
  
  // NgProfiler::Print(stdout);
  
  mesh.CalcSurfacesOfNode();
}



void STLSurfaceOptimization (STLGeometry & geom,
			     class Mesh & mesh,
			     MeshingParameters & meshparam)
{
  PrintFnStart("optimize STL Surface");

  MeshOptimizeSTLSurface optmesh(geom);

  optmesh.SetFaceIndex (0);
  optmesh.SetImproveEdges (0);
  optmesh.SetMetricWeight (meshparam.elsizeweight);

  PrintMessage(5,"optimize string = ", meshparam.optimize2d, " elsizew = ", meshparam.elsizeweight);

  for (int i = 1; i <= meshparam.optsteps2d; i++)
    for (size_t j = 1; j <= strlen(meshparam.optimize2d); j++)
      {
	if (multithread.terminate)
	  break;

	//std::cerr << "optimize, before, step = " << meshparam.optimize2d[j-1] << mesh.Point (3679) <<std::endl;

	mesh.CalcSurfacesOfNode();
	switch (meshparam.optimize2d[j-1])
	  {
	  case 's': 
	    {
	      optmesh.EdgeSwapping (mesh, 0);
	      break;
	    }
	  case 'S': 
	    {
	      optmesh.EdgeSwapping (mesh, 1);
	      break;
	    }
	  case 'm': 
	    {
	      optmesh.ImproveMesh(mesh, mparam);
	      break;
	    }
	  case 'c': 
	    {
	      optmesh.CombineImprove (mesh);
	      break;
	    }
	  }
	//std::cerr << "optimize, after, step = " << meshparam.optimize2d[j-1] << mesh.Point (3679) <<std::endl;
      }

  geom.surfaceoptimized = 1;

  mesh.Compress();
  mesh.CalcSurfacesOfNode();


}



MeshingSTLSurface :: MeshingSTLSurface (STLGeometry & ageom,
					const MeshingParameters & mp)
  : Meshing2(mp, ageom.GetBoundingBox()), geom(ageom)
{
  ;
}

void MeshingSTLSurface :: DefineTransformation (const Point3d & p1, const Point3d & p2,
						const PointGeomInfo * geominfo,
						const PointGeomInfo * geominfo2)
{
  transformationtrig = geominfo[0].trignum;
  
  geom.DefineTangentialPlane(p1, p2, transformationtrig);
}

void MeshingSTLSurface :: TransformToPlain (const Point3d & locpoint, const MultiPointGeomInfo & gi,
					    Point2d & plainpoint, double h, int & zone)
{
  int trigs[10000];

  if (gi.GetNPGI() >= 9999) 
    {
      PrintError("In Transform to plane: increase size of trigs!!!");
    }

  for (int i = 1; i <= gi.GetNPGI(); i++)
    trigs[i-1] = gi.GetPGI(i).trignum;
  trigs[gi.GetNPGI()] = 0;

  //  int trig = gi.trignum;
  //   std::cerr << "locpoint = " << locpoint;

  Point<2> hp2d;
  geom.ToPlane (locpoint, trigs, hp2d, h, zone, 1);
  plainpoint = hp2d;

  //  geom.ToPlane (locpoint, NULL, plainpoint, h, zone, 1);
  /*
  std::cerr << " plainpoint = " << plainpoint
	     << " h = " << h 
	     <<std::endl;
  */
}

/*
int MeshingSTLSurface :: ComputeLineGeoInfo (const Point3d & p1, const Point3d & p2,
					      int & geoinfosize, void *& geoinfo)
{
  static int geomtrig[2] = { 0, 0 };

  Point3d hp;
  hp = p1;
  geomtrig[0] = geom.Project (hp);

  hp = p2;
  geomtrig[1] = geom.Project (hp);
  
  geoinfosize = sizeof (geomtrig);
  geoinfo = &geomtrig;

  if (geomtrig[0] == 0)
    {
      return 1;
    }
  return 0;
}
*/


int MeshingSTLSurface :: ComputePointGeomInfo (const Point3d & p, PointGeomInfo & gi)
{
  // compute triangle of point,
  // if non-unique: 0

  Point<3> hp = p;
  gi.trignum = geom.Project (hp);

  if (!gi.trignum)
    {
      return 1;
    }

  return 0;
}


int MeshingSTLSurface :: 
ChooseChartPointGeomInfo (const MultiPointGeomInfo & mpgi, 
			  PointGeomInfo & pgi)
{
  for (int i = 1; i <= mpgi.GetNPGI(); i++)
    if (geom.TrigIsInOC (mpgi.GetPGI(i).trignum, geom.meshchart))
      {
	pgi = mpgi.GetPGI(i);
	return 0;
      }
  /*
  for (i = 0; i < mpgi.cnt; i++)
    {
      //      std::cerr << "d" <<std::endl;
      if (geom.TrigIsInOC (mpgi.mgi[i].trignum, geom.meshchart))
	{
	  pgi = mpgi.mgi[i];
	  return 0;
	}
    }
  */
  PrintMessage(7,"INFORM: no gi on chart");
  pgi.trignum = 1;
  return 1;
}



int MeshingSTLSurface :: 
IsLineVertexOnChart (const Point3d & p1, const Point3d & p2,
		     int endpoint, const PointGeomInfo & gi)
{
  int lineendtrig = gi.trignum;
  return geom.TrigIsInOC (lineendtrig, geom.meshchart);

  // Vec3d baselinenormal = geom.meshtrignv;
  //  Vec3d linenormal = geom.GetTriangleNormal (lineendtrig);
  //  return ( (baselinenormal * linenormal) > cos (30 * (M_PI/180)) );
}

void MeshingSTLSurface :: 
GetChartBoundary (Array<Point2d > & points, 
		  Array<Point3d > & points3d,
		  Array<INDEX_2> & lines, double h) const
{
  points.resize (0);
  points3d.resize (0);
  lines.resize (0);
  geom.GetMeshChartBoundary (points, points3d, lines, h);
}




int MeshingSTLSurface :: TransformFromPlain (Point2d & plainpoint,
					     Point3d & locpoint, 
					     PointGeomInfo & gi, 
					     double h)
{
  //return 0, wenn alles OK
  Point<3> hp3d;
  int res = geom.FromPlane (plainpoint, hp3d, h);
  locpoint = hp3d;
  ComputePointGeomInfo (locpoint, gi);
  return res;
}


int MeshingSTLSurface :: 
BelongsToActiveChart (const Point3d & p, 
		      const PointGeomInfo & gi)
{
  return (geom.TrigIsInOC(gi.trignum, geom.meshchart) != 0);
}



double MeshingSTLSurface :: CalcLocalH (const Point3d & p, double gh) const
{
  return gh;
}

double MeshingSTLSurface :: Area () const
{
  return geom.Area();
}






MeshOptimizeSTLSurface :: MeshOptimizeSTLSurface (STLGeometry & ageom)
  : MeshOptimize2d(), geom(ageom)
{
  ;
}


void MeshOptimizeSTLSurface :: SelectSurfaceOfPoint (const Point<3> & p,
						     const PointGeomInfo & gi)
{
  //  std::cerr << "sel char: " << gi.trignum <<std::endl;
  
  geom.SelectChartOfTriangle (gi.trignum);
  //  geom.SelectChartOfPoint (p);
}


void MeshOptimizeSTLSurface :: ProjectPoint (INDEX surfind, Point<3> & p) const
{
  if (!geom.Project (p))
    {
      PrintMessage(7,"project failed");
      
      if (!geom.ProjectOnWholeSurface(p)) 
	{
	  PrintMessage(7, "project on whole surface failed");
	}
    }

  //  geometry.GetSurface(surfind)->Project (p);
}

void MeshOptimizeSTLSurface :: ProjectPoint2 (INDEX surfind, INDEX surfind2, Point<3> & p) const
{
  /*
  ProjectToEdge ( geometry.GetSurface(surfind), 
		  geometry.GetSurface(surfind2), p);
  */
}

int  MeshOptimizeSTLSurface :: CalcPointGeomInfo(PointGeomInfo& gi, const Point<3> & p3) const
{
  Point<3> hp = p3;
  gi.trignum = geom.Project (hp);

  if (gi.trignum) return 1;

  return 0;
  
}

void MeshOptimizeSTLSurface :: GetNormalVector(INDEX surfind, const Point<3> & p, Vec<3> & n) const
{
  n = geom.GetChartNormalVector();
}
  









RefinementSTLGeometry :: RefinementSTLGeometry (const STLGeometry & ageom)
  : Refinement(), geom(ageom)
{
  ;
}

RefinementSTLGeometry :: ~RefinementSTLGeometry ()
{
  ;
}
  
void RefinementSTLGeometry :: 
PointBetween  (const Point<3> & p1, const Point<3> & p2, double secpoint,
	       int surfi, 
	       const PointGeomInfo & gi1, 
	       const PointGeomInfo & gi2,
	       Point<3> & newp, PointGeomInfo & newgi) const
{
  newp = p1+secpoint*(p2-p1);

  /*
  std::cerr << "surf-between: p1 = " << p1 << ", p2 = " << p2
	     << ", gi = " << gi1 << " - " << gi2 <<std::endl;
  */

  if (gi1.trignum > 0)
    {
      //      ((STLGeometry&)geom).SelectChartOfTriangle (gi1.trignum);

      Point<3> np1 = newp;
      Point<3> np2 = newp;
      ((STLGeometry&)geom).SelectChartOfTriangle (gi1.trignum);
      int tn1 = geom.Project (np1);

      ((STLGeometry&)geom).SelectChartOfTriangle (gi2.trignum);
      int tn2 = geom.Project (np2);

      newgi.trignum = tn1; //urspruengliche version
      newp = np1;          //urspruengliche version

      if (!newgi.trignum) 
	{ newgi.trignum = tn2; newp = np2; }
      if (!newgi.trignum) newgi.trignum = gi1.trignum;

      /*    
      if (tn1 != 0 && tn2 != 0 && ((STLGeometry&)geom).GetAngle(tn1,tn2) < M_PI*0.05)	{
	  newgi.trignum = tn1;
	  newp = np1;
	}
      else
	{
	  newp = ((STLGeometry&)geom).PointBetween(p1, gi1.trignum, p2, gi2.trignum);
	  tn1 = ((STLGeometry&)geom).Project(newp);
	  newgi.trignum = tn1;

	  if (!tn1) 
	    {
	      newp = Center (p1, p2);
	      newgi.trignum = 0;
	      
	    }
	}
      */
    }
  else
    {
      //      std::cerr << "WARNING: PointBetween got geominfo = 0" <<std::endl;
      newp =  p1+secpoint*(p2-p1);
      newgi.trignum = 0;
    }
     
  //  std::cerr << "newp = " << newp << ", ngi = " << newgi <<std::endl;
}

void RefinementSTLGeometry ::
PointBetween (const Point<3> & p1, const Point<3> & p2, double secpoint,
	      int surfi1, int surfi2, 
	      const EdgePointGeomInfo & gi1, 
	      const EdgePointGeomInfo & gi2,
	      Point<3> & newp, EdgePointGeomInfo & newgi) const
{
  /*
  std::cerr << "edge-between: p1 = " << p1 << ", p2 = " << p2
	     << ", gi1,2 = " << gi1 << ", " << gi2 <<std::endl;
  */
  /*
  newp = Center (p1, p2);
  ((STLGeometry&)geom).SelectChartOfTriangle (gi1.trignum);
  newgi.trignum = geom.Project (newp);
  */
  int hi;
  newgi.dist = (1.0-secpoint) * gi1.dist + secpoint*gi2.dist;
  newgi.edgenr = gi1.edgenr;

  /*
  std::cerr << "p1 = " << p1 << ", p2 = " << p2 <<std::endl;
  std::cerr << "refedge: " << gi1.edgenr
	     << " d1 = " << gi1.dist << ", d2 = " << gi2.dist <<std::endl;
  */
  newp = geom.GetLine (gi1.edgenr)->GetPointInDist (geom.GetPoints(), newgi.dist, hi);

  //  std::cerr << "newp = " << newp <<std::endl;
}


void RefinementSTLGeometry :: ProjectToSurface (Point<3> & p, int surfi) const
{
  std::cout << "RefinementSTLGeometry :: ProjectToSurface not implemented!" <<std::endl;
};


void RefinementSTLGeometry :: ProjectToSurface (Point<3> & p, int surfi,
						PointGeomInfo & gi) const
{
  ((STLGeometry&)geom).SelectChartOfTriangle (gi.trignum);
  gi.trignum = geom.Project (p);
  //  if (!gi.trignum) 
  //    std::cout << "projectSTL failed" <<std::endl;
};

 
}
