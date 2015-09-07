#include <sstream>
#include <meshit.hpp>
#include "validate.hpp"
#include "bisect.hpp"
#include "improve2.hpp"
#include "improve3.hpp"
#include "global.hpp"

namespace meshit
{
  void GetPureBadness(Mesh & mesh, Array<double> & pure_badness,
		      const BitArray & isnewpoint)
  {
    //const int ne = mesh.GetNE();
    const int np = mesh.GetNP();

    pure_badness.resize(np+PointIndex::BASE+1);
    pure_badness = -1;

    Array< Point<3>* > backup(np);

    for(int i=0; i<np; i++)
      {
	backup[i] = new Point<3>(mesh.Point(i+1));

	if(isnewpoint.Test(i+PointIndex::BASE) &&
	   mesh.mlbetweennodes[i+PointIndex::BASE][0] > 0)
	  {
	    mesh.Point(i+1) = Center(mesh.Point(mesh.mlbetweennodes[i+PointIndex::BASE][0]),
				     mesh.Point(mesh.mlbetweennodes[i+PointIndex::BASE][1]));
	  }
      }
    for (ElementIndex i = 0; i < mesh.GetNE(); i++)
      {
	double bad = mesh[i].CalcJacobianBadness (mesh.Points());
	for(int j=0; j<mesh[i].GetNP(); j++)
	  if(bad > pure_badness[mesh[i][j]])
	    pure_badness[mesh[i][j]] = bad;

	// save maximum
	if(bad > pure_badness.Last())
	  pure_badness.Last() = bad; 
      }
    
    for(int i=0; i<np; i++)
      {
	mesh.Point(i+1) = *backup[i];
	delete backup[i];
      }
  }


  double Validate(const Mesh & mesh, Array<ElementIndex> & bad_elements,
		  const Array<double> & pure_badness,
		  double max_worsening, const bool uselocalworsening,
		  Array<double> * quality_loss)
  {
    PrintMessage(3,"!!!! Validating !!!!");
    //if(max_worsening > 0)
    //  std::cerr << "badness " << counter++ <<std::endl;

    bad_elements.resize(0);

    double loc_pure_badness = -1;

    if(!uselocalworsening)
      loc_pure_badness = pure_badness.Last(); // maximum is saved at last position


    double worsening = -1;
    ElementIndex ind;

    if(quality_loss != NULL)
      quality_loss->resize(mesh.GetNE());

    for (ElementIndex i = 0; i < mesh.GetNE(); i++)
      {
	if(uselocalworsening)
	  {
	    loc_pure_badness = -1;
	    for(int j=0; j<mesh[i].GetNP(); j++)
	      if(pure_badness[mesh[i][j]] > loc_pure_badness)
		loc_pure_badness = pure_badness[mesh[i][j]];
	  }


	double bad = mesh[i].CalcJacobianBadness (mesh.Points());
	if (bad > 1e10 || 
	    (max_worsening > 0 && bad > loc_pure_badness*max_worsening))
	  bad_elements.push_back(i);
	  

	if(max_worsening > 0)
	  {
	    double actw = bad/loc_pure_badness;
	    if(quality_loss != NULL)
	      (*quality_loss)[i] = actw;

	    if(actw > worsening)
	      {
		worsening = actw;
		ind = i;
	      }
	  }
      }
    return worsening;
  }


  void GetWorkingArea(BitArray & working_elements, BitArray & working_points,
		      const Mesh & mesh, const Array<ElementIndex> & bad_elements,
		      const int width)
  {
    working_elements.Clear();
    working_points.Clear();

    for(int i=0; i<bad_elements.size(); i++)
      {
	working_elements.Set(bad_elements[i]);
	const Element & el = mesh[bad_elements[i]];
	for(int j=1; j<=el.GetNP(); j++)
	  working_points.Set(el.PNum(j));
      }
    

    for(int i=0; i<width; i++)
      {
	for(ElementIndex j=0; j<mesh.GetNE(); j++)
	  {
	    if(!working_elements.Test(j))
	      {  
		const Element & el = mesh[j];
		bool set_active = false;
		
		for(int k=1; !set_active && k<=el.GetNP(); k++)
		  set_active = working_points.Test(el.PNum(k));
		
		if(set_active)
		  working_elements.Set(j);
	      }
	  }

	for(ElementIndex j=0; j<mesh.GetNE(); j++)
	  {
	    if(working_elements.Test(j))
	      {
		const Element & el = mesh[j];
		for(int k=1; k<=el.GetNP(); k++)
		  working_points.Set(el.PNum(k));
	      }
	  }
      }
  }



  void RepairBisection(Mesh & mesh, Array<ElementIndex> & bad_elements, 
		       const BitArray & isnewpoint, const Refinement & refinement,
		       const Array<double> & pure_badness, 
		       double max_worsening, const bool uselocalworsening,
		       const Array< Array<int,PointIndex::BASE>* > & idmaps)
  {
    std::ostringstream ostrstr;

    const int maxtrials = 100;

    //bool doit;
    //std::cout << "DOIT: " << flush;
    //cin >> doit;

    int ne = mesh.GetNE();
    int np = mesh.GetNP();

    int numbadneighbours = 3;
    const int numtopimprove = 3;

    PrintMessage(1,"repairing");

    PushStatus("Repair Bisection");

    Array<Point<3>* > should(np);
    Array<Point<3>* > can(np);
    Array<Vec<3>* > nv(np);
    for(int i=0; i<np; i++)
      {
	nv[i] = new Vec<3>;
	should[i] = new Point<3>;
	can[i] = new Point<3>;
      }
    
    BitArray isboundarypoint(np),isedgepoint(np);
    isboundarypoint.Clear();
    isedgepoint.Clear();

    for(int i = 1; i <= mesh.GetNSeg(); i++)
      {
	const Segment & seg = mesh.LineSegment(i);
	isedgepoint.Set(seg[0]);
	isedgepoint.Set(seg[1]);
      }

    Array<int> surfaceindex(np);
    surfaceindex = -1;
    
    for (int i = 1; i <= mesh.GetNSE(); i++)
      {
	const Element2d & sel = mesh.SurfaceElement(i);
	for (int j = 1; j <= sel.GetNP(); j++)
	  if(!isedgepoint.Test(sel.PNum(j)))
	    {
	      isboundarypoint.Set(sel.PNum(j));
	      surfaceindex[sel.PNum(j) - PointIndex::BASE] = 
		mesh.GetFaceDescriptor(sel.GetIndex()).SurfNr();
	    }
      }



    Validate(mesh,bad_elements,pure_badness,
	     ((uselocalworsening) ?  (0.8*(max_worsening-1.) + 1.) : (0.1*(max_worsening-1.) + 1.)),
	     uselocalworsening); // -> larger working area
    BitArray working_elements(ne);
    BitArray working_points(np);

    GetWorkingArea(working_elements,working_points,mesh,bad_elements,numbadneighbours);
    //working_elements.Set();
    //working_points.Set();

    ostrstr.str("");
    ostrstr << "worsening: " <<
      Validate(mesh,bad_elements,pure_badness,max_worsening,uselocalworsening);
    PrintMessage(4,ostrstr.str());

    

    int auxnum=0;
    for(int i=1; i<=np; i++)
      if(working_points.Test(i))
	auxnum++;
    
    ostrstr.str("");
    ostrstr << "Percentage working points: " << 100.*double(auxnum)/np;
    PrintMessage(5,ostrstr.str());
    

    BitArray isworkingboundary(np);
    for(int i=1; i<=np; i++)
      if(working_points.Test(i) && isboundarypoint.Test(i))
	isworkingboundary.Set(i);
      else
	isworkingboundary.Clear(i);


    for(int i=0; i<np; i++)
      *should[i] = mesh.Point(i+1);

    
    for(int i=0; i<np; i++)
      {
	if(isnewpoint.Test(i+PointIndex::BASE) && 
	   //working_points.Test(i+PointIndex::BASE) && 
	   mesh.mlbetweennodes[i+PointIndex::BASE][0] > 0)
	  *can[i] = Center(*can[mesh.mlbetweennodes[i+PointIndex::BASE][0]-PointIndex::BASE],
			   *can[mesh.mlbetweennodes[i+PointIndex::BASE][1]-PointIndex::BASE]);
	else
	  *can[i] = mesh.Point(i+1);
      }


    int cnttrials = 1;
    
    double lamedge = 0.5;
    double lamface = 0.5;
    
    double facokedge = 0;
    double facokface = 0;
    double factryedge;
    double factryface = 0;

    double oldlamedge,oldlamface;

    MeshOptimize2d * optimizer2d = refinement.Get2dOptimizer();
    if(!optimizer2d)
      {
	std::cerr << "No 2D Optimizer!" <<std::endl;
	return;
      }    

    while ((facokedge < 1.-1e-8 || facokface < 1.-1e-8) && 
	   cnttrials < maxtrials &&
	   multithread.terminate != 1)
      {
	std::cerr << "   facokedge " << facokedge << " facokface " << facokface << " cnttrials " << cnttrials <<std::endl
		   << " perc. " << 95. * max2( min2(facokedge,facokface),
					       double(cnttrials)/double(maxtrials)) <<std::endl;

	ostrstr.str("");
	ostrstr << "max. worsening " << max_worsening;
	PrintMessage(5,ostrstr.str());
	oldlamedge = lamedge;
	lamedge *= 6;
	if (lamedge > 2)
	  lamedge = 2;
	   
	if(1==1 || facokedge < 1.-1e-8)
	  {
	    for(int i=0; i<nv.size(); i++)
	      *nv[i] = Vec<3>(0,0,0);
	    for (int i = 1; i <= mesh.GetNSE(); i++)
	      {
		const Element2d & sel = mesh.SurfaceElement(i);
		Vec<3> auxvec = Cross(mesh.Point(sel.PNum(2))-mesh.Point(sel.PNum(1)),
                                      mesh.Point(sel.PNum(3))-mesh.Point(sel.PNum(1)));
		auxvec.Normalize();
		for (int j = 1; j <= sel.GetNP(); j++)
		  if(!isedgepoint.Test(sel.PNum(j)))
		    *nv[sel.PNum(j) - PointIndex::BASE] += auxvec;
	      }
	    for(int i=0; i<nv.size(); i++)
	      nv[i]->Normalize();
	    
	    
	    do  // move edges
	      {
		lamedge *= 0.5;
		cnttrials++;
		if(cnttrials % 10 == 0)
		  max_worsening *= 1.1;
		
		
		factryedge = lamedge + (1.-lamedge) * facokedge;

		ostrstr.str("");
		ostrstr << "lamedge = " << lamedge << ", trying: " << factryedge;
		PrintMessage(5,ostrstr.str());
		

		for (int i = 1; i <= np; i++)
		  {
		    if (isedgepoint.Test(i))
		      {
			for (int j = 0; j < 3; j++)
			  mesh.Point(i)(j) = 
			    lamedge * (*should.Get(i))(j) +
			    (1.-lamedge) * (*can.Get(i))(j);
		      }
		    else
		      mesh.Point(i) = *can.Get(i);
		  }
		if(facokedge < 1.-1e-8)
		  {
		    ostrstr.str("");
		    ostrstr << "worsening: " <<
		      Validate(mesh,bad_elements,pure_badness,max_worsening,uselocalworsening);

		    PrintMessage(5,ostrstr.str());
		  }
		else
		  Validate(mesh,bad_elements,pure_badness,-1,uselocalworsening);


		ostrstr.str("");
		ostrstr << bad_elements.size() << " bad elements";
		PrintMessage(5,ostrstr.str());
	      }
	    while (bad_elements.size() > 0 && 
		   cnttrials < maxtrials &&
		   multithread.terminate != 1);
	  }

	if(cnttrials < maxtrials &&
	   multithread.terminate != 1)
	  {
	    facokedge = factryedge;
	    
	    // smooth faces
	    mesh.CalcSurfacesOfNode();
	    
	    MeshingParameters dummymp;
	    mesh.ImproveMeshJacobianOnSurface(dummymp,isworkingboundary,nv,OPT_QUALITY, &idmaps);
	    
	    for (int i = 1; i <= np; i++)
	      *can.Elem(i) = mesh.Point(i);
	    
	    if(optimizer2d)
	      optimizer2d->ProjectBoundaryPoints(surfaceindex,can,should);
	  }


	oldlamface = lamface;
	lamface *= 6;
	if (lamface > 2)
	  lamface = 2;


	if(cnttrials < maxtrials &&
	   multithread.terminate != 1)
	  {

	    do  // move faces
	      {
		lamface *= 0.5;
		cnttrials++;
		if(cnttrials % 10 == 0)
		  max_worsening *= 1.1;
		factryface = lamface + (1.-lamface) * facokface;

		ostrstr.str("");
		ostrstr << "lamface = " << lamface << ", trying: " << factryface;
		PrintMessage(5,ostrstr.str());
		
		
		for (int i = 1; i <= np; i++)
		  {
		    if (isboundarypoint.Test(i))
		      {
			for (int j = 0; j < 3; j++)
			  mesh.Point(i)(j) = 
			    lamface * (*should.Get(i))(j) +
			    (1.-lamface) * (*can.Get(i))(j);
		      }
		    else
		      mesh.Point(i) = *can.Get(i);
		  }

		ostrstr.str("");
		ostrstr << "worsening: " <<
		  Validate(mesh,bad_elements,pure_badness,max_worsening,uselocalworsening);
		PrintMessage(5,ostrstr.str());
	

		ostrstr.str("");
		ostrstr << bad_elements.size() << " bad elements";
		PrintMessage(5,ostrstr.str());
	      }
	    while (bad_elements.size() > 0 && 
		   cnttrials < maxtrials &&
		   multithread.terminate != 1);
	  }



	if(cnttrials < maxtrials &&
	   multithread.terminate != 1)
	  {
	    facokface = factryface;
	    // smooth interior
	    
	    mesh.CalcSurfacesOfNode();
	    
	    MeshingParameters dummymp;
	    mesh.ImproveMeshJacobian (dummymp, OPT_QUALITY,&working_points);
	    //mesh.ImproveMeshJacobian (OPT_WORSTCASE,&working_points);
	  

	    for (int i = 1; i <= np; i++)
	      *can.Elem(i) = mesh.Point(i);
	  }
	  
	//!
	if((facokedge < 1.-1e-8 || facokface < 1.-1e-8) && 
	   cnttrials < maxtrials &&
	   multithread.terminate != 1)
	  {
	    MeshingParameters dummymp;
	    MeshOptimize3d optmesh(dummymp);
	    for(int i=0; i<numtopimprove; i++)
	      {
		optmesh.SwapImproveSurface(mesh,OPT_QUALITY,&working_elements,&idmaps);
		optmesh.SwapImprove(mesh,OPT_QUALITY,&working_elements);
		
	      }	    

	    //	    mesh.mglevels = 1;
	    
		
	    ne = mesh.GetNE();
	    working_elements.SetSize(ne);
	    
	    
	    for (int i = 1; i <= np; i++)
	      mesh.Point(i) = *should.Elem(i);
	    
	    Validate(mesh,bad_elements,pure_badness,
		     ((uselocalworsening) ?  (0.8*(max_worsening-1.) + 1.) : (0.1*(max_worsening-1.) + 1.)),
		     uselocalworsening);
	    
	    if(lamedge < oldlamedge || lamface < oldlamface)
	      numbadneighbours++;
	    GetWorkingArea(working_elements,working_points,mesh,bad_elements,numbadneighbours);
	    for(int i=1; i<=np; i++)
	      if(working_points.Test(i) && isboundarypoint.Test(i))
		isworkingboundary.Set(i);
	      else
		isworkingboundary.Clear(i);
	    auxnum=0;
	    for(int i=1; i<=np; i++)
	      if(working_points.Test(i))
		auxnum++;

	    
	    ostrstr.str("");
	    ostrstr << "Percentage working points: " << 100.*double(auxnum)/np;
	    PrintMessage(5,ostrstr.str());
	    
	    for (int i = 1; i <= np; i++)
	      mesh.Point(i) = *can.Elem(i);
	  }
	//!

      }

    MeshingParameters dummymp;
    MeshOptimize3d optmesh(dummymp);
    for(int i=0; i<numtopimprove && multithread.terminate != 1; i++)
      {
	optmesh.SwapImproveSurface(mesh,OPT_QUALITY,NULL,&idmaps);
	optmesh.SwapImprove(mesh,OPT_QUALITY);
	//mesh.UpdateTopology();
      }
    mesh.UpdateTopology();
    /*
    if(cnttrials < 100)
      {
	nv = Vec3d(0,0,0);
	for (int i = 1; i <= mesh.GetNSE(); i++)
	  {
	    const Element2d & sel = mesh.SurfaceElement(i);
	    Vec3d auxvec = Cross(mesh.Point(sel.PNum(2))-mesh.Point(sel.PNum(1)),
				 mesh.Point(sel.PNum(3))-mesh.Point(sel.PNum(1)));
	    auxvec.Normalize();
	    for (int j = 1; j <= sel.GetNP(); j++)
	      if(!isedgepoint.Test(sel.PNum(j)))
		nv[sel.PNum(j) - PointIndex::BASE] += auxvec;
	  }
	for(int i=0; i<nv.Size(); i++)
	  nv[i].Normalize();
	

	mesh.ImproveMeshJacobianOnSurface(isboundarypoint,nv,OPT_QUALITY);
	mesh.CalcSurfacesOfNode();
	    // smooth interior
	    
	
	for (int i = 1; i <= np; i++)
	  if(isboundarypoint.Test(i))
	    can.Elem(i) = mesh.Point(i);
	    
	if(optimizer2d)
	  optimizer2d->ProjectBoundaryPoints(surfaceindex,can,should);

	
	for (int i = 1; i <= np; i++)
	  if(isboundarypoint.Test(i))
	    for(int j=1; j<=3; j++)
	      mesh.Point(i).X(j) = should.Get(i).X(j);
      }
    */


    if(cnttrials == maxtrials)
      {
	for (int i = 1; i <= np; i++)
	  mesh.Point(i) = *should.Get(i);

	Validate(mesh,bad_elements,pure_badness,max_worsening,uselocalworsening);
	
	for(int i=0; i<bad_elements.size(); i++)
	  {
	    ostrstr.str("");
	    ostrstr << "bad element:" <<std::endl
		    << mesh[bad_elements[i]][0] << ": " << mesh.Point(mesh[bad_elements[i]][0]) <<std::endl
		    << mesh[bad_elements[i]][1] << ": " << mesh.Point(mesh[bad_elements[i]][1]) <<std::endl
		    << mesh[bad_elements[i]][2] << ": " << mesh.Point(mesh[bad_elements[i]][2]) <<std::endl
		    << mesh[bad_elements[i]][3] << ": " << mesh.Point(mesh[bad_elements[i]][3]);
	    PrintMessage(5,ostrstr.str());
	  }
	for (int i = 1; i <= np; i++)
	  mesh.Point(i) = *can.Get(i);
      }

    for(int i=0; i<np; i++)
      {
	delete nv[i];
	delete can[i];
	delete should[i];
      }

    PopStatus();
  }
}
