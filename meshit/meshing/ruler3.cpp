#include <climits>
#include <meshit.hpp>
#include "meshing3.hpp"


namespace meshit
{
extern double minother;
extern double minwithoutother;


static double CalcElementBadness (const Array<Point3d> & points,
				  const Element & elem)
{
  double vol, l, l4, l5, l6;
  if (elem.GetNP() != 4) 
    {
      if (elem.GetNP() == 5)
	{
	  double z = points.Get(elem.PNum(5)).Z();
	  if (z > -1e-8) return 1e8;
	  return (-1 / z) - z; //  - 2;
	}
      return 0;
    }
  
  Vec3d v1 = points.Get(elem.PNum(2)) - points.Get(elem.PNum(1));
  Vec3d v2 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(1));
  Vec3d v3 = points.Get(elem.PNum(4)) - points.Get(elem.PNum(1));
  
  vol = - (Cross (v1, v2) * v3);
  l4 = Dist (points.Get(elem.PNum(2)), points.Get(elem.PNum(3)));
  l5 = Dist (points.Get(elem.PNum(2)), points.Get(elem.PNum(4)));
  l6 = Dist (points.Get(elem.PNum(3)), points.Get(elem.PNum(4)));

  l = v1.Length() + v2.Length() + v3.Length() + l4 + l5 + l6;
  
  //  testout << "vol = " << vol << " l = " << l <<std::endl;
  if (vol < 1e-8) return 1e10;
  //  std::cerr << "l^3/vol = " << (l*l*l / vol) <<std::endl;
  
  double err = pow (l*l*l/vol, 1.0/3.0) / 12;
  return err;
}






int Meshing3 :: ApplyRules 
(
 Array<Point3d> & lpoints,     // in: local points, out: old+new local points
 Array<int> & allowpoint,      // in: 2 .. it is allowed to use pointi, 1..will be allowed later, 0..no means
 Array<MiniElement2d> & lfaces,    // in: local faces, out: old+new local faces
 INDEX lfacesplit,	       // for local faces in outer radius
 INDEX_2_HASHTABLE<int> & connectedpairs,  // connected pairs for prism-meshing
 Array<Element> & elements,    // out: new elements
 Array<INDEX> & delfaces,      // out: face indices of faces to delete
 int tolerance,                // quality class: 1 best 
 double sloppy,                // quality strength
 int rotind1,                  // how to rotate base element
 float & retminerr             // element error 
 )

{
  int i, j, k, ri, nfok, npok, incnpok, refpi, locpi, locfi, locfr;
  float hf, err, minerr, teterr, minteterr;
  char ok, found, hc;
  vnetrule * rule;
  Vector oldu, newu, newu1, newu2, allp;
  Vec3d ui;
  Point3d np;
  int oldnp, noldlp, noldlf;
  const MiniElement2d * locface = NULL;

  Array<int> pused;        // point is already mapped
  Array<char> fused;       // face is already mapped
  Array<int> pmap;         // map of reference point to local point
  Array<char> pfixed;      // point mapped by face-map
  Array<int> fmapi;        // face in reference is mapped to face nr ...
  Array<int> fmapr;        // face in reference is rotated to map 
  Array<Point3d> transfreezone;  // transformed free-zone
  INDEX_2_CLOSED_HASHTABLE<int> ledges(100); // edges in local environment
  
  Array<Point3d> tempnewpoints;
  Array<MiniElement2d> tempnewfaces;
  Array<int> tempdelfaces;
  Array<Element> tempelements;
  Array<Box3d> triboxes;         // bounding boxes of local faces

  Array<int, PointIndex::BASE> pnearness;
  Array<int> fnearness;

  static int cnt = 0;
  cnt++;
  
  delfaces.resize (0);
  elements.resize (0);

  // determine topological distance of faces and points to
  // base element

  pnearness.resize (lpoints.size());
  fnearness.resize (lfacesplit);

  pnearness = INT_MAX/10;
  for (j = 0; j < lfaces[0].GetNP(); j++)
    pnearness[lfaces[0][j]] = 0;

  for (int loop = 0; loop < 2; loop++)
    {

      for (i = 0; i < lfacesplit; i++)
	{
	  const MiniElement2d & hface = lfaces[i];

	  int minn = INT_MAX-1;
	  for (j = 0; j < hface.GetNP(); j++)
	    {
	      int hi = pnearness[hface[j]];
	      if (hi < minn) minn = hi;
	    }
	  if (minn < INT_MAX/10)
	    for (j = 0; j < hface.GetNP(); j++)
	      if (pnearness[hface[j]] > minn+1)
		pnearness[hface[j]] = minn+1;
	}

      for (i = 1; i <= connectedpairs.GetNBags(); i++)
	for (j = 1; j <= connectedpairs.GetBagSize(i); j++)
	  {
	    INDEX_2 edge;
	    int val;
	    connectedpairs.GetData (i, j, edge, val);

	    if (pnearness[edge.I1()] > pnearness[edge.I2()] + 1)
	      pnearness[edge.I1()] = pnearness[edge.I2()] + 1;

	    if (pnearness[edge.I2()] > pnearness[edge.I1()] + 1)
	      pnearness[edge.I2()] = pnearness[edge.I1()] + 1;
	  }

    }

  for (i = 0; i < fnearness.size(); i++)
    {
      int sum = 0;
      for (j = 0; j < lfaces[i].GetNP(); j++)
	sum += pnearness[lfaces[i][j]];
      fnearness[i] = sum;
    }

  
  // find bounding boxes of faces

  triboxes.resize (lfaces.size());
  for (i = 0; i < lfaces.size(); i++)
    {
      const MiniElement2d & face = lfaces[i];
      triboxes[i].SetPoint (lpoints.Get(face[0]));
      for (j = 1; j < face.GetNP(); j++)
	triboxes[i].AddPoint (lpoints.Get(face[j]));
    }
  
  bool useedges = 0;
  for (ri = 0; ri < rules.size(); ri++)
    if (rules[ri]->GetNEd()) useedges = 1;

  if (useedges)
    {
      ledges.SetSize (5 * lfacesplit);
      
      for (j = 0; j < lfacesplit; j++)
	// if (fnearness[j] <= 5) 
	  {
	    const MiniElement2d & face = lfaces[j];
	    int newp, oldp;
	    
	    newp = face[face.GetNP()-1];
	    for (k = 0; k < face.GetNP(); k++)
	      {
		oldp = newp;
		newp = face[k];
		ledges.Set (INDEX_2::Sort(oldp, newp), 1);
	      }
	  }
    }

  pused.resize (lpoints.size());
  fused.resize (lfaces.size());

  found = 0;
  minerr = tolfak * tolerance * tolerance;
  minteterr = sloppy * tolerance;

  // impossible, if no rule can be applied at any tolerance class
  bool impossible = 1;


  // check each rule:

  for (ri = 1; ri <= rules.size(); ri++)
    { 
      // sprintf (problems.Elem(ri), "");
      *problems.Elem(ri) = '\0';

      rule = rules.Get(ri);
      
      if (rule->GetNP(1) != lfaces[0].GetNP())
	continue;

      if (rule->GetQuality() > tolerance)
	{
	  if (rule->GetQuality() < 100) impossible = 0;

	  continue;
	}
      
      pmap.resize (rule->GetNP());
      fmapi.resize (rule->GetNF());
      fmapr.resize (rule->GetNF());
      
      fused = 0;
      pused = 0;
      pmap = 0;
      fmapi = 0;
      for (i = 1; i <= fmapr.size(); i++)
	fmapr.Set(i, rule->GetNP(i));
      
      fused[0] = 1;
      fmapi[0] = 1;
      fmapr[0] = rotind1;

      
      for (j = 1; j <= lfaces.Get(1).GetNP(); j++)
	{
	  locpi = lfaces[0].PNumMod (j+rotind1);
	  pmap.Set (rule->GetPointNr (1, j), locpi);
	  pused.Elem(locpi)++;
	}

      /*
	map all faces
	nfok .. first nfok-1 faces are mapped properly
	*/

      nfok = 2;
      while (nfok >= 2)
	{
	  
	  if (nfok <= rule->GetNOldF())
	    {
	      // not all faces mapped

	      ok = 0;
	      locfi = fmapi.Get(nfok);
	      locfr = fmapr.Get(nfok);

	      int actfnp = rule->GetNP(nfok);

	      while (!ok)
		{
		  locfr++;
		  if (locfr == actfnp + 1)
		    {
		      locfr = 1;
		      locfi++;
		      if (locfi > lfacesplit) break;
		    }
		  
		  
		  if (fnearness.Get(locfi) > rule->GetFNearness (nfok) ||
		      fused.Get(locfi) ||
		      actfnp != lfaces.Get(locfi).GetNP() )
		    {
		      // face not feasible in any rotation

		      locfr = actfnp;
		    }
		  else
		    {
		      
		      ok = 1;
		      
		      locface = &lfaces.Get(locfi);

		      
		      // reference point already mapped differently ?
		      for (j = 1; j <= actfnp && ok; j++)
			{
			  locpi = pmap.Get(rule->GetPointNr (nfok, j));
			  
			  if (locpi && locpi != locface->PNumMod(j+locfr))
			    ok = 0;
			}
		      
		      // local point already used or point outside tolerance ?
		      for (j = 1; j <= actfnp && ok; j++)
			{
			  refpi = rule->GetPointNr (nfok, j);
			  
			  if (pmap.Get(refpi) == 0)
			    {
			      locpi = locface->PNumMod (j + locfr);

			      if (pused.Get(locpi))
				ok = 0;
			      else
				{
				  const Point3d & lp = lpoints.Get(locpi);
				  const Point3d & rp = rule->GetPoint(refpi);

				  if ( Dist2 (lp, rp) * rule->PointDistFactor(refpi) > minerr)
				    {
				      impossible = 0;
				      ok = 0;
				    }
				}
			    }
			}
		    }
		}
	      
	      
	      if (ok)
		{
		  // map face nfok

		  fmapi.Set (nfok, locfi);
		  fmapr.Set (nfok, locfr);
		  fused.Set (locfi, 1);
		  
		  for (j = 1; j <= rule->GetNP (nfok); j++)
		    {
		      locpi = locface->PNumMod(j+locfr);
		      
		      if (rule->GetPointNr (nfok, j) <= 3 &&
			  pmap.Get(rule->GetPointNr(nfok, j)) != locpi)
			std::cerr << "change face1 point, mark1" <<std::endl;
		      
		      pmap.Set(rule->GetPointNr (nfok, j), locpi);
		      pused.Elem(locpi)++;
		    }
		  
		  nfok++;
		}
	      else
		{
		  // backtrack one face
		  fmapi.Set (nfok, 0);
		  fmapr.Set (nfok, rule->GetNP(nfok));
		  nfok--;
		  
		  fused.Set (fmapi.Get(nfok), 0);
		  for (j = 1; j <= rule->GetNP (nfok); j++)
		    {
		      refpi = rule->GetPointNr (nfok, j);
		      pused.Elem(pmap.Get(refpi))--;
		      
		      if (pused.Get(pmap.Get(refpi)) == 0)
			{
			  pmap.Set(refpi, 0);
			}
		    }
		}
	    }
	  
	  else
	    
	    { 
	      // all faces are mapped
	      // now map all isolated points:
	      
	      npok = 1;
	      incnpok = 1;
	      
	      pfixed.resize (pmap.size());
	      for (i = 1; i <= pmap.size(); i++)
		pfixed.Set(i, (pmap.Get(i) != 0) );
	      
	      while (npok >= 1)
		{
		  
		  if (npok <= rule->GetNOldP())
		    {
		      
		      if (pfixed.Get(npok))
			
			{
			  if (incnpok)
			    npok++;
			  else
			    npok--;
			}
		      
		      else
			
			{
			  locpi = pmap.Elem(npok);
			  ok = 0;
			  
			  if (locpi)
			    pused.Elem(locpi)--;
			  
			  while (!ok && locpi < lpoints.size())
			    {
			      ok = 1;
			      locpi++;
			      
			      if (pused.Get(locpi) || 
				  pnearness.Get(locpi) > rule->GetPNearness(npok))
				{
				  ok = 0;
				}
			      else if (allowpoint.Get(locpi) != 2)
				{
				  ok = 0;
				  if (allowpoint.Get(locpi) == 1)
				    impossible = 0;
				}
			      else
				{
				  const Point3d & lp = lpoints.Get(locpi);
				  const Point3d & rp = rule->GetPoint(npok);

				  if ( Dist2 (lp, rp) * rule->PointDistFactor(npok) > minerr)
				    {
				      ok = 0;
				      impossible = 0;
				    }
				}
			    }
			  
			  
			  if (ok)
			    {
			      pmap.Set (npok, locpi);
			      
			      if (npok <= 3)
				std::cerr << "set face1 point, mark3" <<std::endl;
			      
			      pused.Elem(locpi)++;
			      npok++;
			      incnpok = 1;
			    }
			  
			  else
			    
			    {
			      pmap.Set (npok, 0);
			      
			      if (npok <= 3)
				std::cerr << "set face1 point, mark4" <<std::endl;
			      
			      npok--;
			      incnpok = 0;
			    }
			}
		    }
		  
		  else
		    
		    {
		      // all points are mapped
		      
		      ok = 1;
		      
		      
		      // check mapedges:
		      for (i = 1; i <= rule->GetNEd(); i++)
			{
			  INDEX_2 in2(pmap.Get(rule->GetEdge(i).i1),
				      pmap.Get(rule->GetEdge(i).i2));
			  in2.Sort();
			  if (!ledges.Used (in2)) ok = 0;
			}


		      // check prism edges:
		      for (i = 1; i <= rule->GetNE(); i++)
			{
			  const Element & el = rule->GetElement (i);
			  if (el.GetType() == PRISM) 
			    { 
			      for (j = 1; j <= 3; j++)
				{
				  INDEX_2 in2(pmap.Get(el.PNum(j)),
					      pmap.Get(el.PNum(j+3)));      
				  in2.Sort();
				  if (!connectedpairs.Used (in2)) ok = 0;
				}
			    }
			  if (el.GetType() == PYRAMID) 
			    { 
			      for (j = 1; j <= 2; j++)
				{
				  INDEX_2 in2;
				  if (j == 1)
				    {
				      in2.I1() = pmap.Get(el.PNum(2));
				      in2.I2() = pmap.Get(el.PNum(3));
				    }
				  else
				    {
				      in2.I1() = pmap.Get(el.PNum(1));
				      in2.I2() = pmap.Get(el.PNum(4));
				    }
				  in2.Sort();
				  if (!connectedpairs.Used (in2)) 
				    {
				      ok = 0;
				    }
				}
			    }

			}
		      

		      
		      for (i = rule->GetNOldF() + 1; i <= rule->GetNF(); i++)
			fmapi.Set(i, 0);
		      

		      if (ok)
			{
			  foundmap.Elem(ri)++;
			}

		      


		      // deviation of existing points

		      oldu.SetSize (3 * rule->GetNOldP());
		      newu.SetSize (3 * (rule->GetNP() - rule->GetNOldP()));
		      allp.SetSize (3 * rule->GetNP());
		      
		      for (i = 1; i <= rule->GetNOldP(); i++)
			{
			  const Point3d & lp = lpoints.Get(pmap.Get(i));
			  const Point3d & rp = rule->GetPoint(i);
			  oldu (3*i-3) = lp.X()-rp.X();
                          oldu (3*i-2) = lp.Y()-rp.Y();
			  oldu (3*i-1) = lp.Z()-rp.Z();
			  
			  allp (3*i-3) = lp.X();
                          allp (3*i-2) = lp.Y();
                          allp (3*i-1) = lp.Z();
			}

		      if (rule->GetNP() > rule->GetNOldP())
			{
			  newu.SetSize (rule->GetOldUToNewU().Height());
			  rule->GetOldUToNewU().Mult (oldu, newu);
			}

		      //		      int idiff = 3 * (rule->GetNP()-rule->GetNOldP());
		      int idiff = 3 * rule->GetNOldP();
		      for (i = rule->GetNOldP()+1; i <= rule->GetNP(); i++)
			{
			  const Point3d & rp = rule->GetPoint(i);
			  allp (3*i-3) = rp.X() + newu(3*i-3 - idiff);
                          allp (3*i-2) = rp.Y() + newu(3*i-2 - idiff);
                          allp (3*i-1) = rp.Z() + newu(3*i-1 - idiff);
			}
		      
		      rule->SetFreeZoneTransformation (allp, 
						       tolerance + int(sloppy));

		      if (!rule->ConvexFreeZone())
			{
			  ok = 0;
			  sprintf (problems.Elem(ri), "Freezone not convex");

			}

		      // check freezone:
		      
		      for (i = 1; i <= lpoints.size(); i++)
			{
			  if ( !pused.Get(i) )
			    {
			      const Point3d & lp = lpoints.Get(i);

			      if (rule->fzbox.IsIn (lp))
				{
				  if (rule->IsInFreeZone(lp))
				    {
				      ok = 0;
				      break;
				    }
				}
			    }
			}

		      for (i = 1; i <= lfaces.size() && ok; i++)
			{
			  static Array<int> lpi(4);

			  if (!fused.Get(i))
			    { 
			      int triin;
			      const MiniElement2d & lfacei = lfaces.Get(i);

			      if (!triboxes.Elem(i).Intersect (rule->fzbox))
				triin = 0;
			      else
				{
				  int li, lj;
				  for (li = 1; li <= lfacei.GetNP(); li++)
				    {
				      int lpii = 0;
				      int pi = lfacei.PNum(li);
				      for (lj = 1; lj <= rule->GetNOldP(); lj++)
					if (pmap.Get(lj) == pi)
					  lpii = lj;
				      lpi.Elem(li) = lpii;
				    }


				  if (lfacei.GetNP() == 3)
				    {
				      triin = rule->IsTriangleInFreeZone 
					(
					 lpoints.Get(lfacei.PNum(1)),
					 lpoints.Get(lfacei.PNum(2)),
					 lpoints.Get(lfacei.PNum(3)), lpi, 1
					 );
				    }
				  else
				    {
				      triin = rule->IsQuadInFreeZone 
					(
					 lpoints.Get(lfacei.PNum(1)),
					 lpoints.Get(lfacei.PNum(2)),
					 lpoints.Get(lfacei.PNum(3)), 
					 lpoints.Get(lfacei.PNum(4)), 
					 lpi, 1
					 );
				    }
				}


			      if (triin == -1)
				{
				  ok = 0;
				}
			      
			      if (triin == 1)
				{
				  hc = 0;
				  for (k = rule->GetNOldF() + 1; k <= rule->GetNF(); k++)
				    {
				      if (rule->GetPointNr(k, 1) <= rule->GetNOldP() &&
					  rule->GetPointNr(k, 2) <= rule->GetNOldP() &&
					  rule->GetPointNr(k, 3) <= rule->GetNOldP())
					{
					  for (j = 1; j <= 3; j++)
					    if (lfaces.Get(i).PNumMod(j  ) == pmap.Get(rule->GetPointNr(k, 1)) &&
						lfaces.Get(i).PNumMod(j+1) == pmap.Get(rule->GetPointNr(k, 3)) &&
						lfaces.Get(i).PNumMod(j+2) == pmap.Get(rule->GetPointNr(k, 2)))
					      {
						fmapi.Elem(k) = i;
						hc = 1;

						strcpy (problems.Elem(ri), "other");
					      }
					}
				    }
				  
				  if (!hc)
				    {
				      ok = 0;
				    }
				}
			    }
			   
			}

		      
		      if (ok)
			{
			  err = 0;
			  for (i = 1; i <= rule->GetNOldP(); i++)
			    {
			      hf = rule->CalcPointDist (i, lpoints.Get(pmap.Get(i)));
			      if (hf > err) err = hf;
			    }
			  
			  

			  //			  newu = rule->GetOldUToNewU() * oldu;

			  // set new points:
			  
			  oldnp = rule->GetNOldP();
			  noldlp = lpoints.size();
			  noldlf = lfaces.size();
			  
			  
			  for (i = oldnp + 1; i <= rule->GetNP(); i++)
			    {
			      np = rule->GetPoint(i);
			      np.X() += newu (3 * (i-oldnp) - 3);
			      np.Y() += newu (3 * (i-oldnp) - 2);
			      np.Z() += newu (3 * (i-oldnp) - 1);
			      
			      pmap.Elem(i) = lpoints.push_back (np);
			    }
			  
			  // Set new Faces:
			  
			  for (i = rule->GetNOldF() + 1; i <= rule->GetNF(); i++)
			    if (!fmapi.Get(i))
			      {
				MiniElement2d nface(rule->GetNP(i));
				for (j = 1; j <= nface.GetNP(); j++)
				  nface.PNum(j) = pmap.Get(rule->GetPointNr (i, j));
				
				lfaces.push_back (nface);
			      }
			  
			  
			  // Delete old Faces:

			  for (i = 1; i <= rule->GetNDelF(); i++)
			    delfaces.push_back (fmapi.Get(rule->GetDelFace(i)));
			  for (i = rule->GetNOldF()+1; i <= rule->GetNF(); i++)
			    if (fmapi.Get(i))
			      {
				delfaces.push_back (fmapi.Get(i));
				fmapi.Elem(i) = 0;
			      }
			  

			  // check orientation
			  for (i = 1; i <= rule->GetNO() && ok; i++)
			    {
			      const fourint * fouri;
			      
			      fouri = &rule->GetOrientation(i);
			      Vec3d v1 (lpoints.Get(pmap.Get(fouri->i1)), 
					lpoints.Get(pmap.Get(fouri->i2)));
			      Vec3d v2 (lpoints.Get(pmap.Get(fouri->i1)), 
					lpoints.Get(pmap.Get(fouri->i3)));
			      Vec3d v3 (lpoints.Get(pmap.Get(fouri->i1)), 
					lpoints.Get(pmap.Get(fouri->i4)));

			      Vec3d n;
			      Cross (v1, v2, n);
			      //if (n * v3 >= -1e-7*n.Length()*v3.Length()) // OR -1e-7???
			      if (n * v3 >= -1e-9)
				{
				  ok = 0;
				}
			    }

			  

			  // new points in free-zone ?
			  for (i = rule->GetNOldP() + 1; i <= rule->GetNP() && ok; i++)
			    if (!rule->IsInFreeZone (lpoints.Get(pmap.Get(i))))
			      {
				ok = 0;
				
			      }
			  
			  // insert new elements
			  
			  for (i = 1; i <= rule->GetNE(); i++)
			    {
			      elements.push_back (rule->GetElement(i));
			      for (j = 1; j <= elements.Get(i).NP(); j++)
				elements.Elem(i).PNum(j) = pmap.Get(elements.Get(i).PNum(j));
			    }
			  

			  // Calculate Element badness
			  
			  teterr = 0;
			  for (i = 1; i <= elements.size(); i++)
			    {
			      hf = CalcElementBadness (lpoints, elements.Get(i));
			      if (hf > teterr) teterr = hf;
			    }

			  /*
			    // keine gute Erfahrung am 25.1.2000, js
			  if (ok && teterr < 100 &&
			      (rule->TestFlag('b') || tolerance > 10) )
			    {
			      (*mystd::cout) << "Reset teterr " 
				   << rule->Name() 
				   << " err = " << teterr 
				   <<std::endl;
			      teterr = 1;
			    }
			  */

			  // compare edgelength
			  if (rule->TestFlag('l'))
			    {
			      double oldlen = 0;
			      double newlen = 0;

			      for (i = 1; i <= rule->GetNDelF(); i++)
				{
				  const Element2d & face = 
				    rule->GetFace (rule->GetDelFace(i));
				  for (j = 1; j <= 3; j++)
				    {
				      const Point3d & p1 =
					lpoints.Get(pmap.Get(face.PNumMod(j)));
				      const Point3d & p2 =
					lpoints.Get(pmap.Get(face.PNumMod(j+1)));
				      oldlen += Dist(p1, p2);
				    }
				}

			      for (i = rule->GetNOldF()+1; i <= rule->GetNF(); i++)
				{
				  const Element2d & face = rule->GetFace (i);
				  for (j = 1; j <= 3; j++)
				    {
				      const Point3d & p1 =
					lpoints.Get(pmap.Get(face.PNumMod(j)));
				      const Point3d & p2 =
					lpoints.Get(pmap.Get(face.PNumMod(j+1)));
				      newlen += Dist(p1, p2);
				    }
				}

			      if (oldlen < newlen) 
				{
				  ok = 0;
				}
			    }
			  

			  if (ok && teterr < tolerance)
			    {
			      canuse.Elem(ri) ++;
			      /*
			      std::cerr << "can use rule " << rule->Name() 
					 << ", err = " << teterr <<std::endl;
			      for (i = 1; i <= pmap.Size(); i++)
				std::cerr << pmap.Get(i) << " ";
			      std::cerr <<std::endl;
			      */

			      if (strcmp (problems.Elem(ri), "other") == 0)
				{
				  if (teterr < minother)
				    minother = teterr;
				}
			      else
				{
				  if (teterr < minwithoutother)
				    minwithoutother = teterr;
				}
			    }


			  if (teterr > minteterr) impossible = 0;

			  if (ok && teterr < minteterr)
			    {

			      found = ri;
			      minteterr = teterr;
			      
			      tempnewpoints.resize (0);
			      for (i = noldlp+1; i <= lpoints.size(); i++)
				tempnewpoints.push_back (lpoints.Get(i));
			      
			      tempnewfaces.resize (0);
			      for (i = noldlf+1; i <= lfaces.size(); i++)
				tempnewfaces.push_back (lfaces.Get(i));

			      tempdelfaces.resize (0);
			      for (i = 1; i <= delfaces.size(); i++)
				tempdelfaces.push_back (delfaces.Get(i));
			      
			      tempelements.resize (0);
			      for (i = 1; i <= elements.size(); i++)
				tempelements.push_back (elements.Get(i));
			    }
			  

			  lpoints.resize (noldlp);
			  lfaces.resize (noldlf);
			  delfaces.resize (0);
			  elements.resize (0);
			}
		      
		      npok = rule->GetNOldP();
		      incnpok = 0;
		    }
		}
	      
	      nfok = rule->GetNOldF();
	      
	      for (j = 1; j <= rule->GetNP (nfok); j++)
		{
		  refpi = rule->GetPointNr (nfok, j);
		  pused.Elem(pmap.Get(refpi))--;
		  
		  if (pused.Get(pmap.Get(refpi)) == 0)
		    {
		      pmap.Set(refpi, 0);
		    }
		}
	      
	    }
	}
    }

  if (found)
    {
      for (i = 1; i <= tempnewpoints.size(); i++)
	lpoints.push_back (tempnewpoints.Get(i));
      for (i = 1; i <= tempnewfaces.size(); i++)
	if (tempnewfaces.Get(i).PNum(1))
	  lfaces.push_back (tempnewfaces.Get(i));
      for (i = 1; i <= tempdelfaces.size(); i++)
	delfaces.push_back (tempdelfaces.Get(i));
      for (i = 1; i <= tempelements.size(); i++)
	elements.push_back (tempelements.Get(i));
    }
  
  retminerr = minerr;


  if (impossible && found == 0)
    return -1;

  return found;
}
}
