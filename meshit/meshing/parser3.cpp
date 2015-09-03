#include <climits>
#include <fstream>
#include <sstream>
#include <meshit.hpp>
#include "meshing3.hpp"

#ifdef WIN32
#define COMMASIGN ':'
#else
#define COMMASIGN ','
#endif


namespace meshit
{

extern const char * tetrules[];

void LoadVMatrixLine (std::istream & ist, DenseMatrix & m, int line)
{
  char ch;
  int pnum;
  float f;
  
  ist >> ch;
  while (ch != '}')
    {
      ist.putback (ch);
      ist >> f;
      ist >> ch;
      ist >> pnum;
      
      if (ch == 'x' || ch == 'X')
	m.Elem(line, 3 * pnum - 2) = f;
      if (ch == 'y' || ch == 'Y')
	m.Elem(line, 3 * pnum - 1) = f;
      if (ch == 'z' || ch == 'Z')
	m.Elem(line, 3 * pnum    ) = f;

      if (ch == 'p' || ch == 'P')
	{
	  m.Elem(line  , 3 * pnum-2) = f;
	  m.Elem(line+1, 3 * pnum-1) = f;
	  m.Elem(line+2, 3 * pnum  ) = f;
	}

      ist >> ch;
      if (ch == COMMASIGN)
	ist >> ch;
    }
}





int vnetrule :: NeighbourTrianglePoint (const threeint & t1, const threeint & t2) const
{
  Array<int> tr1(3);
  Array<int> tr2(3);
  tr1.Elem(1)=t1.i1;
  tr1.Elem(2)=t1.i2;
  tr1.Elem(3)=t1.i3;
  tr2.Elem(1)=t2.i1;
  tr2.Elem(2)=t2.i2;
  tr2.Elem(3)=t2.i3;


  int ret=0;

  for (int i=1; i<=3; i++)
    {
      for (int j=1; j<=3; j++)
	{
	  if ((tr1.Get(i)==tr2.Get(j) && tr1.Get((i%3)+1)==tr2.Get((j%3)+1)) ||
              (tr1.Get(i)==tr2.Get((j%3)+1) && tr1.Get((i%3)+1)==tr2.Get(j)))
	    {ret = tr2.Get((j+1)%3+1);}
	}      
    }

  return ret;

}

void vnetrule :: LoadRule (std::istream & ist)
{
  char buf[256];
  char ch, ok;
  Point3d p;
  Element2d face;
  int i, j, i1, i2, i3, fs, ii, ii1, ii2, ii3;
  twoint edge;
  DenseMatrix tempoldutonewu(30, 20), 
    tempoldutofreezone(30, 20),
    tempoldutofreezonelimit(30, 20),
    tfz(20, 20),
    tfzl(20, 20);

  tempoldutonewu = 0;
  tempoldutofreezone = 0;
  tfz = 0;
  tfzl = 0;


  noldp = 0;
  noldf = 0;

  ist.get (buf, sizeof(buf), '"');
  ist.get (ch);
  ist.get (buf, sizeof(buf), '"');
  ist.get (ch);

  delete [] name;
  name = new char[strlen (buf) + 1];
  strcpy (name, buf);
  //  (*mystd::cout) << "Rule " << name << " found." <<std::endl;

  do
    {
      ist >> buf;

      if (strcmp (buf, "quality") == 0)

	{
	  ist >> quality;
	}

      else if (strcmp (buf, "flags") == 0)
	{
	  ist >> ch;
	  while (ch != ';')
	    {
	      flags.push_back (ch);
	      ist >> ch;
	    }
	}

      else if (strcmp (buf, "mappoints") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ','
	      ist >> p.Z();
	      ist >> ch;    // ')'

	      points.push_back (p);
	      noldp++;

	      tolerances.resize (noldp);
	      tolerances.Elem(noldp) = 1;

	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      ist >> tolerances.Elem(noldp);
		      ist >> ch;  // '}'
		    }

		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}


      else if (strcmp (buf, "mapfaces") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      face.SetType(TRIG);
	      ist >> face.PNum(1);
	      ist >> ch;    // ','
	      ist >> face.PNum(2);
	      ist >> ch;    // ','
	      ist >> face.PNum(3);
	      ist >> ch;    // ')' or ','
	      if (ch == COMMASIGN)
		{
		  face.SetType(QUAD);
		  ist >> face.PNum(4);
		  ist >> ch;    // ')' 
		}
	      faces.push_back (face);
	      noldf++;

	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == 'd')
		    {
		      delfaces.push_back (noldf);
		      ist >> ch; // 'e'
		      ist >> ch; // 'l'
		    }

		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "mapedges") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> edge.i1;
	      ist >> ch;    // ','
	      ist >> edge.i2;
	      ist >> ch;    // ')'

	      edges.push_back (edge);

	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}


      else if (strcmp (buf, "newpoints") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ','
	      ist >> p.Z();
	      ist >> ch;    // ')'

	      points.push_back (p);

	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      LoadVMatrixLine (ist, tempoldutonewu,
				       3 * (points.size()-noldp) - 2);

		      ist >> ch; // '{'
		      LoadVMatrixLine (ist, tempoldutonewu,
				       3 * (points.size()-noldp) - 1);

		      ist >> ch; // '{'
		      LoadVMatrixLine (ist, tempoldutonewu,
				       3 * (points.size()-noldp)    );
		    }

		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "newfaces") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      face.SetType(TRIG);
	      ist >> face.PNum(1);
	      ist >> ch;    // ','
	      ist >> face.PNum(2);
	      ist >> ch;    // ','
	      ist >> face.PNum(3);
	      ist >> ch;    // ')' or ','
	      if (ch == COMMASIGN)
		{
		  face.SetType(QUAD);
		  ist >> face.PNum(4);
		  ist >> ch;    // ')' 
		}
	      faces.push_back (face);

	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "freezone") == 0)
	{
	  ist >> ch;
	
	  while (ch == '(')
	    {
	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ','
	      ist >> p.Z();
	      ist >> ch;    // ')'
	    
	      freezone.push_back (p);
	    
	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      LoadVMatrixLine (ist, tempoldutofreezone,
				       3 * freezone.size() - 2);
		    
		      ist >> ch; // '{'
		      LoadVMatrixLine (ist, tempoldutofreezone,
				       3 * freezone.size() - 1);
		    
		      ist >> ch; // '{'
		      LoadVMatrixLine (ist, tempoldutofreezone,
				       3 * freezone.size()    );
		    }
		
		  ist >> ch;
		}
	    
	      ist >> ch;
	    }
	
	  ist.putback (ch);
	}
      else if (strcmp (buf, "freezone2") == 0)
	{
	  int k, nfp;

	  nfp = 0;
	  ist >> ch;

	  DenseMatrix hm1(3, 50), hm2(50, 50), hm3(50, 50);
	  hm3 = 0;

	  while (ch == '{')
	    {
	      hm1 = 0;
	      nfp++;
	      LoadVMatrixLine (ist, hm1, 1);

	      for (i = 1; i <= points.size(); i++)
		tfz.Elem(nfp, i) = hm1.Get(1, 3*i-2);


	      p.X() = p.Y() = p.Z() = 0;
	      for (i = 1; i <= points.size(); i++)
		{
		  p.X() += hm1.Get(1, 3*i-2) * points.Get(i).X();
		  p.Y() += hm1.Get(1, 3*i-2) * points.Get(i).Y();
		  p.Z() += hm1.Get(1, 3*i-2) * points.Get(i).Z();
		}
	      freezone.push_back (p);
	      freezonelimit.push_back (p);
	    
	      hm2 = 0;
	      for (i = 1; i <= 3 * noldp; i++)
		hm2.Elem(i, i) = 1;
	      for (i = 1; i <= 3 * noldp; i++)
		for (j = 1; j <= 3 * (points.size() - noldp); j++)
		  hm2.Elem(j + 3 * noldp, i) = tempoldutonewu.Get(j, i);
		  
	      for (i = 1; i <= 3; i++)
		for (j = 1; j <= 3 * noldp; j++)
		  {
		    double sum = 0;
		    for (k = 1; k <= 3 * points.size(); k++)
		      sum += hm1.Get(i, k) * hm2.Get(k, j);
		  
		    hm3.Elem(i + 3 * (nfp-1), j) = sum;
		  }

	      //	    std::cerr << "freepoint: " << p <<std::endl;

	      while (ch != ';')
		ist >> ch; 

	      ist >> ch;
	    }

	  tfzl = tfz;

	  tempoldutofreezone = hm3;
	  tempoldutofreezonelimit = hm3;
	  ist.putback(ch);
	}

      else if (strcmp (buf, "freezonelimit") == 0)
	{
	  int k, nfp;
	  nfp = 0;
	  ist >> ch;

	  DenseMatrix hm1(3, 50), hm2(50, 50), hm3(50, 50);
	  hm3 = 0;

	  while (ch == '{')
	    {
	      hm1 = 0;
	      nfp++;
	      LoadVMatrixLine (ist, hm1, 1);

	      for (i = 1; i <= points.size(); i++)
		tfzl.Elem(nfp, i) = hm1.Get(1, 3*i-2);


	      p.X() = p.Y() = p.Z() = 0;
	      for (i = 1; i <= points.size(); i++)
		{
		  p.X() += hm1.Get(1, 3*i-2) * points.Get(i).X();
		  p.Y() += hm1.Get(1, 3*i-2) * points.Get(i).Y();
		  p.Z() += hm1.Get(1, 3*i-2) * points.Get(i).Z();
		}
	      freezonelimit.Elem(nfp) = p;
	    
	      hm2 = 0;
	      for (i = 1; i <= 3 * noldp; i++)
		hm2.Elem(i, i) = 1;
	      for (i = 1; i <= 3 * noldp; i++)
		for (j = 1; j <= 3 * (points.size() - noldp); j++)
		  hm2.Elem(j + 3 * noldp, i) = tempoldutonewu.Get(j, i);
		  
	      for (i = 1; i <= 3; i++)
		for (j = 1; j <= 3 * noldp; j++)
		  {
		    double sum = 0;
		    for (k = 1; k <= 3 * points.size(); k++)
		      sum += hm1.Get(i, k) * hm2.Get(k, j);
		  
		    hm3.Elem(i + 3 * (nfp-1), j) = sum;
		  }

	      //	    std::cerr << "freepoint: " << p <<std::endl;

	      while (ch != ';')
		ist >> ch; 

	      ist >> ch;
	    }

	  tempoldutofreezonelimit = hm3;
	  ist.putback(ch);
	}

      else if (strcmp (buf, "freeset") == 0)
	{
	  freesets.push_back (new Array<int>);

	  ist >> ch;

	  while (ch != ';')
	    {
	      ist.putback (ch);
	      ist >> i;
	      freesets.Last()->push_back(i);
	      ist >> ch;
	    }
	}

      else if (strcmp (buf, "elements") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      elements.push_back (Element(TET));

	      //	      elements.Last().SetNP(1);
	      ist >> elements.Last().PNum(1);
	      ist >> ch;    // ','

	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(2);
		  ist >> elements.Last().PNum(2);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(3);
		  ist >> elements.Last().PNum(3);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(4);
		  elements.Last().SetType(TET);
		  ist >> elements.Last().PNum(4);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(5);
		  elements.Last().SetType(PYRAMID);
		  ist >> elements.Last().PNum(5);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(6);
		  elements.Last().SetType(PRISM);
		  ist >> elements.Last().PNum(6);
		  ist >> ch;    // ','
		}

	      /*
	      orientations.Append (fourint());
	      orientations.Last().i1 = elements.Last().PNum(1);
	      orientations.Last().i2 = elements.Last().PNum(2);
	      orientations.Last().i3 = elements.Last().PNum(3);
	      orientations.Last().i4 = elements.Last().PNum(4);
	      */

	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "orientations") == 0)

	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      //        fourint a = fourint();
	      orientations.push_back (fourint());

	      ist >> orientations.Last().i1;
	      ist >> ch;    // ','
	      ist >> orientations.Last().i2;
	      ist >> ch;    // ','
	      ist >> orientations.Last().i3;
	      ist >> ch;    // ','
	      ist >> orientations.Last().i4;
	      ist >> ch;    // ','


	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}


      else if (strcmp (buf, "endrule") != 0)
	{
	  PrintSysError ("Parser3d, unknown token " , buf);
	}
    }
  while (!ist.eof() && strcmp (buf, "endrule") != 0);


  //  std::cerr <<std::endl;
  //  std::cerr << Name() <<std::endl;
  //  std::cerr << "no1 = " << GetNO() <<std::endl;

  oldutonewu.SetSize (3 * (points.size() - noldp), 3 * noldp);
  oldutonewu = 0;

  for (i = 1; i <= oldutonewu.Height(); i++)
    for (j = 1; j <= oldutonewu.Width(); j++)
      oldutonewu.Elem(i, j) = tempoldutonewu.Elem(i, j);


  /*
    oldutofreezone = new SparseMatrixFlex (3 * freezone.Size(), 3 * noldp);
    oldutofreezonelimit = new SparseMatrixFlex (3 * freezone.Size(), 3 * noldp);

    oldutofreezone -> SetSymmetric(0);
    oldutofreezonelimit -> SetSymmetric(0);
    */

  /*
    oldutofreezone = new DenseMatrix (3 * freezone.Size(), 3 * noldp);
    oldutofreezonelimit = new DenseMatrix (3 * freezone.Size(), 3 * noldp);
  
    for (i = 1; i <= oldutofreezone->Height(); i++)
    for (j = 1; j <= oldutofreezone->Width(); j++)
    //      if (j == 4 || j >= 7)
    {
    if (tempoldutofreezone.Elem(i, j))
    (*oldutofreezone)(i, j) = tempoldutofreezone(i, j);
    if (tempoldutofreezonelimit.Elem(i, j))
    (*oldutofreezonelimit)(i, j) = tempoldutofreezonelimit(i, j);
    }
    */




  oldutofreezone = new DenseMatrix (freezone.size(), points.size());
  oldutofreezonelimit = new DenseMatrix (freezone.size(), points.size());
  //  oldutofreezone = new SparseMatrixFlex (freezone.Size(), points.Size());
  //  oldutofreezonelimit = new SparseMatrixFlex (freezone.Size(), points.Size());

  for (i = 1; i <= freezone.size(); i++)
    for (j = 1; j <= points.size(); j++)
      {
	if (tfz.Elem(i, j))
	  (*oldutofreezone).Elem(i, j) = tfz.Elem(i, j);
	if (tfzl.Elem(i, j))
	  (*oldutofreezonelimit).Elem(i, j) = tfzl.Elem(i, j);
      }
  
  /*
  std::cerr << "Rule " << Name() <<std::endl;
  std::cerr << "oldutofreezone = " << (*oldutofreezone) <<std::endl;
  std::cerr << "oldutofreezonelimit = " << (*oldutofreezonelimit) <<std::endl;
  */

  freezonepi.resize (freezone.size());
  for (i = 1; i <= freezonepi.size(); i++)
    freezonepi.Elem(i) = 0;
  for (i = 1; i <= freezone.size(); i++)
    for (j = 1; j <= noldp; j++)
      if (Dist (freezone.Get(i), points.Get(j)) < 1e-8)
	freezonepi.Elem(i) = j;



  
  for (i = 1; i <= elements.size(); i++)
    {
      if (elements.Elem(i).GetNP() == 4)
	{
	  orientations.push_back (fourint());
	  orientations.Last().i1 = elements.Get(i).PNum(1);
	  orientations.Last().i2 = elements.Get(i).PNum(2);
	  orientations.Last().i3 = elements.Get(i).PNum(3);
	  orientations.Last().i4 = elements.Get(i).PNum(4);
	}
      if (elements.Elem(i).GetNP() == 5)
	{
	  orientations.push_back (fourint());
	  orientations.Last().i1 = elements.Get(i).PNum(1);
	  orientations.Last().i2 = elements.Get(i).PNum(2);
	  orientations.Last().i3 = elements.Get(i).PNum(3);
	  orientations.Last().i4 = elements.Get(i).PNum(5);

	  orientations.push_back (fourint());
	  orientations.Last().i1 = elements.Get(i).PNum(1);
	  orientations.Last().i2 = elements.Get(i).PNum(3);
	  orientations.Last().i3 = elements.Get(i).PNum(4);
	  orientations.Last().i4 = elements.Get(i).PNum(5);
	}
    }



  if (freesets.size() == 0)
    {
      freesets.push_back (new Array<int>);
      for (i = 1; i <= freezone.size(); i++)
	freesets.Elem(1)->push_back(i);
    }


  //  testout << "Freezone: " <<std::endl;

  //  for (i = 1; i <= freezone.Size(); i++)
  //    std::cerr << "freepoint: " << freezone.Get(i) <<std::endl;
  Vector vp(points.size()), vfp(freezone.size());


  if (quality < 100)
    {
      for (int i = 1; i <= 3; i++)
	{
	  for (int j = 1; j <= points.size(); j++)
	    vp(j-1) = points.Get(j).X(i);
	  oldutofreezone->Mult(vp, vfp);
	  for (int j = 1; j <= freezone.size(); j++)
	    freezone.Elem(j).X(i) = vfp(j-1);
	}
      //      for (i = 1; i <= freezone.Size(); i++)
      //	std::cerr << "freepoint: " << freezone.Get(i) <<std::endl;
    }


  for (fs = 1; fs <= freesets.size(); fs++)
    {
      freefaces.push_back (new Array<threeint>);

      Array<int> & freeset = *freesets.Elem(fs);
      Array<threeint> & freesetfaces = *freefaces.Last();

      for (ii1 = 1; ii1 <= freeset.size(); ii1++)
	for (ii2 = 1; ii2 <= freeset.size(); ii2++)
	  for (ii3 = 1; ii3 <= freeset.size(); ii3++)
	    if (ii1 < ii2 && ii1 < ii3 && ii2 != ii3)
	      {
		i1 = freeset.Get(ii1);
		i2 = freeset.Get(ii2);
		i3 = freeset.Get(ii3);

		Vec3d v1, v2, n;

		v1 = freezone.Get(i3) - freezone.Get(i1);
		v2 = freezone.Get(i2) - freezone.Get(i1);
		n = Cross (v1, v2);
		n /= n.Length();
		//		std::cerr << "i1,2,3 = " << i1 << ", " << i2 << ", " << i3 <<std::endl;
		//		std::cerr << "v1 = " << v1 << " v2 = " << v2 << " n = " << n <<std::endl;
		ok = 1;
		for (ii = 1; ii <= freeset.size(); ii++)
		  {
		    i = freeset.Get(ii);
		    //		    std::cerr << "i = " << i <<std::endl;
		    if (i != i1 && i != i2 && i != i3)
		      if ( (freezone.Get(i) - freezone.Get(i1)) * n < 0 ) ok = 0;
		  }

		if (ok)
		  {
		    freesetfaces.push_back (threeint());
		    freesetfaces.Last().i1 = i1;
		    freesetfaces.Last().i2 = i2;
		    freesetfaces.Last().i3 = i3;
		  }
	      }
    }

  for (fs = 1; fs <= freesets.size(); fs++)
    {
      freefaceinequ.push_back (new DenseMatrix (freefaces.Get(fs)->size(), 4));
    }


  {
    int minn;
    //    Array<int> pnearness (noldp);
    pnearness.resize (noldp);

    for (i = 1; i <= pnearness.size(); i++)
      pnearness.Elem(i) = INT_MAX/10;

    for (j = 1; j <= GetNP(1); j++)
      pnearness.Elem(GetPointNr (1, j)) = 0;

    do
      {
	ok = 1;

	for (i = 1; i <= noldf; i++)
	  {
	    minn = INT_MAX/10;
	    for (j = 1; j <= GetNP(i); j++)
	      minn = min2 (minn, pnearness.Get(GetPointNr (i, j)));

	    for (j = 1; j <= GetNP(i); j++)
	      if (pnearness.Get(GetPointNr (i, j)) > minn+1)
		{
		  ok = 0;
		  pnearness.Elem(GetPointNr (i, j)) = minn+1;
		}
	  }

	for (i = 1; i <= edges.size(); i++)
	  {
	    int pi1 = edges.Get(i).i1;
	    int pi2 = edges.Get(i).i2;

	    if (pnearness.Get(pi1) > pnearness.Get(pi2)+1)
	      {
		ok = 0;
		pnearness.Elem(pi1) = pnearness.Get(pi2)+1;
	      }
	    if (pnearness.Get(pi2) > pnearness.Get(pi1)+1)
	      {
		ok = 0;
		pnearness.Elem(pi2) = pnearness.Get(pi1)+1;
	      }
	  }
	

	for (i = 1; i <= elements.size(); i++)
	  if (elements.Get(i).GetNP() == 6)  // prism rule
	    {
	      for (j = 1; j <= 3; j++)
		{
		  int pi1 = elements.Get(i).PNum(j);
		  int pi2 = elements.Get(i).PNum(j+3);

		  if (pnearness.Get(pi1) > pnearness.Get(pi2)+1)
		    {
		      ok = 0;
		      pnearness.Elem(pi1) = pnearness.Get(pi2)+1;
		    }
		  if (pnearness.Get(pi2) > pnearness.Get(pi1)+1)
		    {
		      ok = 0;
		      pnearness.Elem(pi2) = pnearness.Get(pi1)+1;
		    }
		}
	    }
      }
    while (!ok);

    maxpnearness = 0;
    for (i = 1; i <= pnearness.size(); i++)
      maxpnearness = max2 (maxpnearness, pnearness.Get(i));


    fnearness.resize (noldf);

    for (i = 1; i <= noldf; i++)
      {
	fnearness.Elem(i) = 0;
	for (j = 1; j <= GetNP(i); j++)
	  fnearness.Elem(i) += pnearness.Get(GetPointNr (i, j));
      }

    // std::cerr << "rule " << name << ", pnear = " << pnearness <<std::endl;
  }

  
  //Table of edges:
  for (fs = 1; fs <= freesets.size(); fs++)
    {
      freeedges.push_back (new Array<twoint>);
      
      //      Array<int> & freeset = *freesets.Get(fs);
      Array<twoint> & freesetedges = *freeedges.Last();
      Array<threeint> & freesetfaces = *freefaces.Get(fs);
      int k,l;
      INDEX ind;
      
      for (k = 1; k <= freesetfaces.size(); k++)
	{
          // threeint tr = freesetfaces.Get(k);

	  for (l = k+1; l <= freesetfaces.size(); l++)
	    {
	      ind = NeighbourTrianglePoint(freesetfaces.Get(k), freesetfaces.Get(l));
	      if (!ind) continue;

	      INDEX_3 f1(freesetfaces.Get(k).i1, 
			 freesetfaces.Get(k).i2, 
			 freesetfaces.Get(k).i3);
	      INDEX_3 f2(freesetfaces.Get(l).i1, 
			 freesetfaces.Get(l).i2, 
			 freesetfaces.Get(l).i3);
	      INDEX_2 ed(0, 0);
	      for (int f11 = 1; f11 <= 3; f11++)
		for (int f12 = 1; f12 <= 3; f12++)
		  if (f11 != f12)
		    for (int f21 = 1; f21 <= 3; f21++)
		      for (int f22 = 1; f22 <= 3; f22++)		    
			if (f1.I(f11) == f2.I(f21) && f1.I(f12) == f2.I(f22))
			{
			  ed.I(1) = f1.I(f11);
			  ed.I(2) = f1.I(f12);
			}
	      //	      std::cerr << "ed = " << ed.I(1) << "-" << ed.I(2) <<std::endl;
	      //	      std::cerr << "ind = " << ind << " ed = " << ed <<std::endl;
	      for (int eli = 1; eli <= GetNOldF(); eli++)
		{
		  if (GetNP(eli) == 4)
		    {
		      for (int elr = 1; elr <= 4; elr++)
			{
			  if (GetPointNrMod (eli, elr) == ed.I(1) &&
			      GetPointNrMod (eli, elr+2) == ed.I(2))
			    {
			      /*
			      std::cerr << "ed is diagonal of rectangle" <<std::endl;
			      std::cerr << "ed = " << ed.I(1) << "-" << ed.I(2) <<std::endl;
			      std::cerr << "ind = " << ind <<std::endl;
			      */
			      ind = 0;
			    }

			}
		    }
		}

	      if (ind)
		{
		  /*
		  std::cerr << "new edge from face " << k 
			     << " = (" << freesetfaces.Get(k).i1 
			     << ", " << freesetfaces.Get(k).i2 
			     << ", " << freesetfaces.Get(k).i3
			     << "), point " << ind <<std::endl;
			     */
		  freesetedges.push_back(twoint(k,ind));
		}
	    }	
	}
    }
    
}





void Meshing3 :: LoadRules (const char * filename, const char ** prules)
{
  char buf[256];
  std::istream * ist;
  char *tr1 = NULL;

  if (filename)
    {
      PrintMessage (3, "rule-filename = ", filename);
      ist = new std::ifstream (filename);
    }
  else 
    {
      /* connect tetrules to one string */
      PrintMessage (3, "Use internal rules");
      if (!prules) prules = tetrules;

      const char ** hcp = prules; 
      size_t len = 0;
      while (*hcp)
	{
	  len += strlen (*hcp);
	  hcp++;
	}
      tr1 = new char[len+1];
      tr1[0] = 0;
      hcp = prules; //  tetrules;


      char * tt1 = tr1;
      while (*hcp)
	{
	  strcat (tt1, *hcp);
	  tt1 += strlen (*hcp);	  
	  hcp++;
	}


#ifdef WIN32
      // VC++ 2005 workaround
      for(size_t i=0; i<len; i++)
	if(tr1[i] == ',')
	  tr1[i] = ':';
#endif

      ist = new std::istringstream (tr1);
    }
  
  if (!ist->good())
    {
      std::cerr << "Rule description file " << filename << " not found" <<std::endl;
      delete ist;
      exit (1);
    }
    
  while (!ist->eof())
    {
      buf[0] = 0;
      (*ist) >> buf;
	
      if (strcmp (buf, "rule") == 0)
	{
	  vnetrule * rule = new vnetrule;
	  rule -> LoadRule(*ist);
	  rules.push_back (rule);
	  if (!rule->TestOk())
	    {
	      PrintSysError ("Parser3d: Rule ", rules.size(), " not ok");
	      exit (1);
	    }
	}
      else if (strcmp (buf, "tolfak") == 0)
	{
	  (*ist) >> tolfak;
	}
    }
  delete ist;
  delete [] tr1;
}
}
