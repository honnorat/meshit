//#include <mystdlib.h>


//#include <myadt.hpp>
// class DenseMatrix;
//#include <gprim.hpp>
#include <cstring>
#include "adtree.hpp"

namespace meshit
{


  /* ******************************* ADTree ******************************* */


  ADTreeNode :: ADTreeNode(int adim)
  {
    pi = -1;

    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
    dim = adim;
    data = new double [dim];
    boxmin = NULL;
    boxmax = NULL;
  }




  ADTreeNode :: ~ADTreeNode()
  {
    delete data;
  }


  ADTree :: ADTree (int adim, const double * acmin, 
		    const double * acmax)
    : ela(0), stack(1000), stackdir(1000)
  {
    dim = adim;
    cmin = new double [dim];
    cmax = new double [dim];
    memcpy (cmin, acmin, dim * sizeof(double));
    memcpy (cmax, acmax, dim * sizeof(double));

    root = new ADTreeNode (dim);
    root->sep = (cmin[0] + cmax[0]) / 2;
    root->boxmin = new double [dim];
    root->boxmax = new double [dim];
    memcpy (root->boxmin, cmin, dim * sizeof(double));
    memcpy (root->boxmax, cmax, dim * sizeof(double));
  }

  ADTree :: ~ADTree ()
  {
    ;
  }

  void ADTree :: Insert (const double * p, int pi)
  {
    ADTreeNode *node(NULL);
    ADTreeNode *next;
    int dir;
    int lr(1);

    double * bmin = new double [dim];
    double * bmax = new double [dim];
  
    memcpy (bmin, cmin, dim * sizeof(double));
    memcpy (bmax, cmax, dim * sizeof(double));


    next = root;
    dir = 0;
    while (next)
      {
	node = next;

	if (node->pi == -1)
	  {    
	    memcpy (node->data, p, dim * sizeof(double));
	    node->pi = pi;

	    if (ela.size() < pi+1)
	      ela.resize (pi+1);
	    ela[pi] = node;

	    return;
	  }

	if (node->sep > p[dir])
	  {
	    next = node->left;
	    bmax[dir] = node->sep;
	    lr = 0;
	  }
	else
	  {
	    next = node->right;
	    bmin[dir] = node->sep;
	    lr = 1;
	  }

	dir++;
	if (dir == dim)
	  dir = 0;
      }


    next = new ADTreeNode(dim);
    memcpy (next->data, p, dim * sizeof(double));
    next->pi = pi;
    next->sep = (bmin[dir] + bmax[dir]) / 2;
    next->boxmin = bmin;
    next->boxmax = bmax;

    if (ela.size() < pi+1)
      ela.resize (pi+1);
    ela[pi] = next;


    if (lr)
      node->right = next;
    else
      node->left = next;
    next -> father = node;

    while (node)
      {
	node->nchilds++;
	node = node->father;
      }
  }

  void ADTree :: DeleteElement (int pi)
  {
    ADTreeNode * node = ela[pi];

    node->pi = -1;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }


  void ADTree :: SetCriterion (ADTreeCriterion & acriterion)
  {
    criterion = & acriterion;
  }


  void ADTree :: Reset ()
  {
    stack.Elem(1) = root;
    stackdir.Elem(1) = 0;
    stackindex = 1;
  }


  int ADTree:: Next ()
  {
    ADTreeNode *node;
    int dir;

    if (stackindex == 0)
      return -1;

    do 
      {
	node = stack.Get(stackindex);
	dir = stackdir.Get(stackindex);
	stackindex --;

	if (criterion -> Eval(node))
	  {
	    int ndir = dir + 1;
	    if (ndir == dim)
	      ndir = 0;

	    if (node -> left && criterion -> Eval (node->left))
	      {
		stackindex ++;
		stack.Elem(stackindex) = node -> left;
		stackdir.Elem(stackindex) = ndir;
	      }
	    if (node->right && criterion -> Eval (node -> right))
	      {
		stackindex++;
		stack.Elem(stackindex) = node->right;
		stackdir.Elem(stackindex) = ndir;
	      }
	  
	    if (node -> pi != -1)
	      return node->pi;
	  }
      }
    while (stackindex > 0);

    return -1;
  }


  void ADTree :: GetMatch (Array <int> & matches)
  {
    int nodenr;

    Reset();

    while ( (nodenr = Next()) != -1)
      matches.push_back (nodenr);
  }


  void ADTree :: PrintRec (std::ostream & ost, const ADTreeNode * node) const
  {
  
    if (node->data)
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (int i = 0; i < dim; i++)
	  ost << node->data[i] << " ";
	ost << std::endl;
      }
    if (node->left)
      {
	ost << "l ";
	PrintRec (ost, node->left);
      }
    if (node->right)
      {
	ost << "r ";
	PrintRec (ost, node->right);
      }
  }


  /* ******************************* ADTree3 ******************************* */


  ADTreeNode3 :: ADTreeNode3()
  {
    pi = -1;

    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
  }

  void ADTreeNode3 :: DeleteChilds ()
  {
    if (left)
      {
	left->DeleteChilds();
	delete left;
	left = NULL;
      }
    if (right)
      {
	right->DeleteChilds();
	delete right;
	right = NULL;
      }
  }


  BlockAllocator ADTreeNode3 :: ball(sizeof (ADTreeNode3));


  void * ADTreeNode3 :: operator new(size_t s)
  {
    return ball.Alloc();
  }

  void ADTreeNode3 :: operator delete (void * p)
  {
    ball.Free (p);
  }







  ADTree3 :: ADTree3 (const double * acmin, 
		      const double * acmax)
    : ela(0)
  {
    memcpy (cmin, acmin, 3 * sizeof(double));
    memcpy (cmax, acmax, 3 * sizeof(double));

    root = new ADTreeNode3;
    root->sep = (cmin[0] + cmax[0]) / 2;
  }

  ADTree3 :: ~ADTree3 ()
  {
    root->DeleteChilds();
    delete root;
  }


  void ADTree3 :: Insert (const double * p, int pi)
  {
    ADTreeNode3 *node(NULL);
    ADTreeNode3 *next;
    int dir;
    int lr(0);

    double bmin[3];
    double bmax[3];
  
    memcpy (bmin, cmin, 3 * sizeof(double));
    memcpy (bmax, cmax, 3 * sizeof(double));

    next = root;
    dir = 0;
    while (next)
      {
	node = next;

	if (node->pi == -1)
	  {    
	    memcpy (node->data, p, 3 * sizeof(double));
	    node->pi = pi;

	    if (ela.size() < pi+1)
	      ela.resize (pi+1);
	    ela[pi] = node;

	    return;
	  }

	if (node->sep > p[dir])
	  {
	    next = node->left;
	    bmax[dir] = node->sep;
	    lr = 0;
	  }
	else
	  {
	    next = node->right;
	    bmin[dir] = node->sep;
	    lr = 1;
	  }

	dir++;
	if (dir == 3)
	  dir = 0;
      }


    next = new ADTreeNode3;
    memcpy (next->data, p, 3 * sizeof(double));
    next->pi = pi;
    next->sep = (bmin[dir] + bmax[dir]) / 2;


    if (ela.size() < pi+1)
      ela.resize (pi+1);
    ela[pi] = next;		


    if (lr)
      node->right = next;
    else
      node->left = next;
    next -> father = node;

    while (node)
      {
	node->nchilds++;
	node = node->father;
      }
  }

  void ADTree3 :: DeleteElement (int pi)
  {
    ADTreeNode3 * node = ela[pi];

    node->pi = -1;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }

  void ADTree3 :: GetIntersecting (const double * bmin, 
				   const double * bmax,
				   Array<int> & pis) const
  {
    static Array<ADTreeNode3*> stack(1000);
    static Array<int> stackdir(1000);
    ADTreeNode3 * node;
    int dir, stacks;

    stack.resize (1000);
    stackdir.resize(1000);
    pis.resize(0);

    stack.Elem(1) = root;
    stackdir.Elem(1) = 0;
    stacks = 1;

    while (stacks)
      {
	node = stack.Get(stacks);
	dir = stackdir.Get(stacks); 
	stacks--;

	if (node->pi != -1)
	  {
	    if (node->data[0] >= bmin[0] && node->data[0] <= bmax[0] &&
		node->data[1] >= bmin[1] && node->data[1] <= bmax[1] &&
		node->data[2] >= bmin[2] && node->data[2] <= bmax[2])

	      pis.push_back (node->pi);
	  }


	int ndir = dir+1;
	if (ndir == 3)
	  ndir = 0;

	if (node->left && bmin[dir] <= node->sep)
	  {
	    stacks++;
	    stack.Elem(stacks) = node->left;
	    stackdir.Elem(stacks) = ndir;
	  }
	if (node->right && bmax[dir] >= node->sep)
	  {
	    stacks++;
	    stack.Elem(stacks) = node->right;
	    stackdir.Elem(stacks) = ndir;
	  }
      }
  }

  void ADTree3 :: PrintRec (std::ostream & ost, const ADTreeNode3 * node) const
  {
  
    if (node->data)
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (int i = 0; i < 3; i++)
	  ost << node->data[i] << " ";
	ost << std::endl;
      }
    if (node->left)
      PrintRec (ost, node->left);
    if (node->right)
      PrintRec (ost, node->right);
  }


  /* ******************************* ADTree6 ******************************* */


  ADTreeNode6 :: ADTreeNode6()
  {
    pi = -1;

    left = NULL;
    right = NULL;
    father = NULL;
    nchilds = 0;
  }

  void ADTreeNode6 :: DeleteChilds ()
  {
    if (left)
      {
	left->DeleteChilds();
	delete left;
	left = NULL;
      }
    if (right)
      {
	right->DeleteChilds();
	delete right;
	right = NULL;
      }
  }


  BlockAllocator ADTreeNode6 :: ball (sizeof (ADTreeNode6));
  void * ADTreeNode6 :: operator new(size_t s)
  {
    return ball.Alloc();
  }

  void ADTreeNode6 :: operator delete (void * p)
  {
    ball.Free (p);
  }





  ADTree6 :: ADTree6 (const double * acmin, 
		      const double * acmax)
    : ela(0)
  {
    memcpy (cmin, acmin, 6 * sizeof(double));
    memcpy (cmax, acmax, 6 * sizeof(double));

    root = new ADTreeNode6;
    root->sep = (cmin[0] + cmax[0]) / 2;
  }

  ADTree6 :: ~ADTree6 ()
  {
    root->DeleteChilds();
    delete root;
  }

  void ADTree6 :: Insert (const double * p, int pi)
  {
    ADTreeNode6 *node(NULL);
    ADTreeNode6 *next;
    int dir;
    int lr(0);

    double bmin[6];
    double bmax[6];

  
    memcpy (bmin, cmin, 6 * sizeof(double));
    memcpy (bmax, cmax, 6 * sizeof(double));

    next = root;
    dir = 0;
    while (next)
      {
	node = next;

	if (node->pi == -1)
	  {    
	    memcpy (node->data, p, 6 * sizeof(double));
	    node->pi = pi;

	    if (ela.size() < pi+1)
	      ela.resize (pi+1);
	    ela[pi] = node;

	    return;
	  }

	if (node->sep > p[dir])
	  {
	    next = node->left;
	    bmax[dir] = node->sep;
	    lr = 0;
	  }
	else
	  {
	    next = node->right;
	    bmin[dir] = node->sep;
	    lr = 1;
	  }

	dir++;
	if (dir == 6) dir = 0;
      }


    next = new ADTreeNode6;
    memcpy (next->data, p, 6 * sizeof(double));
    next->pi = pi;
    next->sep = (bmin[dir] + bmax[dir]) / 2;

    if (ela.size() < pi+1)
      ela.resize (pi+1);
    ela[pi] = next;

    if (lr)
      node->right = next;
    else
      node->left = next;
    next -> father = node;

    while (node)
      {
	node->nchilds++;
	node = node->father;
      }
  }

  void ADTree6 :: DeleteElement (int pi)
  {
    ADTreeNode6 * node = ela[pi];

    node->pi = -1;

    node = node->father;
    while (node)
      {
	node->nchilds--;
	node = node->father;
      }
  }

  void ADTree6 :: PrintMemInfo (std::ostream & ost) const
  {
    ost << Elements() << " elements a " << sizeof(ADTreeNode6) 
	<< " Bytes = "
	<< Elements() * sizeof(ADTreeNode6) << std::endl;
    ost << "maxind = " << ela.size() << " = " << sizeof(ADTreeNode6*) * ela.size() << " Bytes" << std::endl;
  }



  class inttn6 {
  public:
    int dir;
    ADTreeNode6 * node;
  };




  void ADTree6 :: GetIntersecting (const double * bmin, 
				   const double * bmax,
				   Array<int> & pis) const
  {
    // static Array<inttn6> stack(10000);
    // stack.SetSize (10000);
    ArrayMem<inttn6,10000> stack(10000);
    pis.resize(0);

    stack[0].node = root;
    stack[0].dir = 0;
    int stacks = 0;

    while (stacks >= 0)
      {
	ADTreeNode6 * node = stack[stacks].node;
	int dir = stack[stacks].dir; 

	stacks--;
	if (node->pi != -1)
	  {
	    if (node->data[0] > bmax[0] || 
		node->data[1] > bmax[1] || 
		node->data[2] > bmax[2] || 
		node->data[3] < bmin[3] || 
		node->data[4] < bmin[4] || 
		node->data[5] < bmin[5])
	      ;
	    else
              {
                pis.push_back (node->pi);
              }
	  }

	int ndir = (dir+1) % 6;

	if (node->left && bmin[dir] <= node->sep)
	  {
	    stacks++;
	    stack[stacks].node = node->left;
	    stack[stacks].dir = ndir;
	  }
	if (node->right && bmax[dir] >= node->sep)
	  {
	    stacks++;
	    stack[stacks].node = node->right;
	    stack[stacks].dir = ndir;
	  }
      }
  }

  void ADTree6 :: PrintRec (std::ostream & ost, const ADTreeNode6 * node) const
  {
  
    if (node->data)
      {
	ost << node->pi << ": ";
	ost << node->nchilds << " childs, ";
	for (int i = 0; i < 6; i++)
	  ost << node->data[i] << " ";
	ost << std::endl;
      }
    if (node->left)
      PrintRec (ost, node->left);
    if (node->right)
      PrintRec (ost, node->right);
  }


  int ADTree6 :: DepthRec (const ADTreeNode6 * node) const
  {
    int ldepth = 0;
    int rdepth = 0;

    if (node->left)
      ldepth = DepthRec(node->left);
    if (node->right)
      rdepth = DepthRec(node->right);
    return 1 + std::max(ldepth, rdepth);
  }

  int ADTree6 :: ElementsRec (const ADTreeNode6 * node) const
  {
    int els = 1;
    if (node->left)
      els += ElementsRec(node->left);
    if (node->right)
      els += ElementsRec(node->right);
    return els;
  }



  /* ************************************* Point3dTree ********************** */



  Point3dTree :: Point3dTree (const Point<3> & pmin, const Point<3> & pmax)
  {
    double pmi[3], pma[3];
    for (int i = 0; i < 3; i++)
      {
	pmi[i] = pmin(i);
	pma[i] = pmax(i);
      }
    tree = new ADTree3 (pmi, pma);
  }

  Point3dTree :: ~Point3dTree ()
  {
    delete tree;
  }



  void Point3dTree :: Insert (const Point<3> & p, int pi)
  {
    double pd[3];
    pd[0] = p(0);
    pd[1] = p(1);
    pd[2] = p(2);
    tree->Insert (pd, pi);
  }

  void Point3dTree :: GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, 
				       Array<int> & pis) const
  {
    double pmi[3], pma[3];
    for (int i = 0; i < 3; i++)
      {
	pmi[i] = pmin(i);
	pma[i] = pmax(i);
      }
    tree->GetIntersecting (pmi, pma, pis);
  }









  Box3dTree :: Box3dTree (const Box<3> & abox)
  {
    boxpmin = abox.PMin();
    boxpmax = abox.PMax();
    double tpmin[6], tpmax[6];
    for (int i = 0; i < 3; i++)
      {
	tpmin[i] = tpmin[i+3] = boxpmin(i);
	tpmax[i] = tpmax[i+3] = boxpmax(i);
      }
    tree = new ADTree6 (tpmin, tpmax);
  }

  Box3dTree :: Box3dTree (const Point<3> & apmin, const Point<3> & apmax)
  {
    boxpmin = apmin;
    boxpmax = apmax;
    double tpmin[6], tpmax[6];
    for (int i = 0; i < 3; i++)
      {
	tpmin[i] = tpmin[i+3] = boxpmin(i);
	tpmax[i] = tpmax[i+3] = boxpmax(i);
      }
    tree = new ADTree6 (tpmin, tpmax);
  }

  Box3dTree :: ~Box3dTree ()
  {
    delete tree;
  }

  void Box3dTree :: Insert (const Point<3> & bmin, const Point<3> & bmax, int pi)
  {
    double tp[6];

    for (int i = 0; i < 3; i++)
      {
	tp[i] = bmin(i);
	tp[i+3] = bmax(i);
      }

    tree->Insert (tp, pi);
  }

  void Box3dTree ::GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, 
				    Array<int> & pis) const
  {
    double tpmin[6];
    double tpmax[6];

    for (int i = 0; i < 3; i++)
      {
	tpmin[i] = boxpmin(i);
	tpmax[i] = pmax(i);
      
	tpmin[i+3] = pmin(i);
	tpmax[i+3] = boxpmax(i);
      }

    tree->GetIntersecting (tpmin, tpmax, pis);
  }

}
