/***************************************************************************/
/*                                                                         */
/* Vorlesung Optimierung I, Gfrerer, WS94/95                               */
/* BFGS-Verfahren zur Lösung freier nichtlinearer Optimierungsprobleme     */
/*                                                                         */
/* Programmautor:  Joachim Schöberl                                        */
/* Matrikelnummer: 9155284                                                 */
/*                                                                         */
/***************************************************************************/

#include <meshgen.hpp>

#include "opti.hpp"

namespace netgen
{

void Cholesky (const DenseMatrix & a,
	       DenseMatrix & l, Vector & d)
{
  // Factors   A = L D L^T

  double x;

  int n = a.Height();
  
  //  std::cerr << "a = " << a <<std::endl;

  l = a;

  for (int i = 1; i <= n; i++)
    {
      for (int j = i; j <= n; j++)
	{
	  x = l.Get(i, j);

	  for (int k = 1; k < i; k++)
	    x -= l.Get(i, k) * l.Get(j, k) * d(k-1); 
          
	  if (i == j)
	    {
	      d(i-1) = x;
	    }
	  else
	    {
	      l.Elem(j, i) = x / d(i-1);
	    }
	}
    }

  for (int i = 1; i <= n; i++)
    {
      l.Elem(i, i) = 1;
      for (int j = i+1; j <= n; j++)
	l.Elem(i, j) = 0;
    }

  /*
  // Multiply:
  std::cerr << "multiplied factors: " <<std::endl;
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      {
	x = 0;
	for (k = 1; k <= n; k++)
	  x += l.Get(i, k) * l.Get(j, k) * d.Get(k);
	std::cerr << x << " ";
      }
  std::cerr <<std::endl;
  */
}


void MultLDLt (const DenseMatrix & l, const Vector & d, const Vector & g, Vector & p)
{
  /*
  int i, j, n;
  double val;

  n = l.Height();
  p = g;
  for (i = 1; i <= n; i++)
    {
      val = 0;
      for (j = i; j <= n; j++)
	val += p.Get(j) * l.Get(j, i);
      p.Set(i, val);
    }
  for (i = 1; i <= n; i++)
    p.Elem(i) *= d.Get(i);

  for (i = n; i >= 1; i--)
    {
      val = 0;
      for (j = 1; j <= i; j++)
	val += p.Get(j) * l.Get(i, j);
      p.Set(i, val);
    }
  */



  double val;

  int n = l.Height();
  p = g;
  
  for (int i = 0; i < n; i++)
    {
      val = 0;
      for (int j = i; j < n; j++)
	val += p(j) * l(j, i);
      p(i) = val;
    }

  for (int i = 0; i < n; i++)
    p(i) *= d(i);

  for (int i = n-1; i >= 0; i--)
    {
      val = 0;
      for (int j = 0; j <= i; j++)
	val += p(j) * l(i, j);
      p(i) = val;
    }
}

void SolveLDLt (const DenseMatrix & l, const Vector & d, const Vector & g, Vector & p)
{
  double val;

  int n = l.Height();
  p = g;

  for (int i = 0; i < n; i++)
    {
      val = 0;
      for (int j = 0; j < i; j++)
	val += p(j) * l(i,j);
      p(i) -= val;
    }

  for (int i = 0; i < n; i++)
    p(i) /= d(i);
  
  for (int i = n-1; i >= 0; i--)
    {
      val = 0;
      for (int j = i+1; j < n; j++)
	val += p(j) * l(j, i);
      p(i) -= val;
    }
}

int LDLtUpdate (DenseMatrix & l, Vector & d, double a, const Vector & u)
{
  // Bemerkung: Es wird a aus R erlaubt
  // Rueckgabewert: 0 .. D bleibt positiv definit
  //                1 .. sonst

  int n = l.Height();

  Vector v(n);
  double t, told, xi;

  told = 1;
  v = u;

  for (int j = 1; j <= n; j++)
    {
      t = told + a * (v(j-1)*v(j-1)) / d(j-1);

      if (t <= 0) 
	{
	  std::cout << "update err, t = " << t << std::endl;
	  return 1;
	}

      xi = a * v(j-1) / (d(j-1) * t);

      d(j-1) *= t / told;

      for (int i = j + 1; i <= n; i++)
	{
	  v(i-1) -= v(j-1) * l.Elem(i, j);
	  l.Elem(i, j) += xi * v(i-1);
	}

      told = t;
    }

  return 0;
}


double BFGS (
	     Vector & x,         // i: Startwert
	     // o: Loesung, falls IFAIL = 0
	     const MinFunction & fun,
	     const OptiParameters & par,
	     double eps
	     )


{
  int n = x.Size();
  long it;
  char a1crit, a3acrit;


  Vector d(n), g(n), p(n), temp(n), bs(n), xneu(n), y(n), s(n), x0(n);
  DenseMatrix l(n);
  DenseMatrix hesse(n);

  double /* normg, */ alphahat, hd, fold;
  double a1, a2;
  const double mu1 = 0.1, sigma = 0.1, xi1 = 1, xi2 = 10;
  const double tau = 0.1, tau1 = 0.1, tau2 = 0.6;

  Vector typx(x.Size());      // i: typische Groessenordnung der Komponenten
  double f, f0;
  double typf;               // i: typische Groessenordnung der Loesung
  double fmin = -1e5;           // i: untere Schranke fuer Funktionswert
  //  double eps = 1e-8;            // i: Abbruchschranke fuer relativen Gradienten
  double tauf = 0.1;            // i: Abbruchschranke fuer die relative Aenderung der
                                //    Funktionswerte
  int ifail;                    // o:  0 .. Erfolg
                                //    -1 .. Unterschreitung von fmin
                                //     1 .. kein Erfolg bei Liniensuche
                                //     2 .. Überschreitung von itmax

  typx = par.typx;
  typf = par.typf;


  l = 0;
  for (int i = 1; i <= n; i++)
    l.Elem(i, i) = 1;

  f = fun.FuncGrad (x, g);
  f0 = f;
  x0 = x;

  it = 0;
  do
    {
      // Restart
      // std::cout << "it " << it << "f = " << f <<std::endl;
      if (it % (5 * n) == 0)
	{

	  for (int i = 1; i <= n; i++)
	    d(i-1) = typf/ (typx(i-1)*typx(i-1));   // 1;
	  for (int i = 2; i <= n; i++)
	    for (int j = 1; j < i; j++)
	      l.Elem(i, j) = 0;

	  /*
	  hesse = 0;
	  for (i = 1; i <= n; i++)
	    hesse.Elem(i, i) = typf / sqr (typx.Get(i));  

	  fun.ApproximateHesse (x, hesse);

	  Cholesky (hesse, l, d);
	  */
	}

      it++;
      if (it > par.maxit_bfgs)
	{
	  ifail = 2;
	  break;
	}


      // Solve with factorized B

      SolveLDLt (l, d, g, p);

 //      std::cerr << "l " << l <<std::endl
// 		 << "d " << d <<std::endl
// 		 << "g " << g <<std::endl
// 		 << "p " << p <<std::endl;


      p *= -1;
      y = g;

      fold = f;

      // line search

      alphahat = 1;
      lines (x, xneu, p, f, g, fun, par, alphahat, fmin,
	     mu1, sigma, xi1, xi2, tau, tau1, tau2, ifail);

      if(ifail == 1)
	std::cerr << "no success with linesearch" << std::endl;

       /*
      // if (it > par.maxit_bfgs/2)
	{
	  std::cerr << "x = " << x <<std::endl;
	  std::cerr << "xneu = " << xneu <<std::endl;
	  std::cerr << "f = " << f <<std::endl;
	  std::cerr << "g = " << g <<std::endl;
	}
      */

      //      std::cerr << "it = " << it << " f = " << f <<std::endl;
      //      if (ifail != 0) break;

      s.Set2 (1, xneu, -1, x);
      y *= -1;
      y.Add (1,g); // y += g;

      x = xneu;

      // BFGS Update

      MultLDLt (l, d, s, bs);

      a1 = y * s;
      a2 = s * bs;

      if (a1 > 0 && a2 > 0)
	{
	  if (LDLtUpdate (l, d, 1 / a1, y) != 0)
	    {
              std::cerr << "BFGS update error1" << std::endl;
	      std::cerr << "BFGS update error1" << std::endl;
	      std::cerr << "l " << std::endl << l << std::endl
			 << "d " << d << std::endl;
	      ifail = 1;
	      break;
	    }

	  if (LDLtUpdate (l, d, -1 / a2, bs) != 0)
	    {
              std::cerr << "BFGS update error2" << std::endl;
	      std::cerr << "l " << std::endl << l << std::endl
			<< "d " << d << std::endl;
	      ifail = 1;
	      break;
	    }
	}

      // Calculate stop conditions

      hd = eps * std::max (typf, fabs (f));
      a1crit = 1;
      for (int i = 1; i <= n; i++)
	if ( fabs (g(i-1)) * std::max (typx(i-1), fabs (x(i-1))) > hd)
	  a1crit = 0;


      a3acrit = (fold - f <= tauf * std::max (typf, fabs (f)));

      //    testout << "g = " << g <<std::endl;
      //    testout << "a1crit, a3crit = " << int(a1crit) << ", " << int(a3acrit) <<std::endl;

      /*
	// Output for tests

	normg = sqrt (g * g);

	testout << "it =" << setw (5) << it
	<< " f =" << setw (12) << setprecision (5) << f
	<< " |g| =" << setw (12) << setprecision (5) << normg;

	testout << " x = (" << setw (12) << setprecision (5) << x.Elem(1);
	for (i = 2; i <= n; i++)
	testout << "," << setw (12) << setprecision (5) << x.Elem(i);
	testout << ")" <<std::endl;
	*/

      //std::cerr << "it = " << it << " f = " << f << " x = " << x <<std::endl
      //	 << " g = " << g << " p = " << p <<std::endl <<std::endl;

      //      std::cerr << "|g| = " << g.L2Norm() <<std::endl;

      if (g.L2Norm() < fun.GradStopping (x)) break;

    }
  while (!a1crit || !a3acrit);

  /*
  std::cerr << "it = " << it << " g = " << g << " f = " << f 
	     << " fail = " << ifail <<std::endl;
  */
  if (f0 < f || (ifail == 1))
    {
      std::cerr << "fail, f = " << f << " f0 = " << f0 << std::endl;
      f = f0;
      x = x0;
    }

  // std::cout <<std::endl;

  //  std::cerr << "x = " << x << ", x0 = " << x0 <<std::endl;
  return f;
}

}
