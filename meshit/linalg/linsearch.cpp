/***************************************************************************/
/*                                                                         */
/* Problem:        Liniensuche                                             */
/*                                                                         */
/* Programmautor:  Joachim Schöberl                                        */
/* Matrikelnummer: 9155284                                                 */
/*                                                                         */
/* Algorithmus nach:                                                       */
/*                                                                         */
/*   Optimierung I, Gfrerer, WS94/95                                       */
/*   Algorithmus 2.1: Liniensuche Problem (ii)                             */
/*                                                                         */
/***************************************************************************/


#include <algorithm>

#include "opti.hpp"


namespace meshit {

    const double eps0 = 1E-15;

    double MinFunction::Func(const Vector& /* x */) const
    {
        std::cerr << "Func of MinFunction called" << std::endl;
        return 0;
    }

    void MinFunction::Grad(const Vector& /* x */, Vector& /* g */) const
    {
        std::cerr << "Grad of MinFunction called" << std::endl;
    }

    double MinFunction::FuncGrad(const Vector& x, Vector& g) const
    {
        std::cerr << "Grad of MinFunction called" << std::endl;
        return 0;
    }

    double MinFunction::FuncDeriv(const Vector& x, const Vector& dir, double& deriv) const
    {
        Vector g(x.Size());
        double f = FuncGrad(x, g);
        deriv = (g * dir);

        return f;
    }

    /// Line search, modified Mangasarien conditions
    void lines(Vector& x,       // i: initial point of line-search
               Vector& xneu,    // o: solution, if successful
               Vector& p,       // i: search direction
               double& f,       // i: function-value at x, o: function-value at xneu, iff ifail = 0
               Vector& g,       // i: gradient at x, o: gradient at xneu, iff ifail = 0
               const MinFunction& fun,     // function to minimize
               const OptiParameters& par,
               double& alphahat,  // i: initial value for alpha_hat, o: solution alpha iff ifail = 0
               double fmin,     // i: lower bound for f
               double mu1,      // i: Parameter mu_1 of Alg.2.1
               double sigma,    // i: Parameter sigma of Alg.2.1
               double xi1,      // i: Parameter xi_1 of Alg.2.1
               double xi2,      // i: Parameter xi_1 of Alg.2.1
               double tau,      // i: Parameter tau of Alg.2.1
               double tau1,     // i: Parameter tau_1 of Alg.2.1
               double tau2,     // i: Parameter tau_2 of Alg.2.1
               int& ifail)      // o: 0 on success
    //    -1 bei termination because lower limit fmin
    //     1 bei illegal termination due to different reasons
    {
        double phi0, phi0prime, phi1, phi1prime, phihatprime;
        double alpha1, alpha2, alphaincr, c;
        bool flag = true;
        long it;

        alpha1 = 0;
        alpha2 = 1e50;
        phi0 = phi1 = f;

        phi0prime = g * p;

        if (phi0prime > 0) {
            ifail = 1;
            return;
        }

        ifail = 1;  // Markus

        phi1prime = phi0prime;
        it = 0;

        while (it++ <= par.maxit_linsearch) {
            xneu.Set2(1, x, alphahat, p);

            //    f = fun.FuncGrad (xneu, g);
            f = fun.Func(xneu);

            if (f < fmin) {
                ifail = -1;
                break;
            }

            if (alpha2 - alpha1 < eps0 * alpha2) {
                ifail = 0;
                break;
            }

            if (f - phi0 > mu1 * alphahat * phi1prime + eps0 * fabs(phi0)) {
                flag = false;
                alpha2 = alphahat;

                c = (f - phi1 - phi1prime * (alphahat - alpha1)) /
                    ((alphahat - alpha1) * (alphahat - alpha1));

                alphahat = alpha1 - 0.5 * phi1prime / c;

                if (alphahat > alpha2) {
                    alphahat = alpha1 + 1 / (4 * c) *
                                        ((sigma + mu1) * phi0prime - 2 * phi1prime
                                         + sqrt((phi1prime - mu1 * phi0prime) * (phi1prime - mu1 * phi0prime) -
                                                4 * (phi1 - phi0 - mu1 * alpha1 * phi0prime) * c));
                }
                alphahat = std::max(alphahat, alpha1 + tau * (alpha2 - alpha1));
                alphahat = std::min(alphahat, alpha2 - tau * (alpha2 - alpha1));
            } else {
                /*
                f = fun.FuncGrad (xneu, g);
                phihatprime = g * p;
                 */
                f = fun.FuncDeriv(xneu, p, phihatprime);

                if (phihatprime < sigma * phi0prime * (1 + eps0)) {
                    if (phi1prime < phihatprime) {
                        // Approximationsfunktion ist konvex
                        alphaincr = (alphahat - alpha1) * phihatprime /
                                    (phi1prime - phihatprime);
                    } else {
                        alphaincr = 1e99;  // MAXDOUBLE;
                    }
                    if (flag) {
                        alphaincr = std::max(alphaincr, xi1 * (alphahat - alpha1));
                        alphaincr = std::min(alphaincr, xi2 * (alphahat - alpha1));
                    } else {
                        alphaincr = std::max(alphaincr, tau1 * (alpha2 - alphahat));
                        alphaincr = std::min(alphaincr, tau2 * (alpha2 - alphahat));
                    }
                    alpha1 = alphahat;
                    alphahat += alphaincr;
                    phi1 = f;
                    phi1prime = phihatprime;
                } else {
                    ifail = 0;  // Erfolg !!
                    break;
                }
            }
        }
        fun.FuncGrad(xneu, g);

        if (it < 0)
            ifail = 1;
    }

}  // namespace meshit
