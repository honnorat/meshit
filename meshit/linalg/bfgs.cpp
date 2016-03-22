/***************************************************************************/
/*                                                                         */
/* Vorlesung Optimierung I, Gfrerer, WS94/95                               */
/* BFGS-Verfahren zur Lösung freier nichtlinearer Optimierungsprobleme     */
/*                                                                         */
/* Programmautor:  Joachim Schöberl                                        */
/* Matrikelnummer: 9155284                                                 */
/*                                                                         */
/***************************************************************************/

#include <algorithm>

#include "../meshing/msghandler.hpp"
#include "../meshing/smoothing2.hpp"

namespace meshit {

    void MultLDLt(const DenseMatrix& l, const double* d, const double* g, double* p)
    {
        p[0] = g[0] * l(0, 0) + g[1] * l(1, 0);
        p[1] = g[1] * l(1, 1);

        p[0] *= d[0];
        p[1] *= d[1];

        p[1] = p[0] * l(1, 0) + p[1] * l(1, 1);
        p[0] = p[0] * l(0, 0);
    }

    void SolveLDLt(const DenseMatrix& l, const double* d, const double* g, double* p)
    {
        p[0] = g[0];
        p[1] = g[1] - g[0] * l(1, 0);

        p[0] /= d[0];
        p[1] /= d[1];

        p[0] -= p[1] * l(1, 0);
    }

    int LDLtUpdate(DenseMatrix& l, double* d, double a, const double* u)
    {
        double t0 = 1.0 + a * (u[0] * u[0]) / d[0];
        if (t0 <= 0.0) return 1;

        double v1 = u[1] - u[0] * l.Elem(1, 0);
        l.Elem(1, 0) += v1 * a * u[0] / (d[0] * t0);

        double t1 = t0 + a * (v1 * v1) / d[1];
        if (t1 <= 0.0) return 1;

        d[0] *= t0;
        d[1] *= t1 / t0;

        return 0;
    }

    /// Line search, modified Mangasarien conditions
    int lines(double* x,       // i: initial point of line-search
              double* xneu,    // o: solution, if successful
              double* p,       // i: search direction
              double& f,       // i: function-value at x, o: function-value at xneu, iff ifail = 0
              double* g,       // i: gradient at x, o: gradient at xneu, iff ifail = 0
              const MinFunction_2d& fun,     // function to minimize
              const OptiParameters& par,
              double& alphahat,  // i: initial value for alpha_hat, o: solution alpha iff ifail = 0
              double fmin,     // i: lower bound for f
              double mu1,      // i: Parameter mu_1 of Alg.2.1
              double sigma,    // i: Parameter sigma of Alg.2.1
              double xi1,      // i: Parameter xi_1 of Alg.2.1
              double xi2,      // i: Parameter xi_1 of Alg.2.1
              double tau,      // i: Parameter tau of Alg.2.1
              double tau1,     // i: Parameter tau_1 of Alg.2.1
              double tau2)     // i: Parameter tau_2 of Alg.2.1
    // o: 0 on success
    //   -1 bei termination because lower limit fmin
    //    1 bei illegal termination due to different reasons
    {
        constexpr double eps0 = 1E-15;
        double alpha1 = 0;
        double alpha2 = 1e50;
        double phi0 = f;
        double phi1 = f;
        double phi0prime = g[0] * p[0] + g[1] * p[1];
        double phi1prime = phi0prime;

        if (phi0prime > 0) {
            return 1;
        }

        bool flag = true;
        size_t it = 0;
        int ifail = 1;

        while (it++ <= par.maxit_linsearch) {
            xneu[0] = x[0] + alphahat * p[0];
            xneu[1] = x[1] + alphahat * p[1];

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

                double c = (f - phi1 - phi1prime * (alphahat - alpha1)) /
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
                double phihatprime, alphaincr;
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
        return ifail;
    }

    double BFGS_2d(
            double* x,  // i: Startwert, o: Loesung, falls IFAIL = 0
            const MinFunction_2d& fun,
            const OptiParameters& par,
            double eps)
    {
        bool a1crit;
        bool a3acrit;

        double /* normg, */ alphahat, hd, fold;
        double a1, a2;
        const double mu1 = 0.1, sigma = 0.1, xi1 = 1, xi2 = 10;
        const double tau = 0.1, tau1 = 0.1, tau2 = 0.6;

        double typx[2];  // i: typische Groessenordnung der Komponenten
        double f, f0;
        double typf;  // i: typische Groessenordnung der Loesung
        double fmin = -1e5;  // i: untere Schranke fuer Funktionswert
        double tauf = 0.1;  // i: Abbruchschranke fuer die relative Aenderung der
        //    Funktionswerte
        int ifail;   // o:  0 .. Erfolg
        //    -1 .. Unterschreitung von fmin
        //     1 .. kein Erfolg bei Liniensuche
        //     2 .. Überschreitung von itmax

        typx[0] = typx[1] = par.typx;
        typf = par.typf;

        DenseMatrix l(2);
        l.Elem(0, 0) = 1;
        l.Elem(0, 1) = 0;
        l.Elem(1, 0) = 0;
        l.Elem(1, 1) = 1;

        double d[2], g[2], p[2], bs[2], xneu[2], y[2], s[2], x0[2];

        f = fun.FuncGrad(x, g);
        f0 = f;
        x0[0] = x[0];
        x0[1] = x[1];

        size_t it = 0;
        do {
            // Restart
            if (it % (5 * 2) == 0) {
                d[0] = typf / (typx[0] * typx[0]);
                d[1] = typf / (typx[1] * typx[1]);
                for (size_t i = 1; i < 2; i++) {
                    for (size_t j = 0; j < 1; j++) {
                        l.Elem(i, j) = 0;
                    }
                }
            }

            it++;
            if (it > par.maxit_bfgs) {
                ifail = 2;
                break;
            }

            // Solve with factorized B
            SolveLDLt(l, d, g, p);

            p[0] = -p[0];
            p[1] = -p[1];
            y[0] = g[0];
            y[1] = g[1];

            fold = f;

            // line search

            alphahat = 1;
            ifail = lines(x, xneu, p, f, g, fun, par, alphahat, fmin,
                  mu1, sigma, xi1, xi2, tau, tau1, tau2);

            if (ifail == 1)
                std::cerr << "no success with linesearch" << std::endl;

            s[0] = xneu[0] - x[0];
            s[1] = xneu[1] - x[1];

            y[0] = g[0] - y[0];
            y[1] = g[1] - y[1];

            x[0] = xneu[0];
            x[1] = xneu[1];

            // BFGS Update
            MultLDLt(l, d, s, bs);

            a1 = y[0] * s[0] + y[1] * s[1];
            a2 = s[0] * bs[0] + s[1] * bs[1];

            if (a1 > 0 && a2 > 0) {
                if (LDLtUpdate(l, d, 1 / a1, y) != 0) {
                    MESHIT_LOG_ERROR("BFGS update error1");
                    ifail = 1;
                    break;
                }

                if (LDLtUpdate(l, d, -1 / a2, bs) != 0) {
                    MESHIT_LOG_ERROR("BFGS update error2");
                    ifail = 1;
                    break;
                }
            }

            // Calculate stop conditions
            hd = eps * std::max(typf, fabs(f));
            a1crit = true;
            a1crit &= (fabs(g[0]) * std::max(typx[0], fabs(x[0])) <= hd);
            a1crit &= (fabs(g[1]) * std::max(typx[1], fabs(x[1])) <= hd);
            a3acrit = (fold - f <= tauf * std::max(typf, fabs(f)));
        } while (!a1crit || !a3acrit);

        if (f0 < f || (ifail == 1)) {
            MESHIT_LOG_DEBUG("fail, f = " << f << " f0 = " << f0 << " diff = " << f0 - f << "  ifail = " << ifail);
            f = f0;
            x = x0;
        }

        return f;
    }
}  // namespace meshit
