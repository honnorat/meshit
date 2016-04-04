#ifndef FILE_SMOOTHING_HPP
#define FILE_SMOOTHING_HPP

/**************************************************************************/
/* File:   opti.hpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include "meshtype.hpp"
#include "../linalg/densemat.hpp"

namespace meshit
{
    class Opti2dLocalData
    {
     public:
        MeshPoint sp1;
        std::vector<SurfaceElementIndex> loc_elements;
        std::vector<double> lochs;
        std::vector<MeshPoint> loc_pnts2;
        std::vector<MeshPoint> loc_pnts3;
        double loc_metric_weight;

     public:
        Opti2dLocalData()
            : loc_metric_weight{0.0} { }
    };

    class Mesh;

    /**
        Function to be minimized.
     */
    class MinFunction_2d
    {
     protected:
        const Mesh& mesh;
        Opti2dLocalData& ld;

     public:
        MinFunction_2d(const Mesh& amesh, Opti2dLocalData& ald)
            : mesh{amesh}, ld(ald) { }

        double Func(const double* x);

        // function + gradient
        double FuncGrad(const double* x, double* g);

        // directional derivative
        double FuncDeriv(const double* x, const double* dir, double& deriv);
    };

    class OptiParameters
    {
     public:
        size_t maxit_bfgs;
        size_t maxit_linsearch;
        double typf;
        double typx;

        OptiParameters()
            : maxit_bfgs{100},
              maxit_linsearch{100},
              typf{1.0}, typx{1.0} { }
    };


    /** Implementation of BFGS method.
        Efficient method for non-linear minimiztion problems.
        @param x initial value and solution 
        @param fun function to be minimized
     */
    double BFGS_2d(double* x, MinFunction_2d& fun, const OptiParameters& par, double eps = 1e-8);

    int lines(
        double* x,      // i: Ausgangspunkt der Liniensuche
        double* xneu,   // o: Loesung der Liniensuche bei Erfolg
        double* p,      // i: Suchrichtung
        double& f,      // i: Funktionswert an der Stelle x, o: Funktionswert an der Stelle xneu, falls ifail = 0
        double* g,      // i: Gradient an der Stelle x, o: Gradient an der Stelle xneu, falls ifail = 0
        MinFunction_2d& fun,  // function to minmize
        const OptiParameters& par,  // parameters
        double& alphahat,   // i: Startwert für alpha_hat, o: Loesung falls ifail = 0
        double fmin,    // i: untere Schranke für f
        double mu1,     // i: Parameter mu_1 aus Alg.2.1
        double sigma,   // i: Parameter sigma aus Alg.2.1
        double xi1,     // i: Parameter xi_1 aus Alg.2.1
        double xi2,     // i: Parameter xi_1 aus Alg.2.1
        double tau,     // i: Parameter tau aus Alg.2.1
        double tau1,    // i: Parameter tau_1 aus Alg.2.1
        double tau2);   // i: Parameter tau_2 aus Alg.2.1
    // o:  0 bei erfolgreicher Liniensuche
    //    -1 bei Abbruch wegen Unterschreiten von fmin
    //     1 bei Abbruch, aus sonstigen Grnden

}  // namespace meshit

#endif

