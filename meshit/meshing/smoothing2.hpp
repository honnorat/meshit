#ifndef MESHIT_SMOOTHING_HPP
#define MESHIT_SMOOTHING_HPP
/**
 * meshit - a 2d mesh generator
 *
 * Copyright © 1995-2015 Joachim Schoeberl <joachim.schoeberl@tuwien.ac.at>
 * Copyright © 2015-2016 Marc Honnorat <marc.honnorat@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library in the file LICENSE.LGPL; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 */

#include "../linalg/densemat.hpp"
#include "mesh_types.hpp"

namespace meshit {

class Opti2dLocalData
{
 public:
    explicit Opti2dLocalData(double metric_weight)
        : loc_metric_weight{metric_weight} { }

 public:
    MeshPoint sp1;
    std::vector<ElementIndex> loc_elements;
    std::vector<double> lochs;
    std::vector<MeshPoint> loc_pnts2;
    std::vector<MeshPoint> loc_pnts3;
    double loc_metric_weight;
};

class Mesh;

/**
    Function to be minimized.
 */
class MinFunction_2d
{
 public:
    MinFunction_2d(const Mesh& amesh, Opti2dLocalData& ald)
        : mesh_{amesh}, ld(ald) { }

    double Func(const double* x);

    // function + gradient
    double FuncGrad(const double* x, double* g);

    // directional derivative
    double FuncDeriv(const double* x, const double* dir, double& deriv);

 protected:
    const Mesh& mesh_;
    Opti2dLocalData& ld;
};

struct OptiParameters
{
    OptiParameters()
        : maxit_bfgs{100}, maxit_linsearch{100}, eps{1e-2}, typf{1.0}, typx{1.0} { }

    size_t maxit_bfgs;
    size_t maxit_linsearch;
    double eps;
    double typf;
    double typx;
};

/** Implementation of BFGS method.
    Efficient method for non-linear minimiztion problems.
    @param x initial value and solution
    @param fun function to be minimized
 */
double BFGS_2d(double* x, MinFunction_2d& fun, const OptiParameters& par);

int lines(
    double* x,     // i: Ausgangspunkt der Liniensuche
    double* xneu,  // o: Loesung der Liniensuche bei Erfolg
    double* p,     // i: Suchrichtung
    double& f,     // i: Funktionswert an der Stelle x, o: Funktionswert an der Stelle xneu, falls ifail = 0
    double* g,     // i: Gradient an der Stelle x, o: Gradient an der Stelle xneu, falls ifail = 0
    MinFunction_2d& fun,        // function to minmize
    const OptiParameters& par,  // parameters
    double& alphahat,           // i: Startwert für alpha_hat, o: Loesung falls ifail = 0
    double fmin,                // i: untere Schranke für f
    double mu1,                 // i: Parameter mu_1 aus Alg.2.1
    double sigma,               // i: Parameter sigma aus Alg.2.1
    double xi1,                 // i: Parameter xi_1 aus Alg.2.1
    double xi2,                 // i: Parameter xi_1 aus Alg.2.1
    double tau,                 // i: Parameter tau aus Alg.2.1
    double tau1,                // i: Parameter tau_1 aus Alg.2.1
    double tau2);               // i: Parameter tau_2 aus Alg.2.1
// o:  0 bei erfolgreicher Liniensuche
//    -1 bei Abbruch wegen Unterschreiten von fmin
//     1 bei Abbruch, aus sonstigen Grnden

}  // namespace meshit

#endif
