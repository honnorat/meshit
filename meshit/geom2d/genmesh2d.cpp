#include <sstream>

#include "geometry2d.hpp"
#include "../meshing/meshing2.hpp"
#include "../meshing/global.hpp"

namespace meshit
{
    void Optimize2d(Mesh& mesh, MeshingParameters& mp)
    {
        mesh.CalcSurfacesOfNode();

        const char* optstr = mp.optimize2d;
        int optsteps = mp.optsteps2d;

        for (int i = 1; i <= optsteps; i++) {
            for (size_t j = 0; j < strlen(optstr); j++) {
                switch (optstr[j]) {
                    case 's': {  // topological swap
                        MeshOptimize2d meshopt;
                        meshopt.SetMetricWeight(0);
                        meshopt.EdgeSwapping(mesh, 0);
                        mp.n_steps++;
                        break;
                    }
                    case 'S': {  // metric swap
                        MeshOptimize2d meshopt;
                        meshopt.SetMetricWeight(0);
                        meshopt.EdgeSwapping(mesh, 1);
                        mp.n_steps++;
                        break;
                    }
                    case 'm': {
                        MeshOptimize2d meshopt;
                        meshopt.SetMetricWeight(1);
                        meshopt.ImproveMesh(mesh, mp);
                        mp.n_steps++;
                        break;
                    }
                    case 'c': {
                        MeshOptimize2d meshopt;
                        meshopt.SetMetricWeight(0.2);
                        meshopt.CombineImprove(mesh);
                        mp.n_steps++;
                        break;
                    }
                    case 'p': {
                        // print mesh
                        std::stringstream mesh_name;
                        mesh_name << "mesh_debug_" << std::setfill('0') << std::setw(2) << mp.n_steps << ".msh";
                        mesh.Export(mesh_name.str());
                        break;
                    }
                    default:
                        MESHIT_LOG_ERROR("Optimization code " << optstr[j] << " not defined");
                }
            }
        }
    }

    void CalcPartition(const SplineSegExt& spline,
                       MeshingParameters& mp,
                       Mesh& mesh,
                       double elto0,
                       std::vector<double>& points)
    {
        double fperel, oldf, f;

        size_t n = 1000;

        std::vector<Point<2> > xi(n);
        std::vector<double> hi(n);

        for (size_t i = 0; i < n; i++) {
            xi[i] = spline.GetPoint((i + 0.5) / n);
            hi[i] = mesh.GetH(Point3d(xi[i][0], xi[i][1], 0));
        }

        // limit slope
        double gradh = 1 / elto0;
        for (size_t i = 0; i < n - 1; i++) {
            double hnext = hi[i] + gradh * (xi[i + 1] - xi[i]).Length();
            hi[i + 1] = std::min(hi[i + 1], hnext);
        }
        for (size_t i = n - 1; i > 1; i--) {
            double hnext = hi[i] + gradh * (xi[i - 1] - xi[i]).Length();
            hi[i - 1] = std::min(hi[i - 1], hnext);
        }

        points.resize(0);

        double len = spline.Length();
        double dt = len / n;

        double sum = 0;
        for (size_t i = 0; i < n; i++) {
            sum += dt / hi[i];
        }

        size_t nel = static_cast<size_t>(sum + 1);
        fperel = sum / nel;

        points.push_back(0);

        size_t i = 1;
        oldf = 0;

        for (size_t j = 0; j < n && i < nel; j++) {
            double fun = hi[j];

            f = oldf + dt / fun;

            while (i * fperel < f && i < nel) {
                points.push_back(dt * j + (i * fperel - oldf) * fun);
                i++;
            }
            oldf = f;
        }
        points.push_back(len);
    }

    // partitionizes spline curve

    void Partition(const SplineSegExt& spline,
                   MeshingParameters& mp,
                   double elto0,
                   Mesh& mesh,
                   Point3dTree& searchtree,
                   int segnr)
    {
        size_t n = 100;

        Point<2> mark, oldmark;
        std::vector<double> curvepoints;
        double edgelength, edgelengthold;

        CalcPartition(spline, mp, mesh, elto0, curvepoints);

        double dt = 1.0 / n;

        size_t j = 1;

        Point<2> pold = spline.GetPoint(0);
        double lold = 0.0;
        oldmark = pold;
        edgelengthold = 0;
        std::vector<size_t> locsearch;

        for (size_t i = 1; i <= n; i++) {
            double t = static_cast<double>(i) * dt;
            Point<2> p = spline.GetPoint(t);
            double l = lold + Dist(p, pold);
            while (j < curvepoints.size() && (l >= curvepoints[j] || i == n)) {
                double frac = (curvepoints[j] - l) / (l - lold);
                edgelength = t + frac * dt;
                mark = spline.GetPoint(edgelength);

                PointIndex pi1{-1}, pi2{-1};
                Point3d mark3(mark[0], mark[1], 0);
                Point3d oldmark3(oldmark[0], oldmark[1], 0);

                double h = mesh.GetH(Point3d(oldmark[0], oldmark[1], 0));
                Vec3d v(1e-4 * h, 1e-4 * h, 1e-4 * h);
                searchtree.GetIntersecting(oldmark3 - v, oldmark3 + v, locsearch);

                for (size_t k = 0; k < locsearch.size(); k++) {
                    if (mesh[PointIndex(locsearch[k])].GetLayer() == spline.layer) {
                        pi1 = locsearch[k];
                    }
                }
                searchtree.GetIntersecting(mark3 - v, mark3 + v, locsearch);
                for (size_t k = 0; k < locsearch.size(); k++) {
                    if (mesh[PointIndex(locsearch[k])].GetLayer() == spline.layer) {
                        pi2 = locsearch[k];
                    }
                }
                if (pi1 == -1) {
                    pi1 = mesh.AddPoint(oldmark3, spline.layer);
                    searchtree.Insert(oldmark3, pi1);
                }
                if (pi2 == -1) {
                    pi2 = mesh.AddPoint(mark3, spline.layer);
                    searchtree.Insert(mark3, pi2);
                }

                Segment seg;
                seg.edgenr = segnr;
                seg.si = spline.bc;  // segnr;
                seg[0] = pi1;
                seg[1] = pi2;
                seg.domin = spline.leftdom;
                seg.domout = spline.rightdom;
                seg.epgeominfo[0].edgenr = segnr;
                seg.epgeominfo[0].dist = edgelengthold;
                seg.epgeominfo[1].edgenr = segnr;
                seg.epgeominfo[1].dist = edgelength;
                mesh.AddSegment(seg);

                oldmark = mark;
                edgelengthold = edgelength;
                j++;
            }

            pold = p;
            lold = l;
        }
    }

    void SplineGeometry2d::PartitionBoundary(MeshingParameters& mp, double h, Mesh& mesh2d)
    {
        Box<2> bbox;
        GetBoundingBox(bbox);
        double dist = Dist(bbox.PMin(), bbox.PMax());
        Point3d pmin;
        Point3d pmax;

        pmin[2] = -dist;
        pmax[2] = dist;
        for (int j = 0; j < 2; j++) {
            pmin[j] = bbox.PMin()[j];
            pmax[j] = bbox.PMax()[j];
        }

        Point3dTree searchtree(pmin, pmax);

        for (size_t i = 0; i < splines.size(); i++) {
            for (int side = 0; side <= 1; side++) {
                int dom = (side == 0) ? GetSpline(i).leftdom : GetSpline(i).rightdom;
                if (dom != 0) GetSpline(i).layer = GetDomainLayer(dom);
            }
        }

        // mesh size restrictions ...

        for (size_t i = 0; i < splines.size(); i++) {
            const SplineSegExt& spline = GetSpline(i);
            const GeomPoint<2>& p1 = spline.StartPI();
            const GeomPoint<2>& p2 = spline.EndPI();

            double h1 = std::min(p1.hmax, h / p1.refatpoint);
            double h2 = std::min(p2.hmax, h / p2.refatpoint);
            mesh2d.RestrictLocalH(Point3d(p1[0], p1[1], 0), h1);
            mesh2d.RestrictLocalH(Point3d(p2[0], p2[1], 0), h2);

            double len = spline.Length();
            mesh2d.RestrictLocalHLine(Point3d(p1[0], p1[1], 0),
                                      Point3d(p2[0], p2[1], 0), len / mp.segments_per_edge);

            double hcurve = std::min(spline.hmax, h / spline.reffak);
            double hl = GetDomainMaxh(spline.leftdom);
            double hr = GetDomainMaxh(spline.rightdom);
            if (hl > 0) hcurve = std::min(hcurve, hl);
            if (hr > 0) hcurve = std::min(hcurve, hr);

            int np = 1000;
            for (double t = 0.5 / np; t < 1; t += 1.0 / np) {
                Point<2> x = spline.GetPoint(t);
                double hc = 1.0 / mp.curvature_safety / (1e-99 + spline.CalcCurvature(t));
                mesh2d.RestrictLocalH(Point3d(x[0], x[1], 0), std::min(hc, hcurve));
            }
        }

        for (size_t i = 0; i < splines.size(); i++) {
            Partition(GetSpline(i), mp, elto0, mesh2d, searchtree, i + 1);
        }
    }

}  // namespace meshit
