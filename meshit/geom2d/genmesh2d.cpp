#include "../meshit.hpp"
#include "geometry2d.hpp"
#include "../meshing/meshing2.hpp"
#include "../meshing/global.hpp"

namespace meshit {

    void Optimize2d(Mesh & mesh, MeshingParameters & mp)
    {
        mesh.CalcSurfacesOfNode();

        const char * optstr = mp.optimize2d;
        int optsteps = mp.optsteps2d;

        for (int i = 1; i <= optsteps; i++)
            for (size_t j = 1; j <= strlen(optstr); j++) {

                switch (optstr[j - 1]) {
                    case 's':
                    { // topological swap
                        MeshOptimize2d meshopt;
                        meshopt.SetMetricWeight(0);
                        meshopt.EdgeSwapping(mesh, 0);
                        break;
                    }
                    case 'S':
                    { // metric swap
                        MeshOptimize2d meshopt;
                        meshopt.SetMetricWeight(0);
                        meshopt.EdgeSwapping(mesh, 1);
                        break;
                    }
                    case 'm':
                    {
                        MeshOptimize2d meshopt;
                        meshopt.SetMetricWeight(1);
                        meshopt.ImproveMesh(mesh, mp);
                        break;
                    }
                    case 'c':
                    {
                        MeshOptimize2d meshopt;
                        meshopt.SetMetricWeight(0.2);
                        meshopt.CombineImprove(mesh);
                        break;
                    }
                    default:
                        MESHIT_LOG_ERROR("Optimization code " << optstr[j - 1] << " not defined");
                }
            }
    }

    void CalcPartition(
            const SplineSegExt & spline,
            // double l, 
            MeshingParameters & mp, Mesh & mesh,
            // double h, double h1, double h2, double hcurve, 
            double elto0, Array<double> & points)
    {
        double fperel, oldf, f;

        int n = 1000;

        Array<Point<2> > xi(n);
        Array<double> hi(n);

        for (int i = 0; i < n; i++) {
            xi[i] = spline.GetPoint((i + 0.5) / n);
            hi[i] = mesh.GetH(Point3d(xi[i](0), xi[i](1), 0));
        }

        // limit slope
        double gradh = 1 / elto0;
        for (int i = 0; i < n - 1; i++) {
            double hnext = hi[i] + gradh * (xi[i + 1] - xi[i]).Length();
            hi[i + 1] = std::min(hi[i + 1], hnext);
        }
        for (int i = n - 1; i > 1; i--) {
            double hnext = hi[i] + gradh * (xi[i - 1] - xi[i]).Length();
            hi[i - 1] = std::min(hi[i - 1], hnext);
        }

        points.resize(0);

        double len = spline.Length();
        double dt = len / n;

        double sum = 0;
        for (int i = 1; i <= n; i++) {
            double fun = hi[i - 1];
            sum += dt / fun;
        }

        int nel = int (sum + 1);
        fperel = sum / nel;

        points.push_back(0);

        int i = 1;
        oldf = 0;

        for (int j = 1; j <= n && i < nel; j++) {
            double t = (j - 0.5) * dt;
            double fun = hi[j - 1];

            f = oldf + dt / fun;

            while (i * fperel < f && i < nel) {
                points.push_back(dt * (j - 1) + (i * fperel - oldf) * fun);
                i++;
            }
            oldf = f;
            t += dt;
        }
        points.push_back(len);
    }


    // partitionizes spline curve

    void Partition(const SplineSegExt & spline,
            MeshingParameters & mp, double hxxx, double elto0,
            Mesh & mesh, Point3dTree & searchtree, int segnr)
    {
        int n = 100;

        Point<2> mark, oldmark;
        Array<double> curvepoints;
        double edgelength, edgelengthold;

        CalcPartition(spline, mp, mesh, elto0, curvepoints);

        double dt = 1.0 / n;

        int j = 1;

        Point<2> pold = spline.GetPoint(0);
        double lold = 0;
        oldmark = pold;
        edgelengthold = 0;
        Array<int> locsearch;

        for (int i = 1; i <= n; i++) {
            Point<2> p = spline.GetPoint(i * dt);
            double l = lold + Dist(p, pold);
            while (j < curvepoints.size() && (l >= curvepoints[j] || i == n)) {
                double frac = (curvepoints[j] - lold) / (l - lold);
                edgelength = i * dt + (frac - 1) * dt;
                mark = spline.GetPoint(edgelength);

                {
                    PointIndex pi1 = -1, pi2 = -1;

                    Point3d mark3(mark(0), mark(1), 0);
                    Point3d oldmark3(oldmark(0), oldmark(1), 0);

                    double h = mesh.GetH(Point<3> (oldmark(0), oldmark(1), 0));
                    Vec<3> v(1e-4 * h, 1e-4 * h, 1e-4 * h);
                    searchtree.GetIntersecting(oldmark3 - v, oldmark3 + v, locsearch);

                    for (int k = 0; k < locsearch.size(); k++) {
                        if (mesh[PointIndex(locsearch[k])].GetLayer() == spline.layer) {
                            pi1 = locsearch[k];
                        }
                    }

                    searchtree.GetIntersecting(mark3 - v, mark3 + v, locsearch);
                    for (int k = 0; k < locsearch.size(); k++) {
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
                    seg.si = spline.bc; // segnr;
                    seg[0] = pi1;
                    seg[1] = pi2;
                    seg.domin = spline.leftdom;
                    seg.domout = spline.rightdom;
                    seg.epgeominfo[0].edgenr = segnr;
                    seg.epgeominfo[0].dist = edgelengthold;
                    seg.epgeominfo[1].edgenr = segnr;
                    seg.epgeominfo[1].dist = edgelength;
                    seg.singedge_left = spline.hpref_left;
                    seg.singedge_right = spline.hpref_right;
                    mesh.AddSegment(seg);
                }

                oldmark = mark;
                edgelengthold = edgelength;
                j++;
            }

            pold = p;
            lold = l;
        }
    }

    void SplineGeometry2d::PartitionBoundary(MeshingParameters & mp, double h, Mesh & mesh2d)
    {
        Box<2> bbox;
        GetBoundingBox(bbox);
        double dist = Dist(bbox.PMin(), bbox.PMax());
        Point<3> pmin;
        Point<3> pmax;

        pmin(2) = -dist;
        pmax(2) = dist;
        for (int j = 0; j < 2; j++) {
            pmin(j) = bbox.PMin()(j);
            pmax(j) = bbox.PMax()(j);
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
            const SplineSegExt & spline = GetSpline(i);
            const GeomPoint<2> & p1 = spline.StartPI();
            const GeomPoint<2> & p2 = spline.EndPI();

            double h1 = std::min(p1.hmax, h / p1.refatpoint);
            double h2 = std::min(p2.hmax, h / p2.refatpoint);
            mesh2d.RestrictLocalH(Point3d(p1(0), p1(1), 0), h1);
            mesh2d.RestrictLocalH(Point3d(p2(0), p2(1), 0), h2);

            double len = spline.Length();
            mesh2d.RestrictLocalHLine(
                    Point3d(p1(0), p1(1), 0),
                    Point3d(p2(0), p2(1), 0), len / mp.segmentsperedge);

            double hcurve = std::min(spline.hmax, h / spline.reffak);
            double hl = GetDomainMaxh(spline.leftdom);
            double hr = GetDomainMaxh(spline.rightdom);
            if (hl > 0) hcurve = std::min(hcurve, hl);
            if (hr > 0) hcurve = std::min(hcurve, hr);

            int np = 1000;
            for (double t = 0.5 / np; t < 1; t += 1.0 / np) {
                Point<2> x = spline.GetPoint(t);
                double hc = 1.0 / mp.curvaturesafety / (1e-99 + spline.CalcCurvature(t));
                mesh2d.RestrictLocalH(Point3d(x(0), x(1), 0), std::min(hc, hcurve));
            }
        }

        for (size_t i = 0; i < splines.size(); i++)
            if (GetSpline(i).copyfrom == -1) {
                // astrid - set boundary meshsize to  domain meshsize h
                // if no domain mesh size is given, the max h value from the bounding box is used
                double hl = GetDomainMaxh(GetSpline(i).leftdom);
                double hr = GetDomainMaxh(GetSpline(i).rightdom);

                double useh = h;
                if (hl > 0) useh = std::min(h, hl);
                if (hr > 0) useh = std::min(h, hr);
                Partition(GetSpline(i), mp, useh, elto0, mesh2d, searchtree, i + 1);
            }
            else {
                CopyEdgeMesh(GetSpline(i).copyfrom, i + 1, mesh2d, searchtree);
            }
    }

    void SplineGeometry2d::CopyEdgeMesh(int from, int to, Mesh & mesh, Point3dTree & searchtree)
    {
        Array<int, PointIndex::BASE> mappoints(mesh.GetNP());
        Array<double, PointIndex::BASE> param(mesh.GetNP());
        mappoints = -1;
        param = 0;

        Point3d pmin, pmax;
        mesh.GetBox(pmin, pmax);
        double diam2 = Dist2(pmin, pmax);

        MESHIT_LOG_DEBUG("copy edge, from = " << from << " to " << to);

        for (int i = 1; i <= mesh.GetNSeg(); i++) {
            const Segment & seg = mesh.LineSegment(i);
            if (seg.edgenr == from) {
                mappoints.Elem(seg[0]) = 1;
                param.Elem(seg[0]) = seg.epgeominfo[0].dist;

                mappoints.Elem(seg[1]) = 1;
                param.Elem(seg[1]) = seg.epgeominfo[1].dist;
            }
        }

        bool mapped = false;
        for (int i = 1; i <= mappoints.size(); i++) {
            if (mappoints.Get(i) != -1) {
                Point<2> newp = splines[to+1]->GetPoint(param.Get(i));
                Point<3> newp3(newp(0), newp(1), 0);

                int npi = -1;

                for (PointIndex pi = PointIndex::BASE;
                        pi < mesh.GetNP() + PointIndex::BASE; pi++)
                    if (Dist2(mesh.Point(pi), newp3) < 1e-12 * diam2)
                        npi = pi;

                if (npi == -1) {
                    npi = mesh.AddPoint(newp3);
                    searchtree.Insert(newp3, npi);
                }

                mappoints.Elem(i) = npi;

                mesh.GetIdentifications().Add(i, npi, to);
                mapped = true;
            }
        }
        if (mapped)
            mesh.GetIdentifications().SetType(to, Identifications::PERIODIC);

        // copy segments
        int oldnseg = mesh.GetNSeg();
        for (int i = 1; i <= oldnseg; i++) {
            const Segment & seg = mesh.LineSegment(i);
            if (seg.edgenr == from) {
                Segment nseg;
                nseg.edgenr = to;
                nseg.si = GetSpline(to - 1).bc; // splines.Get(to)->bc;
                nseg[0] = mappoints.Get(seg[0]);
                nseg[1] = mappoints.Get(seg[1]);
                nseg.domin = GetSpline(to - 1).leftdom;
                nseg.domout = GetSpline(to - 1).rightdom;

                nseg.epgeominfo[0].edgenr = to;
                nseg.epgeominfo[0].dist = param.Get(seg[0]);
                nseg.epgeominfo[1].edgenr = to;
                nseg.epgeominfo[1].dist = param.Get(seg[1]);
                mesh.AddSegment(nseg);
            }
        }
    }
}
