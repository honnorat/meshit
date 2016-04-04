#include "geomtest3d.hpp"

#include "../linalg/densemat.hpp"

namespace meshit
{
    int IntersectTriangleLine(const Point3d** tri, const Point3d** line)
    {
        Vec3d vl(*line[0], *line[1]);
        Vec3d vt1(*tri[0], *tri[1]);
        Vec3d vt2(*tri[0], *tri[2]);
        Vec3d vrs(*tri[0], *line[0]);

        // static DenseMatrix a(3), ainv(3);
        // static Vector rs(3), lami(3);
        Mat<3, 3> a, ainv;
        Vec3d rs, lami;
        int i;

        for (i = 0; i < 3; i++) {
            a(i, 0) = -vl.X(i + 1);
            a(i, 1) = vt1.X(i + 1);
            a(i, 2) = vt2.X(i + 1);
            rs.X(i + 1) = vrs.X(i + 1);
        }

        double det = Det(a);
        double arel = vl.Length() * vt1.Length() * vt2.Length();

        // new !!!!
        if (fabs(det) <= 1e-10 * arel) {
            return 0;
        }

        CalcInverse(a, ainv);
        lami = ainv * rs;

        if (lami.X() >= 0 && lami.X() <= 1 &&
            lami.Y() >= 0 && lami.Z() >= 0 && lami.Y() + lami.Z() <= 1) {
            return 1;
        }

        return 0;
    }

    int IntersectTriangleTriangle(const Point3d** tri1, const Point3d** tri2)
    {
        int i, j;
        double diam = Dist(*tri1[0], *tri1[1]);
        double epsrel = 1e-8;
        double eps = diam * epsrel;
        double eps2 = eps * eps;


        int cnt = 0;
        for (i = 0; i <= 2; i++) {
            for (j = 0; j <= 2; j++) {
                if (Dist2(*tri1[j], *tri2[i]) < eps2) {
                    cnt++;
                    break;
                }
            }
        }

        switch (cnt) {
            case 0: {
                const Point3d* line[2];

                for (i = 0; i <= 2; i++) {
                    line[0] = tri2[i];
                    line[1] = tri2[(i + 1) % 3];

                    if (IntersectTriangleLine(tri1, &line[0])) {
                        std::cerr << "int1, line = " << *line[0] << " - " << *line[1] << std::endl;
                        return 1;
                    }
                }

                for (i = 0; i <= 2; i++) {
                    line[0] = tri1[i];
                    line[1] = tri1[(i + 1) % 3];

                    if (IntersectTriangleLine(tri2, &line[0])) {
                        std::cerr << "int2, line = " << *line[0] << " - " << *line[1] << std::endl;
                        return 1;
                    }
                }
                break;
            }
            default:
                return 0;
        }

        return 0;
    }

    void LocalCoordinates(const Vec3d& e1, const Vec3d& e2,
                          const Vec3d& v, double& lam1, double& lam2)
    {
        double m11 = e1 * e1;
        double m12 = e1 * e2;
        double m22 = e2 * e2;
        double rs1 = v * e1;
        double rs2 = v * e2;

        double det = m11 * m22 - m12 * m12;
        lam1 = (rs1 * m22 - rs2 * m12) / det;
        lam2 = (m11 * rs2 - m12 * rs1) / det;
    }

    int CalcSphereCenter(const Point3d** pts, Point3d& c)
    {
        Vec3d row1(*pts[0], *pts[1]);
        Vec3d row2(*pts[0], *pts[2]);
        Vec3d row3(*pts[0], *pts[3]);

        Vec3d rhs(0.5 * (row1 * row1),
                  0.5 * (row2 * row2),
                  0.5 * (row3 * row3));
        Transpose(row1, row2, row3);

        Vec3d sol;
        if (SolveLinearSystem(row1, row2, row3, rhs, sol)) {
            std::cerr << "CalcSphereCenter: degenerated" << std::endl;
            return 1;
        }

        c = *pts[0] + sol;
        return 0;
    }

    double ComputeCylinderRadius(
        const Point3d& p1,
        const Point3d& p2,
        const Point3d& p3,
        const Point3d& p4)
    {
        Vec3d v12(p1, p2);
        Vec3d v13(p1, p3);
        Vec3d v14(p1, p4);

        Vec3d n1 = Cross(v12, v13);
        Vec3d n2 = Cross(v14, v12);

        double n1l = n1.Length();
        double n2l = n2.Length();
        n1 /= n1l;
        n2 /= n2l;

        double v12len = v12.Length();
        double h1 = n1l / v12len;
        double h2 = n2l / v12len;

        /*
        std::cerr << "n1 = " << n1 << " n2 = " << n2 
               << "h1 = " << h1 << " h2 = " << h2 <<std::endl;
         */
        return ComputeCylinderRadius(n1, n2, h1, h2);
    }

    /*
      Two triangles T1 and T2 have normals n1 and n2.
      The height over the common edge is h1, and h2.
     */
    double ComputeCylinderRadius(const Vec3d& n1, const Vec3d& n2,
                                 double h1, double h2)
    {
        Vec3d t1, t2;
        double n11 = n1 * n1;
        double n12 = n1 * n2;
        double n22 = n2 * n2;
        double det = n11 * n22 - n12 * n12;

        if (fabs(det) < 1e-14 * n11 * n22)
            return 1e20;

        // a biorthogonal bases   (ti * nj) = delta_ij:
        t1 = (n22 / det) * n1 + (-n12 / det) * n2;
        t2 = (-n12 / det) * n1 + (n11 / det) * n2;

        // normalize:
        t1 /= t1.Length();
        t2 /= t2.Length();

        /*
          vector to center point has form
          v = lam1 n1 + lam2 n2
          and fulfills
          t2 v = h1/2
          t1 v = h2/2
         */

        double lam1 = 0.5 * h2 / (n1 * t1);
        double lam2 = 0.5 * h1 / (n2 * t2);

        double rad = (lam1 * n1 + lam2 * n2).Length();
        /*
        std::cerr << "n1 = " << n1
               << " n2 = " << n2
               << " t1 = " << t1
               << " t2 = " << t2
               << " rad = " << rad <<std::endl;
         */
        return rad;
    }

    double MinDistLP2(const Point2d& lp1, const Point2d& lp2, const Point2d& p)
    {
        Vec2d v(lp1, lp2);
        Vec2d vlp(lp1, p);

        // dist(lam) = \| vlp \|^2 - 2 lam (v1p, v) + lam^2 \| v \|^2

        // lam = (v * vlp) / (v * v);
        // if (lam < 0) lam = 0;
        // if (lam > 1) lam = 1;

        double num = v * vlp;
        double den = v * v;

        if (num <= 0)
            return Dist2(lp1, p);

        if (num >= den)
            return Dist2(lp2, p);

        if (den > 0) {
            return vlp.Length2() - num * num / den;
        }
        else
            return vlp.Length2();
    }

    double MinDistLP2(const Point3d& lp1, const Point3d& lp2, const Point3d& p)
    {
        Vec3d v(lp1, lp2);
        Vec3d vlp(lp1, p);

        // dist(lam) = \| vlp \|^2 - 2 lam (v1p, v) + lam^2 \| v \|^2

        // lam = (v * vlp) / (v * v);
        // if (lam < 0) lam = 0;
        // if (lam > 1) lam = 1;

        double num = v * vlp;
        double den = v * v;

        if (num <= 0)
            return Dist2(lp1, p);

        if (num >= den)
            return Dist2(lp2, p);

        if (den > 0) {
            return vlp.Length2() - num * num / den;
        }
        else
            return vlp.Length2();
    }

    double MinDistTP2(const Point3d& tp1, const Point3d& tp2,
                      const Point3d& tp3, const Point3d& p)
    {
        double lam1, lam2;
        double res;

        LocalCoordinates(Vec3d(tp1, tp2), Vec3d(tp1, tp3),
                         Vec3d(tp1, p), lam1, lam2);
        int in1 = lam1 >= 0;
        int in2 = lam2 >= 0;
        int in3 = lam1 + lam2 <= 1;

        if (in1 && in2 && in3) {
            Point3d pp = tp1 + lam1 * Vec3d(tp1, tp2) + lam2 * Vec3d(tp1, tp3);
            res = Dist2(p, pp);
        } else {
            res = Dist2(tp1, p);
            if (!in1) {
                double hv = MinDistLP2(tp1, tp3, p);
                if (hv < res) res = hv;
            }
            if (!in2) {
                double hv = MinDistLP2(tp1, tp2, p);
                if (hv < res) res = hv;
            }
            if (!in3) {
                double hv = MinDistLP2(tp2, tp3, p);
                if (hv < res) res = hv;
            }
        }
        return res;
    }

}  // namespace meshit
