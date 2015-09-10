#include <iomanip>
#include <meshit.hpp>
#include "meshtool.hpp"
#include "../gprim/geomtest3d.hpp"

namespace meshit {

    int CheckSurfaceMesh(const Mesh & mesh)
    {
        LOG_DEBUG("Check Surface mesh");

        int nf = mesh.GetNSE();
        INDEX_2_HASHTABLE<int> edges(nf + 2);
        INDEX_2 i2;
        int cnt1 = 0, cnt2 = 0;

        for (int i = 1; i <= nf; i++) {
            for (int j = 1; j <= 3; j++) {
                i2.I1() = mesh.SurfaceElement(i).PNumMod(j);
                i2.I2() = mesh.SurfaceElement(i).PNumMod(j + 1);
                if (edges.Used(i2)) {
                    int hi = edges.Get(i2);
                    if (hi != 1)
                        LOG_ERROR("CheckSurfaceMesh, hi = " << hi);
                    edges.Set(i2, 2);
                    cnt2++;
                }
                else {
                    std::swap(i2.I1(), i2.I2());
                    edges.Set(i2, 1);
                    cnt1++;
                }
            }
        }

        if (cnt1 != cnt2) {
            LOG_ERROR("Surface mesh not consistent : cnt1 = " << cnt1 << " / cnt2 = " << cnt2);
            return 0;
        }
        return 1;
    }

    int CheckSurfaceMesh2(const Mesh & mesh)
    {
        const Point<3> *tri1[3], *tri2[3];

        for (int i = 1; i <= mesh.GetNOpenElements(); i++) {
            for (int j = 1; j < i; j++) {
                for (int k = 1; k <= 3; k++) {
                    tri1[k - 1] = &mesh.Point(mesh.OpenElement(i).PNum(k));
                    tri2[k - 1] = &mesh.Point(mesh.OpenElement(j).PNum(k));
                }
                if (IntersectTriangleTriangle(&tri1[0], &tri2[0])) {
                    LOG_ERROR("Surface elements are intersecting :");
                    for (int k = 0; k <= 2; k++) {
                        LOG_ERROR_CONT(*tri1[k] << "  ");
                    }
                    std::cerr << std::endl;
                    for (int k = 0; k <= 2; k++)
                        std::cerr << *tri2[k] << "   ";
                    std::cerr << std::endl;
                }

            }
        }
        return 0;
    }

    static double TriangleQualityInst(const Point3d & p1, const Point3d & p2,
            const Point3d & p3)
    {
        // quality 0 (worst) .. 1 (optimal)

        Vec3d v1, v2, v3;
        double s1, s2, s3;
        double an1, an2, an3;

        v1 = p2 - p1;
        v2 = p3 - p1;
        v3 = p3 - p2;

        an1 = Angle(v1, v2);
        v1 *= -1;
        an2 = Angle(v1, v3);
        an3 = Angle(v2, v3);

        s1 = sin(an1 / 2);
        s2 = sin(an2 / 2);
        s3 = sin(an3 / 2);

        return 8 * s1 * s2 * s3;
    }

    void MeshQuality2d(const Mesh & mesh)
    {
        int ncl = 20, cl;
        Array<INDEX> incl(ncl);
        INDEX i;
        double qual;

        incl = 0;

        for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++) {
            qual = TriangleQualityInst(
                    mesh[mesh.SurfaceElement(sei)[0]],
                    mesh[mesh.SurfaceElement(sei)[1]],
                    mesh[mesh.SurfaceElement(sei)[2]]);

            cl = int ( (ncl - 1e-3) * qual) + 1;
            incl.Elem(cl)++;
        }

        LOG_INFO("\n\n");
        LOG_INFO("Points:           " << mesh.GetNP());
        LOG_INFO("Surface Elements: " << mesh.GetNSE());
        LOG_INFO("\nElements in qualityclasses:");
        for (i = 1; i <= ncl; i++) {
            LOG_INFO(std::fixed << std::setprecision(2) <<
                    std::setw(4) << double (i - 1) / ncl << " - " <<
                    std::setw(4) << double (i) / ncl << ": " << incl.Get(i));
        }
    }

    void SaveEdges(const Mesh & mesh, const char * geomfile, double h, char * filename)
    {
        std::ofstream of(filename);
        int i;
        const Segment * seg;

        of << "edges" << std::endl;
        of << geomfile << std::endl;
        of << h << std::endl;

        of << mesh.GetNP() << std::endl;
        for (i = 1; i <= mesh.GetNP(); i++)
            of << mesh.Point(i)(0) << " "
            << mesh.Point(i)(1) << " "
            << mesh.Point(i)(2) << "\n";

        of << 2 * mesh.GetNSeg() << std::endl;
        for (i = 1; i <= mesh.GetNSeg(); i++) {

            seg = &mesh.LineSegment(i);

            of << (*seg)[1] << " " << (*seg)[0] << " " << seg->si << "\n";
        }

    }

    void SaveSurfaceMesh(const Mesh & mesh,
            double h,
            char * filename)

    {
        INDEX i;

        std::ofstream outfile(filename);

        outfile << "surfacemesh" << std::endl;
        outfile << h << std::endl;

        outfile << mesh.GetNP() << std::endl;
        for (i = 1; i <= mesh.GetNP(); i++)
            outfile << mesh.Point(i)(0) << " "
            << mesh.Point(i)(1) << " "
            << mesh.Point(i)(2) << std::endl;



        outfile << mesh.GetNSE() << std::endl;
        for (i = 1; i <= mesh.GetNSE(); i++) {
            const Element2d & el = mesh.SurfaceElement(i);

            if (mesh.GetFaceDescriptor(el.GetIndex()).DomainOut() == 0)
                outfile << mesh.SurfaceElement(i).PNum(1) << " "
                << mesh.SurfaceElement(i).PNum(2) << " "
                << mesh.SurfaceElement(i).PNum(3) << std::endl;

            if (mesh.GetFaceDescriptor(el.GetIndex()).DomainIn() == 0)
                outfile << mesh.SurfaceElement(i).PNum(1) << " "
                << mesh.SurfaceElement(i).PNum(3) << " "
                << mesh.SurfaceElement(i).PNum(2) << std::endl;
        }
    }
}
