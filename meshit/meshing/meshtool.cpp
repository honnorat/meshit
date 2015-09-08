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
        int i, j, k;
        const Point<3> *tri1[3], *tri2[3];

        for (i = 1; i <= mesh.GetNOpenElements(); i++) {
            PrintDot();
            for (j = 1; j < i; j++) {
                for (k = 1; k <= 3; k++) {
                    tri1[k - 1] = &mesh.Point(mesh.OpenElement(i).PNum(k));
                    tri2[k - 1] = &mesh.Point(mesh.OpenElement(j).PNum(k));
                }
                if (IntersectTriangleTriangle(&tri1[0], &tri2[0])) {
                    PrintSysError("Surface elements are intersecting");
                    std::cerr << "Intersecting: " << std::endl;
                    for (k = 0; k <= 2; k++)
                        std::cerr << *tri1[k] << "   ";
                    std::cerr << std::endl;
                    for (k = 0; k <= 2; k++)
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

    static double TetElementQuality(
            const Point3d & p1, const Point3d & p2,
            const Point3d & p3, const Point3d & p4)
    {
        double vol, l, l4, l5, l6;

        Vec3d v1 = p2 - p1;
        Vec3d v2 = p3 - p1;
        Vec3d v3 = p4 - p1;

        vol = fabs((Cross(v1, v2) * v3)) / 6;
        l4 = Dist(p2, p3);
        l5 = Dist(p2, p4);
        l6 = Dist(p3, p4);

        l = v1.Length() + v2.Length() + v3.Length() + l4 + l5 + l6;

        if (vol <= 1e-8 * l * l * l) return 1e-10;

        return vol / (l * l * l) * 1832.82; // 6^4 * sqrt(2)
    }

    double CalcTetBadness(
            const Point3d & p1, const Point3d & p2,
            const Point3d & p3, const Point3d & p4, double h,
            const MeshingParameters & mp)
    {
        double vol, l, ll, lll, ll1, ll2, ll3, ll4, ll5, ll6;
        double err;

        Vec3d v1(p1, p2);
        Vec3d v2(p1, p3);
        Vec3d v3(p1, p4);

        vol = Determinant(v1, v2, v3) * (-0.166666666666666);

        ll1 = v1.Length2();
        ll2 = v2.Length2();
        ll3 = v3.Length2();
        ll4 = Dist2(p2, p3);
        ll5 = Dist2(p2, p4);
        ll6 = Dist2(p3, p4);

        ll = ll1 + ll2 + ll3 + ll4 + ll5 + ll6;
        l = sqrt(ll);
        lll = l * ll;

        if (vol <= 1e-24 * lll)
            return 1e24;

        err = 0.0080187537 * lll / vol; // sqrt(216) / (6^4 * sqrt(2))

        if (h > 0)
            err += ll / (h * h) +
            h * h * (1 / ll1 + 1 / ll2 + 1 / ll3 +
                1 / ll4 + 1 / ll5 + 1 / ll6) - 12;

        double teterrpow = mp.opterrpow;
        if (teterrpow < 1) teterrpow = 1;

        if (teterrpow == 1) return err;
        if (teterrpow == 2) return err * err;

        return pow(err, teterrpow);
    }

    double CalcTetBadnessGrad(const Point3d & p1, const Point3d & p2,
            const Point3d & p3, const Point3d & p4, double h,
            int pi, Vec<3> & grad,
            const MeshingParameters & mp)
    {
        double vol, l, ll, lll;
        double err;

        const Point3d *pp1, *pp2, *pp3, *pp4;

        pp1 = &p1;
        pp2 = &p2;
        pp3 = &p3;
        pp4 = &p4;

        switch (pi) {
            case 2:
            {
                std::swap(pp1, pp2);
                std::swap(pp3, pp4);
                break;
            }
            case 3:
            {
                std::swap(pp1, pp3);
                std::swap(pp2, pp4);
                break;
            }
            case 4:
            {
                std::swap(pp1, pp4);
                std::swap(pp3, pp2);
                break;
            }
        }


        Vec3d v1(*pp1, *pp2);
        Vec3d v2(*pp1, *pp3);
        Vec3d v3(*pp1, *pp4);

        Vec3d v4(*pp2, *pp3);
        Vec3d v5(*pp2, *pp4);
        Vec3d v6(*pp3, *pp4);

        vol = Determinant(v1, v2, v3) * (-0.166666666666666);

        Vec3d gradvol;
        Cross(v5, v4, gradvol);
        gradvol *= (-1.0 / 6.0);


        double ll1 = v1.Length2();
        double ll2 = v2.Length2();
        double ll3 = v3.Length2();
        double ll4 = v4.Length2();
        double ll5 = v5.Length2();
        double ll6 = v6.Length2();

        ll = ll1 + ll2 + ll3 + ll4 + ll5 + ll6;
        l = sqrt(ll);
        lll = l * ll;

        if (vol <= 1e-24 * lll) {
            grad = Vec3d(0, 0, 0);
            return 1e24;
        }



        Vec3d gradll1(*pp2, *pp1);
        Vec3d gradll2(*pp3, *pp1);
        Vec3d gradll3(*pp4, *pp1);
        gradll1 *= 2;
        gradll2 *= 2;
        gradll3 *= 2;

        Vec3d gradll(gradll1);
        gradll += gradll2;
        gradll += gradll3;

        /*
        Vec3d gradll;
        gradll = v1+v2+v3;
        gradll *= -2;
         */

        err = 0.0080187537 * lll / vol;


        gradll *= (0.0080187537 * 1.5 * l / vol);
        Vec3d graderr(gradll);
        gradvol *= (-0.0080187537 * lll / (vol * vol));
        graderr += gradvol;

        if (h > 0) {
            /*
            Vec3d gradll1 (*pp2, *pp1);
            Vec3d gradll2 (*pp3, *pp1);
            Vec3d gradll3 (*pp4, *pp1);
            gradll1 *= 2;
            gradll2 *= 2;
            gradll3 *= 2;
             */
            err += ll / (h * h) +
                    h * h * (1 / ll1 + 1 / ll2 + 1 / ll3 +
                    1 / ll4 + 1 / ll5 + 1 / ll6) - 12;

            graderr += (1 / (h * h) - h * h / (ll1 * ll1)) * gradll1;
            graderr += (1 / (h * h) - h * h / (ll2 * ll2)) * gradll2;
            graderr += (1 / (h * h) - h * h / (ll3 * ll3)) * gradll3;
        }

        double errpow;

        double teterrpow = mp.opterrpow;
        if (teterrpow < 1) teterrpow = 1;

        if (teterrpow == 1) {
            errpow = err;
            grad = graderr;
        }
        else if (teterrpow == 2) {
            errpow = err*err;
            grad = (2 * err) * graderr;
        }
        else {

            errpow = pow(err, teterrpow);
            grad = (teterrpow * errpow / err) * graderr;
        }
        return errpow;
    }

    double CalcVolume(
            const Array<Point3d> & points,
            const Array<Element> & elements)
    {
        double vol;
        Vec3d v1, v2, v3;

        vol = 0;
        for (int i = 0; i < elements.size(); i++) {

            v1 = points.Get(elements[i][1]) - points.Get(elements[i][0]);
            v2 = points.Get(elements[i][2]) - points.Get(elements[i][0]);
            v3 = points.Get(elements[i][3]) - points.Get(elements[i][0]);
            vol -= (Cross(v1, v2) * v3) / 6;
        }
        return vol;
    }

    void MeshQuality3d(const Mesh & mesh, Array<int> * inclass)
    {
        int ncl = 20;
        signed int cl;
        Array<INDEX> incl(ncl);
        INDEX i;
        double qual;
        double sum = 0;
        int nontet = 0;

        for (i = 1; i <= incl.size(); i++)
            incl.Elem(i) = 0;

        for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++) {
            if (mesh[ei].GetType() != TET) {
                nontet++;
                continue;
            }

            qual = TetElementQuality(mesh.Point(mesh[ei][0]),
                    mesh.Point(mesh[ei][1]),
                    mesh.Point(mesh[ei][2]),
                    mesh.Point(mesh[ei][3]));

            if (qual > 1) qual = 1;
            cl = int (ncl * qual) + 1;

            if (cl < 1) cl = 1;
            if (cl > ncl) cl = ncl;

            incl.Elem(cl)++;
            if (inclass) (*inclass)[ei] = cl;
            sum += 1 / qual;
        }

        LOG_INFO("\n\n");
        LOG_INFO("Points:           " << mesh.GetNP());
        LOG_INFO("Volume Elements:  " << mesh.GetNE());
        if (nontet > 0)
            LOG_INFO(nontet << " non tetrahedral elements");
        LOG_INFO("Surface Elements: " << mesh.GetNSE());
        LOG_INFO("\nVolume elements in qualityclasses:");
        for (i = 1; i <= ncl; i++) {
            LOG_INFO(std::fixed << std::setprecision(2) <<
                    std::setw(4) << double (i - 1) / ncl << " - " <<
                    std::setw(4) << double (i) / ncl << ": " << incl.Get(i));
        }
        LOG_INFO("Total error: " << sum);
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

    void SaveVolumeMesh(const Mesh & mesh,
            const CSGeometry & geometry,
            char * filename)
    {
        INDEX i;

        std::ofstream outfile(filename);
        outfile << "volumemesh" << std::endl;

        outfile << mesh.GetNSE() << std::endl;
        for (i = 1; i <= mesh.GetNSE(); i++) {
            if (mesh.SurfaceElement(i).GetIndex())
                outfile << mesh.GetFaceDescriptor(mesh.SurfaceElement(i).GetIndex()).SurfNr()
                << "\t";
            else
                outfile << "0" << "\t";
            outfile << mesh.SurfaceElement(i)[0] << " "
                    << mesh.SurfaceElement(i)[1] << " "
                    << mesh.SurfaceElement(i)[2] << std::endl;
        }
        outfile << mesh.GetNE() << std::endl;
        for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
            outfile << mesh[ei].GetIndex() << "\t"
            << mesh[ei][0] << " " << mesh[ei][1] << " "
                << mesh[ei][2] << " " << mesh[ei][3] << std::endl;

        outfile << mesh.GetNP() << std::endl;
        for (i = 1; i <= mesh.GetNP(); i++)
            outfile << mesh.Point(i)(0) << " "
            << mesh.Point(i)(1) << " "
            << mesh.Point(i)(2) << std::endl;
    }


    /* ******************** CheckMesh ******************************* */

    /// Checks, whether mesh contains a valid 3d mesh

    int CheckMesh3D(const Mesh & mesh)
    {
        INDEX_3_HASHTABLE<int> faceused(mesh.GetNE() / 3);
        INDEX i;
        int j, k, l;
        INDEX_3 i3;
        int ok = 1;
        ElementIndex ei;

        for (i = 1; i <= mesh.GetNSE(); i++) {
            const Element2d & el = mesh.SurfaceElement(i);

            if (mesh.GetFaceDescriptor(el.GetIndex()).DomainIn() == 0 ||
                    mesh.GetFaceDescriptor(el.GetIndex()).DomainOut() == 0) {
                for (j = 1; j <= 3; j++)
                    i3.I(j) = el.PNum(j);

                i3.Sort();
                faceused.Set(i3, 1);
            }
        }

        for (ei = 0; ei < mesh.GetNE(); ei++) {
            const Element & el = mesh[ei];

            for (j = 1; j <= 4; j++) {
                l = 0;
                for (k = 1; k <= 4; k++) {
                    if (j != k) {
                        l++;
                        i3.I(l) = el.PNum(k);
                    }
                }

                i3.Sort();
                if (faceused.Used(i3))
                    faceused.Set(i3, faceused.Get(i3) + 1);
                else
                    faceused.Set(i3, 1);
            }
        }


        for (i = 1; i <= mesh.GetNSE(); i++) {
            const Element2d & el = mesh.SurfaceElement(i);

            for (j = 1; j <= 3; j++)
                i3.I(j) = el.PNum(j);

            i3.Sort();
            k = faceused.Get(i3);
            if (k != 2) {
                ok = 0;
                std::cerr << "face " << i << " with points "
                        << i3.I1() << "-" << i3.I2() << "-" << i3.I3()
                        << " has " << k << " elements" << std::endl;
            }
        }

        for (ei = 0; ei < mesh.GetNE(); ei++) {
            const Element & el = mesh[ei];

            for (j = 1; j <= 4; j++) {
                l = 0;
                for (k = 1; k <= 4; k++) {
                    if (j != k) {
                        l++;
                        i3.I(l) = el.PNum(k);
                    }
                }

                i3.Sort();
                k = faceused.Get(i3);
                if (k != 2) {
                    ok = 0;
                    std::cerr << "element " << ei << " with face "
                            << i3.I1() << "-" << i3.I2() << "-"
                            << i3.I3()
                            << " has " << k << " elements" << std::endl;
                }
            }
        }

        if (!ok) {
            std::cerr << "surfelements: " << std::endl;
            for (i = 1; i <= mesh.GetNSE(); i++) {
                const Element2d & el = mesh.SurfaceElement(i);
                std::cerr << std::setw(5) << i << ":"
                        << std::setw(6) << el.GetIndex()
                        << std::setw(6) << el.PNum(1)
                        << std::setw(4) << el.PNum(2)
                        << std::setw(4) << el.PNum(3) << std::endl;
            }
            std::cerr << "volelements: " << std::endl;
            for (ei = 0; ei < mesh.GetNE(); ei++) {

                const Element & el = mesh[ei];
                std::cerr << std::setw(5) << i << ":"
                        << std::setw(6) << el.GetIndex()
                        << std::setw(6) << el[0] << std::setw(4) << el[1]
                        << std::setw(4) << el[2] << std::setw(4) << el[3] << std::endl;
            }
        }

        return ok;
    }

    void RemoveProblem(Mesh & mesh, int domainnr)
    {
        int i, j, k;

        mesh.FindOpenElements(domainnr);
        int np = mesh.GetNP();

        BitArrayChar<PointIndex::BASE> ppoints(np);

        LOG_DEBUG("Elements before Remove: " << mesh.GetNE());
        k = domainnr;
        {
            ppoints.Clear();

            for (i = 1; i <= mesh.GetNOpenElements(); i++) {
                const Element2d & sel = mesh.OpenElement(i);
                if (sel.GetIndex() == k) {
                    for (j = 1; j <= sel.GetNP(); j++)
                        ppoints.Set(sel.PNum(j));
                }
            }

            for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++) {
                const Element & el = mesh[ei];
                if (el.GetIndex() == k) {
                    int todel = 0;
                    for (j = 0; j < el.GetNP(); j++)
                        if (ppoints.Test(el[j]))
                            todel = 1;

                    if (el.GetNP() != 4)
                        todel = 0;

                    if (todel) {
                        mesh[ei].Delete();
                    }
                }
            }
        }

        mesh.Compress();
        LOG_DEBUG("Elements after Remove: " << mesh.GetNE());
    }

}
