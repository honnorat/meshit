/*

2d Spline curve for Mesh generator

 */

#include "geometry2d.hpp"

#include "../general/flags.hpp"

namespace meshit
{
    SplineGeometry::~SplineGeometry()
    {
        for (size_t i = 0; i < splines.size(); i++) {
            delete splines[i];
        }
        for (size_t i = 0; i < materials.size(); i++) {
            delete[] materials[i];
        }
    }

    void SplineGeometry::GetBoundingBox(Box<2>& box) const
    {
        if (!splines.size()) {
            Point<2> auxp{0.0};
            box.Set(auxp);
            return;
        }

        std::vector<Point<2> > points;
        for (size_t i = 0; i < splines.size(); i++) {
            splines[i]->GetPoints(20, points);

            if (i == 0) box.Set(points[0]);
            for (size_t j = 0; j < points.size(); j++) {
                box.Add(points[j]);
            }
        }
    }

    Box<2> SplineGeometry::GetBoundingBox() const
    {
        Box<2> box;
        GetBoundingBox(box);
        return box;
    }

    void SplineGeometry::Load(const std::string& filename)
    {
        std::ifstream infile;
        char buf[50];

        infile.open(filename);

        if (!infile.good()) {
            std::string message = "Input file '" + filename + "' not available!";
            throw std::runtime_error(message);
        }

        TestComment(infile);

        infile >> buf;  // file recognition

        TestComment(infile);
        if (strcmp(buf, "splinecurves2dv2") == 0) {
            LoadData(infile);
        } else {
            MESHIT_LOG_FATAL("Unsupported file format : '" << buf << "'");
            throw std::runtime_error("Unsupported file format");
        }
        infile.close();
    }

    void SplineGeometry::TestComment(std::istream& infile)
    {
        bool comment = true;
        char ch;
        while (comment && !infile.eof()) {
            infile.get(ch);
            if (ch == '#') {
                // skip comments
                while (ch != '\n' && !infile.eof()) {
                    infile.get(ch);
                }
            } else if (ch == '\n') {
                // skip empty lines
            } else if (isspace(ch)) {
                // skip whitespaces
            } else {
                // end of comment
                infile.putback(ch);
                comment = false;
            }
        }
        return;
    }

    void SplineGeometry::LoadData(std::istream& infile)
    {
        MESHIT_LOG_INFO("Load 2D Geometry");
        int nump, leftdom, rightdom;
        Point<2> x;
        int hi1, hi2, hi3;
        double hd;
        char buf[50], ch;
        int pointnr;

        std::string keyword;

        std::vector<GeomPoint<2>> infilepoints;
        std::vector<int> pointnrs;
        nump = 0;
        int numdomains = 0;

        TestComment(infile);
        // refinement factor
        infile >> elto0;
        TestComment(infile);

        while (infile.good()) {
            infile >> keyword;

            if (keyword == "points") {
                MESHIT_LOG_DEBUG("load points");
                infile.get(ch);
                infile.putback(ch);

                while (!isalpha(static_cast<int>(ch))) {
                    TestComment(infile);
                    infile >> pointnr;
                    // pointnrs 1-based
                    if (pointnr > nump) nump = pointnr;
                    pointnrs.push_back(pointnr);

                    for (int j = 0; j < 2; j++) {
                        infile >> x[j];
                    }
                    // hd is now optional, default 1
                    //  infile >> hd;
                    hd = 1;

                    Flags flags;

                    // get flags,
                    ch = 'a';

                    do {
                        infile.get(ch);
                        // if another int-value, set refinement flag to this value
                        // (corresponding to old files)
                        if (isdigit(static_cast<int>(ch))) {
                            infile.putback(ch);
                            infile >> hd;
                            infile.get(ch);
                        }
                    } while (isspace(ch) && ch != '\n');
                    while (ch == '-') {
                        char flag[100];
                        flag[0] = '-';
                        infile >> (flag + 1);
                        flags.SetCommandLineFlag(flag);
                        ch = 'a';
                        do {
                            infile.get(ch);
                        } while (isspace(ch) && ch != '\n');
                    }
                    if (infile.good()) infile.putback(ch);

                    if (hd == 1) hd = flags.GetNumFlag("ref", 1.0);

                    infilepoints.push_back(GeomPoint<2>(x, hd, flags.GetNumFlag("maxh", 1e99)));

                    TestComment(infile);
                    infile.get(ch);
                    infile.putback(ch);
                }

                geompoints.resize(nump);
                for (int i = 0; i < nump; i++) {
                    geompoints[pointnrs[i] - 1] = infilepoints[i];
                }
                TestComment(infile);
            }

            else if (keyword == "segments") {
                MESHIT_LOG_DEBUG("load segments");

                infile.get(ch);
                infile.putback(ch);
                int i = 0;

                while (!isalpha(static_cast<int>(ch))) {
                    i++;
                    TestComment(infile);

                    SplineSeg* spline = nullptr;
                    TestComment(infile);

                    infile >> leftdom >> rightdom;

                    if (leftdom > numdomains) numdomains = leftdom;
                    if (rightdom > numdomains) numdomains = rightdom;

                    infile >> buf;
                    // type of spline segement
                    if (strcmp(buf, "2") == 0) {  // a line
                        infile >> hi1 >> hi2;
                        spline = new LineSeg(geompoints[hi1 - 1],
                                             geompoints[hi2 - 1]);
                    } else if (strcmp(buf, "3") == 0) {  // a rational spline
                        infile >> hi1 >> hi2 >> hi3;
                        spline = new SplineSeg3(geompoints[hi1 - 1],
                                                geompoints[hi2 - 1],
                                                geompoints[hi3 - 1]);
                    } else {
                        MESHIT_LOG_ERROR("Unknown segment type : " << buf);
                        throw std::runtime_error("SplineGeometry::LoadData : unknown segment type");
                    }

                    SplineSegExt* spex = new SplineSegExt(*spline);

                    spex->leftdom = leftdom;
                    spex->rightdom = rightdom;
                    splines.push_back(spex);

                    // hd is now optional, default 1
                    hd = 1;
                    infile >> ch;

                    // get flags
                    Flags flags;
                    while (ch == '-') {
                        char flag[100];
                        flag[0] = '-';
                        infile >> (flag + 1);
                        flags.SetCommandLineFlag(flag);
                        ch = 'a';
                        infile >> ch;
                    }

                    if (infile.good()) infile.putback(ch);

                    spex->bc = static_cast<int>(flags.GetNumFlag("bc", i + 1));
                    spex->reffak = flags.GetNumFlag("ref", 1);
                    spex->hmax = flags.GetNumFlag("maxh", 1e99);
                    if (hd != 1) spex->reffak = hd;

                    TestComment(infile);
                    infile.get(ch);
                    infile.putback(ch);
                }
                infile.get(ch);
                infile.putback(ch);
            }

            else if (keyword == "materials") {
                TestComment(infile);
                int domainnr;
                char material[100];

                if (!infile.good()) return;

                materials.resize(numdomains);
                maxh.resize(numdomains);
                for (int i = 0; i < numdomains; i++) {
                    maxh[i] = 1000;
                }

                TestComment(infile);

                for (int i = 0; i < numdomains; i++) {
                    materials[i] = new char[100];
                }
                for (int i = 0; i < numdomains && infile.good(); i++) {
                    TestComment(infile);
                    infile >> domainnr;
                    infile >> material;

                    strncpy(materials[domainnr - 1], material, 100);

                    Flags flags;
                    ch = 'a';
                    infile >> ch;
                    while (ch == '-') {
                        char flag[100];
                        flag[0] = '-';
                        infile >> (flag + 1);
                        flags.SetCommandLineFlag(flag);
                        ch = 'a';
                        infile >> ch;
                    }

                    if (infile.good()) infile.putback(ch);

                    maxh[domainnr - 1] = flags.GetNumFlag("maxh", 1000);
                }
            }
        }
        return;
    }

    void SplineGeometry::AddLine(
        const std::vector<Point2d>& points, double hmax, int bc, int face_left, int face_right)
    {
        size_t nold_points = geompoints.size();
        size_t nnew_points = points.size();

        std::vector<GeomPoint<2>> gpts;

        gpts.reserve(nnew_points);
        geompoints.reserve(nold_points + nnew_points);

        for (size_t i = 0; i < nnew_points; i++) {
            gpts.push_back(GeomPoint<2>(points[i]));
            geompoints.push_back(gpts[i]);
        }
        for (size_t i0 = 0; i0 < nnew_points; i0++) {
            size_t i1 = (i0 == nnew_points - 1) ? 0 : i0 + 1;
            SplineSeg* spline = new LineSeg(gpts[i0], gpts[i1]);
            SplineSegExt* seg = new SplineSegExt(*spline);
            seg->leftdom = face_left;
            seg->rightdom = face_right;
            seg->bc = bc;
            seg->reffak = 1;  // Refinement factor
            seg->hmax = hmax;
            splines.push_back(seg);
        }
    }

    void SplineGeometry::AddHole(const std::vector<Point2d>& point_list, double hmax, int bc, int face)
    {
        AddLine(point_list, hmax, bc, 0, face);
    }

    void SplineGeometry::AddStructureLine(const std::vector<Point2d>& points, double hmax, int bc, int face)
    {
        AddLine(points, hmax, bc, face, face);
    }

    void SplineGeometry::AddSpline(const std::vector<Point2d>& points,
                                     double hmax, int bc,
                                     int face_left,
                                     int face_right)
    {
        size_t nb_points = points.size();
        size_t ip_0 = geompoints.size();
        size_t ip = ip_0;

        if (nb_points < 4 || nb_points % 2 > 0) {
            throw std::runtime_error("SplineGeometry::AddSpline : wrong number of points.");
        }
        for (size_t i = 0; i < nb_points; i++) {
            geompoints.push_back(meshit::GeomPoint<2>(points[i]));
        }

        size_t nb_splines = nb_points / 2;

        for (size_t i = 0; i < nb_splines; i++) {
            size_t id0 = ip + 0;
            size_t id1 = ip + 1;
            size_t id2 = ip + 2;
            if (i == nb_splines - 1) id2 = ip_0;
            SplineSeg* spline = new SplineSeg3(geompoints[id0], geompoints[id1], geompoints[id2]);
            SplineSegExt* spex = new SplineSegExt(*spline);
            spex->leftdom = face_left;
            spex->rightdom = face_right;
            spex->bc = bc;
            spex->reffak = 1;
            spex->hmax = hmax;
            splines.push_back(spex);
            ip += 2;
        }
    }

    void SplineGeometry::AddCircle(const Point2d& center, double radius,
                                     double hmax, int bc,
                                     int face_left, int face_right)
    {
        std::vector<Point2d> spline_points;
        double c_x = center.X();
        double c_y = center.Y();

        spline_points.reserve(8);
        spline_points.push_back(Point2d(c_x + radius, c_y));
        spline_points.push_back(Point2d(c_x + radius, c_y + radius));
        spline_points.push_back(Point2d(c_x, c_y + radius));
        spline_points.push_back(Point2d(c_x - radius, c_y + radius));
        spline_points.push_back(Point2d(c_x - radius, c_y));
        spline_points.push_back(Point2d(c_x - radius, c_y - radius));
        spline_points.push_back(Point2d(c_x, c_y - radius));
        spline_points.push_back(Point2d(c_x + radius, c_y - radius));

        AddSpline(spline_points, hmax, bc, face_left, face_right);
    }

    void SplineGeometry::FakeData()
    {
        int numdomains = 1;
        materials.resize(numdomains);
        maxh.assign(numdomains, 1000);

        for (int i = 0; i < numdomains; i++) {
            materials[i] = new char[1];
            materials[i][0] = '\0';
        }
    }

    int SplineGeometry::AddFace(const std::string& name, double maxh_f)
    {
        maxh.push_back(maxh_f);
        materials.push_back(new char[100]);
        strncpy(materials.back(), name.c_str(), 100);

        return materials.size();
    }

    void SplineGeometry::GetMaterial(size_t domnr, char*& material)
    {
        if (domnr <= materials.size())
            material = materials[domnr - 1];
        else
            material = nullptr;
    }

    double SplineGeometry::GetDomainMaxh(size_t domnr)
    {
        if (domnr > 0 && domnr <= maxh.size()) {
            return maxh[domnr - 1];
        } else {
            return -1.0;
        }
    }

    int SplineGeometry::GenerateMesh(Mesh*& mesh, MeshingParameters& mp)
    {
        mesh->BuildFromSplineGeometry(*this, mp);
        return 0;
    }

}  // namespace meshit
