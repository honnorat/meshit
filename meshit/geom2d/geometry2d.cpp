/*

2d Spline curve for Mesh generator

 */

#include "geometry2d.hpp"

#include <string>
#include <vector>

#include "../general/flags.hpp"

namespace meshit {

    SplineGeometry2d::~SplineGeometry2d()
    {
        for (size_t i = 0; i < materials.size(); i++) {
            delete[] materials[i];
        }
    }

    void SplineGeometry2d::Load(const std::string& filename)
    {
        std::ifstream infile;
        Point<2> x;
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

    void SplineGeometry2d::TestComment(std::istream& infile)
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

    void SplineGeometry2d::LoadData(std::istream& infile)
    {
        MESHIT_LOG_INFO("Load 2D Geometry");
        int nump, leftdom, rightdom;
        Point<2> x;
        int hi1, hi2, hi3;
        double hd;
        char buf[50], ch;
        int pointnr;

        std::string keyword;

        std::vector<GeomPoint<2> > infilepoints;
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
                    if (infile.good())
                        infile.putback(ch);

                    if (hd == 1)
                        hd = flags.GetNumFlag("ref", 1.0);

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

                    SplineSeg* spline = 0;
                    TestComment(infile);

                    infile >> leftdom >> rightdom;

                    if (leftdom > numdomains) numdomains = leftdom;
                    if (rightdom > numdomains) numdomains = rightdom;

                    infile >> buf;
                    // type of spline segement
                    if (strcmp(buf, "2") == 0) {  // a line
                        infile >> hi1 >> hi2;
                        spline = new LineSeg(
                                geompoints[hi1 - 1],
                                geompoints[hi2 - 1]);
                    }
                    else if (strcmp(buf, "3") == 0) {  // a rational spline
                        infile >> hi1 >> hi2 >> hi3;
                        spline = new SplineSeg3(
                                geompoints[hi1 - 1],
                                geompoints[hi2 - 1],
                                geompoints[hi3 - 1]);
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

                    if (infile.good())
                        infile.putback(ch);

                    spex->bc = static_cast<int>(flags.GetNumFlag("bc", i + 1));
                    spex->copyfrom = static_cast<int>(flags.GetNumFlag("copy", -1));
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

                if (!infile.good())
                    return;

                materials.resize(numdomains);
                maxh.resize(numdomains);
                for (int i = 0; i < numdomains; i++) {
                    maxh[i] = 1000;
                }
                layer.resize(numdomains, 1);

                TestComment(infile);

                for (int i = 0; i < numdomains; i++) {
                    materials[i] = new char[100];
                }
                for (int i = 0; i < numdomains && infile.good(); i++) {
                    TestComment(infile);
                    infile >> domainnr;
                    infile >> material;

                    strcpy(materials[domainnr - 1], material);

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

                    if (infile.good())
                        infile.putback(ch);

                    maxh[domainnr - 1] = flags.GetNumFlag("maxh", 1000);
                    layer[domainnr - 1] = static_cast<int>(flags.GetNumFlag("layer", 1));
                }
            }
        }
        return;
    }

    void SplineGeometry2d::AddLine(const std::vector<Point2d>& point_list,
                                   double hmax,
                                   bool hole,
                                   int bc)
    {
        size_t nold_points = geompoints.size();
        size_t nnew_points = point_list.size();

        std::vector<GeomPoint<2> > gpts;

        gpts.reserve(nnew_points);
        geompoints.reserve(nold_points + nnew_points);

        for (size_t i = 0; i < nnew_points; i++) {
            gpts.push_back(GeomPoint<2>(point_list[i]));
            geompoints.push_back(gpts[i]);
        }
        for (size_t i = 0; i < nnew_points; i++) {
            size_t i0 = (hole) ? nnew_points - i - 1 : i;
            size_t i1 = (hole) ? (i0 == 0) ? nnew_points - 1 : i0 - 1
                               : (i0 == nnew_points - 1) ? 0 : i0 + 1;
            SplineSeg* spline = new LineSeg(gpts[i0], gpts[i1]);
            SplineSegExt* seg = new SplineSegExt(*spline);
            if (hole) {
                seg->leftdom = 0;
                seg->rightdom = 1;
            } else {
                seg->leftdom = 1;
                seg->rightdom = 0;
            }
            seg->bc = bc;
            seg->reffak = 1;  // Refinement factor
            seg->hmax = hmax;
            splines.push_back(seg);
        }
    }

    void SplineGeometry2d::AddStructureLine(const std::vector<Point2d>& point_list,
                                            double hmax,
                                            int bc)
    {

        size_t nold_points = geompoints.size();
        size_t nnew_points = point_list.size();

        std::vector<GeomPoint<2> > gpts;

        gpts.reserve(nnew_points);
        geompoints.reserve(nold_points + nnew_points);

        for (size_t i = 0; i < nnew_points; i++) {
            gpts.push_back(GeomPoint<2>(point_list[i]));
            geompoints.push_back(gpts[i]);
        }
        for (size_t i = 0; i < nnew_points; i++) {
            size_t i0 = i;
            size_t i1 = (i0 == nnew_points - 1) ? 0 : i0 + 1;
            SplineSeg* spline = new LineSeg(gpts[i0], gpts[i1]);
            SplineSegExt* seg = new SplineSegExt(*spline);
            seg->leftdom = 1;
            seg->rightdom = 1;
            seg->bc = bc;
            seg->reffak = 2;  // Refinement factor
            seg->hmax = hmax;
            splines.push_back(seg);
        }
    }

    void SplineGeometry2d::FakeData()
    {
        int numdomains = 1;
        materials.resize(numdomains);
        maxh.assign(numdomains, 1000);
        layer.assign(numdomains, 1);

        for (int i = 0; i < numdomains; i++) {
            materials[i] = new char[1];
            materials[i][0] = '\0';
        }
    }

    void SplineGeometry2d::GetMaterial(const int domnr, char*& material)
    {
        if ((int) materials.size() >= domnr)
            material = materials[domnr - 1];
        else
            material = 0;
    }

    double SplineGeometry2d::GetDomainMaxh(const int domnr)
    {
        if ((int) maxh.size() >= domnr && domnr > 0)
            return maxh[domnr - 1];
        else
            return -1;
    }

    int SplineGeometry2d::GenerateMesh(Mesh*& mesh, MeshingParameters& mp)
    {
        mesh->BuildFromSpline2D(*this, mp);
        return 0;
    }

}  // namespace meshit
