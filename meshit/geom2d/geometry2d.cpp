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

    void SplineGeometry::GetBoundingBox(Box2d& box) const
    {
        if (!splines.size()) {
            Point2d auxp{0.0, 0.0};
            box.SetPoint(auxp);
            return;
        }

        std::vector<Point2d> points;
        for (size_t i = 0; i < splines.size(); i++) {
            splines[i]->GetPoints(20, points);

            if (i == 0) box.SetPoint(points[0]);
            for (size_t j = 0; j < points.size(); j++) {
                box.AddPoint(points[j]);
            }
        }
    }

    Box2d SplineGeometry::GetBoundingBox() const
    {
        Box2d box;
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

    char SplineGeometry::TestComment(std::istream& infile)
    {
        bool comment = true;
        char ch = '\0';
        while (comment && !infile.eof()) {
            infile >> ch;
            if (ch == '#') {
                // skip comments
                while (ch != '\n' && !infile.eof()) {
                    infile.get(ch);
                }
            } else if (ch == '\n') {
                // skip empty lines
            } else if (isblank(ch)) {
                // skip whitespaces
            } else if (ch == '-') {
                // do not unget '-'
                comment = false;
            } else {
                // end of comment
                infile.unget();
                comment = false;
            }
        }
        return ch;
    }

    void SplineGeometry::LoadData(std::istream& infile)
    {
        MESHIT_LOG_INFO("Load 2D Geometry");
        int nump, dom_left, dom_right;
        Point2d x;
        int hi1, hi2, hi3;
        char buf[50], ch;
        int pointnr;

        std::string keyword;
        std::string flag;

        std::vector<GeomPoint> infilepoints;
        std::vector<int> pointnrs;
        nump = 0;
        int numdomains = 0;

        TestComment(infile);
        // refinement factor
        infile >> elto0;

        while (infile.good()) {
            TestComment(infile);
            infile >> keyword;
            ch = TestComment(infile);

            if (keyword == "points") {
                while (!isalpha(static_cast<int>(ch))) {
                    infile >> pointnr;  // pointnrs is 1-based
                    if (pointnr > nump) nump = pointnr;
                    pointnrs.push_back(pointnr);

                    infile >> x.X() >> x.Y();
                    infile >> ch;

                    Flags flags;
                    while (ch == '-') {
                        infile >> flag;
                        flags.SetCommandLineFlag(flag);
                        ch = TestComment(infile);
                    }
                    infile.unget();

                    infilepoints.push_back(GeomPoint(x,
                                                     flags.GetNumFlag("ref", 1.0),
                                                     flags.GetNumFlag("maxh", 1e99)));
                    ch = TestComment(infile);
                }
                geompoints.resize(nump);
                for (int i = 0; i < nump; i++) {
                    geompoints[pointnrs[i] - 1] = infilepoints[i];
                }
            }
            else if (keyword == "segments") {
                int i = 1;
                while (!isalpha(static_cast<int>(ch))) {
                    infile >> dom_left >> dom_right;
                    if (dom_left > numdomains) numdomains = dom_left;
                    if (dom_right > numdomains) numdomains = dom_right;

                    SplineSeg* spline = nullptr;
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

                    spline->dom_left = dom_left;
                    spline->dom_right = dom_right;
                    splines.push_back(spline);

                    Flags flags;

                    infile >> ch;
                    while (ch == '-') {
                        infile >> flag;
                        flags.SetCommandLineFlag(flag);
                        ch = TestComment(infile);
                    }
                    infile.unget();

                    spline->id_ = static_cast<int>(flags.GetNumFlag("id", ++i));
                    spline->reffak = flags.GetNumFlag("ref", 1);
                    spline->hmax = flags.GetNumFlag("maxh", 1e99);
                    ch = TestComment(infile);
                }
            }
            else if (keyword == "materials") {
                int domainnr;
                char material[100];

                materials.resize(numdomains);
                maxh.resize(numdomains);
                for (int i = 0; i < numdomains; i++) {
                    maxh[i] = 1000;
                }

                for (int i = 0; i < numdomains; i++) {
                    materials[i] = new char[100];
                }
                for (int i = 0; i < numdomains && infile.good(); i++) {
                    infile >> domainnr;
                    infile >> material;
                    strncpy(materials[domainnr - 1], material, 100);

                    Flags flags;
                    infile >> ch;
                    while (ch == '-') {
                        infile >> flag;
                        flags.SetCommandLineFlag(flag);
                        ch = TestComment(infile);
                    }
                    infile.unget();
                    maxh[domainnr - 1] = flags.GetNumFlag("maxh", 1000);
                    ch = TestComment(infile);
                }
            }
        }
        return;
    }

    void SplineGeometry::AddLine(
        const std::vector<Point2d>& points, double hmax, int spline_id, int domain_left, int domain_right)
    {
        size_t nold_points = geompoints.size();
        size_t nnew_points = points.size();

        std::vector<GeomPoint> gpts;

        gpts.reserve(nnew_points);
        geompoints.reserve(nold_points + nnew_points);

        for (size_t i = 0; i < nnew_points; i++) {
            gpts.push_back(GeomPoint(points[i]));
            geompoints.push_back(gpts[i]);
        }
        for (size_t i0 = 0; i0 < nnew_points; i0++) {
            size_t i1 = (i0 == nnew_points - 1) ? 0 : i0 + 1;
            SplineSeg* spline = new LineSeg(gpts[i0], gpts[i1]);
            spline->dom_left = domain_left;
            spline->dom_right = domain_right;
            spline->id_ = spline_id;
            spline->reffak = 1;  // Refinement factor
            spline->hmax = hmax;
            splines.push_back(spline);
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
                                   double hmax, int spline_id,
                                   int domain_left,
                                   int domain_right)
    {
        size_t nb_points = points.size();
        size_t ip_0 = geompoints.size();
        size_t ip = ip_0;

        if (nb_points < 4 || nb_points % 2 > 0) {
            throw std::runtime_error("SplineGeometry::AddSpline : wrong number of points.");
        }
        for (size_t i = 0; i < nb_points; i++) {
            geompoints.push_back(meshit::GeomPoint(points[i]));
        }

        size_t nb_splines = nb_points / 2;

        for (size_t i = 0; i < nb_splines; i++) {
            size_t id0 = ip + 0;
            size_t id1 = ip + 1;
            size_t id2 = ip + 2;
            if (i == nb_splines - 1) id2 = ip_0;
            SplineSeg* spline = new SplineSeg3(geompoints[id0], geompoints[id1], geompoints[id2]);
            spline->dom_left = domain_left;
            spline->dom_right = domain_right;
            spline->id_ = spline_id;
            spline->reffak = 1;
            spline->hmax = hmax;
            splines.push_back(spline);
            ip += 2;
        }
    }

    void SplineGeometry::AddCircle(const Point2d& center, double radius,
                                   double hmax, int spline_id,
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

        AddSpline(spline_points, hmax, spline_id, face_left, face_right);
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

    double SplineGeometry::GetDomainMaxh(size_t domain_id)
    {
        if (domain_id > 0 && domain_id <= maxh.size()) {
            return maxh[domain_id - 1];
        } else {
            return -1.0;
        }
    }

    void SplineSegmenter::Partition(const SplineSeg& spline)
    {
        constexpr size_t n = 10000;
        constexpr double dt = 1.0 / n;

        std::vector<double> curve_points;
        CalcPartition(spline, curve_points);

        std::vector<size_t> loc_search;
        Point2d p_old = spline.GetPoint(0);
        Point2d mark_old = p_old;
        double l_old = 0.0;
        size_t j = 1;

        for (size_t i = 1; i <= n; i++) {
            double t = static_cast<double>(i) * dt;
            Point2d p = spline.GetPoint(t);
            double l = l_old + Dist(p, p_old);
            while (j < curve_points.size() && (l >= curve_points[j] || i == n)) {
                PointIndex pi1, pi2;
                double frac = (curve_points[j] - l) / (l - l_old);
                Point2d mark = spline.GetPoint(t + frac * dt);

                double h = mesh_.GetH(mark_old);
                Vec2d v(1e-4 * h, 1e-4 * h);
                searchtree_.GetIntersecting(mark_old - v, mark_old + v, loc_search);
                if (loc_search.size() > 0) {
                    pi1 = loc_search.back();
                } else {
                    pi1 = mesh_.AddPoint(mark_old);
                    searchtree_.Insert(mark_old, pi1);
                }
                searchtree_.GetIntersecting(mark - v, mark + v, loc_search);
                if (loc_search.size() > 0) {
                    pi2 = loc_search.back();
                } else {
                    pi2 = mesh_.AddPoint(mark);
                    searchtree_.Insert(mark, pi2);
                }

                Segment seg;
                seg.edge_id = spline.get_id();
                seg[0] = pi1;
                seg[1] = pi2;
                seg.dom_left = spline.dom_left;
                seg.dom_right = spline.dom_right;
                mesh_.AddSegment(seg);

                mark_old = mark;
                j++;
            }

            p_old = p;
            l_old = l;
        }
    }

    void SplineSegmenter::CalcPartition(const SplineSeg& spline, std::vector<double>& points)
    {
        double fperel, oldf, f;

        size_t n = 10000;

        std::vector<Point2d> xi(n);
        std::vector<double> hi(n);

        for (size_t i = 0; i < n; i++) {
            xi[i] = spline.GetPoint((i + 0.5) / n);
            hi[i] = mesh_.GetH(xi[i]);
        }

        // limit slope
        double gradh = 1 / elto0_;
        for (size_t i = 0; i < n - 1; i++) {
            double hnext = hi[i] + gradh * (xi[i + 1] - xi[i]).Length();
            hi[i + 1] = std::min(hi[i + 1], hnext);
        }
        for (size_t i = n - 1; i > 1; i--) {
            double hnext = hi[i] + gradh * (xi[i - 1] - xi[i]).Length();
            hi[i - 1] = std::min(hi[i - 1], hnext);
        }

        points.clear();

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


}  // namespace meshit
