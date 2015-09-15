#include <meshit.hpp>

#include <stdexcept>

#include <meshit/meshing/meshclass.hpp>
#include <meshit/meshing/meshtool.hpp>
#include <meshit/geom2d/geometry2d.hpp>

#include <sstream>

struct CubicPoly
{
    double c0, c1, c2, c3;

    double eval(double t)
    {
        double t2 = t*t;
        double t3 = t * t* t;
        return c0 + c1 * t + c2 * t2 + c3*t3;
    }
    // compute coefficients for a nonuniform Catmull-Rom spline

    void InitNonuniformCatmullRom(double x0, double x1, double x2, double x3, double dt0, double dt1, double dt2)
    {
        // CF. http://www.cemyuksel.com/research/catmullrom_param/catmullrom.pdf
        // compute tangents when parameterized in [t1,t2]
        double t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1;
        double t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2;

        // rescale tangents for parametrization in [0,1]
        t1 *= dt1;
        t2 *= dt1;

        // Compute coefficients for a cubic polynomial
        //    p(s) = c0 + c1*s + c2*s^2 + c3*s^3
        // such that  p(0) = x1, p(1) = x2  and  p'(0) = t1, p'(1) = t2.
        c0 = x1;
        c1 = t1;
        c2 = -3 * x1 + 3 * x2 - 2 * t1 - t2;
        c3 = 2 * x1 - 2 * x2 + t1 + t2;
    }
};

class BoundaryLine
{
  public:
    typedef meshit::Point2d Point2d;

  protected:
    Point2d _center;
    std::vector<Point2d> _bbox;
    std::vector<Point2d> _points_cart;
    std::vector<Point2d> _points_polar;

  public:

    BoundaryLine() { }

    BoundaryLine(const BoundaryLine& other) :
        _center(other._center),
        _bbox(other._bbox),
        _points_cart(other._points_cart),
        _points_polar(other._points_polar) { }

    BoundaryLine(const std::string& filename)
    {
        read_boundary_file(filename);
    }

    ~BoundaryLine() { }

    size_t size() const
    {
        return _points_cart.size();
    }

    const Point2d& cart(size_t i) const
    {
        return _points_cart[i];
    }

    const Point2d& polar(size_t i) const
    {
        return _points_polar[i];
    }

    const std::vector<Point2d>& cart() const
    {
        return _points_cart;
    }

    std::vector<Point2d>& cart()
    {
        return _points_cart;
    }

    void reverse()
    {
        std::reverse(std::begin(_points_cart), std::end(_points_cart));
        std::reverse(std::begin(_points_polar), std::end(_points_polar));
    }

    void read_boundary_file(const std::string& filename)
    {
        std::ifstream boundary_file(filename.c_str());
        if (!boundary_file.good()) {
            std::stringstream msg;
            msg << "file '" << filename << "' not found !";
            throw std::runtime_error(msg.str());
        }

        _points_cart.clear();
        std::string line;
        while (std::getline(boundary_file, line)) {
            std::istringstream iss(line);
            double x, y;
            if (!(iss >> x >> y)) {
                continue;
            }
            _points_cart.push_back(Point2d(x, y));
        }
        calc_bbox();
        to_polar();
    }

    void calc_bbox()
    {
        size_t npoints = _points_cart.size();
        double x_cen = 0.0;
        double y_cen = 0.0;
        double x_min = std::numeric_limits<double>::max();
        double y_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::lowest();
        double y_max = std::numeric_limits<double>::lowest();

        for (size_t i = 0; i < npoints; i++) {
            double x = _points_cart[i].X();
            double y = _points_cart[i].Y();
            x_cen += x;
            y_cen += y;
            if (x < x_min) x_min = x;
            if (x > x_max) x_max = x;
            if (y < y_min) y_min = y;
            if (y > y_max) y_max = y;
        }
        _center.X() = x_cen / npoints;
        _center.Y() = y_cen / npoints;

        _bbox.push_back(Point2d(x_min, y_min));
        _bbox.push_back(Point2d(x_max, y_min));
        _bbox.push_back(Point2d(x_max, y_max));
        _bbox.push_back(Point2d(x_min, y_max));
    }

    void to_polar()
    {
        size_t nb_points = _points_cart.size();

        _points_polar.clear();
        _points_polar.reserve(nb_points);

        for (size_t i = 0; i < nb_points; i++) {
            double x = _points_cart[i].X() - _center.X();
            double y = _points_cart[i].Y() - _center.Y();
            double r = sqrt(x * x + y * y);
            double t = atan2(y, x);
            _points_polar.push_back(Point2d(r, t));
        }
    }

    void to_cart()
    {
        size_t nb_points = _points_polar.size();

        _points_cart.clear();
        _points_cart.reserve(nb_points);

        for (size_t i = 0; i < nb_points; i++) {
            double r = _points_polar[i].X();
            double t = _points_polar[i].Y();
            double x = _center.X() + r * cos(t);
            double y = _center.Y() + r * sin(t);
            _points_cart.push_back(Point2d(x, y));
        }
    }

    BoundaryLine subsample(unsigned int sub)
    {
        if (sub < 2) {
            return *this;
        }

        BoundaryLine out;
        out._points_cart.reserve(_points_cart.size() / sub + 1);
        for (size_t i = 0; i < _points_cart.size(); i += sub) {
            out._points_cart.push_back(_points_cart[i]);
        }
        out.calc_bbox();
        out.to_polar();
        return out;
    }

    void resample(size_t nb_points)
    {
        if (nb_points < 3) {
            throw std::runtime_error("BoundaryLine::resample : you must use at least 3 points.");
        }

        // Prepare polar data
        to_polar();
        std::vector<double> new_angles;
        new_angles.reserve(nb_points);

        const double CST_2_PI = 2 * M_PI;
        for (size_t i = 0; i < nb_points; i++) {
            double new_t = _points_polar[0].Y() + CST_2_PI * double(i) / nb_points;
            if (new_t > M_PI)
                new_t -= CST_2_PI;
            if (i > 0) {
                double dtheta = new_t - new_angles[i - 1];
                if (dtheta < 0) dtheta += CST_2_PI;
            }
            new_angles.push_back(new_t);
        }

        // Handle first point
        BoundaryLine newline;
        newline.cart().reserve(nb_points);
        newline.cart().push_back(_points_cart[0]);

        // Handle all other points
        size_t old_size = _points_cart.size();
        size_t i0 = old_size - 1;
        size_t i1 = 0, i2 = 1, i3 = 2;

        for (size_t i = 1; i < nb_points; i++) {

            double theta_test = new_angles[i];
            double theta_p1 = _points_polar[i1].Y();
            double theta_p2 = _points_polar[i2].Y();

            if (theta_p2 <= theta_p1) {
                theta_p2 += CST_2_PI;
            }

            bool test = ((theta_p1 < theta_test && theta_test <= theta_p2) ||
                    (theta_test < M_PI && theta_test < theta_p1
                    && (theta_p1 < theta_test + CST_2_PI && theta_test + CST_2_PI <= theta_p2)));

            // Find index of source point
            while (!test) {
                i0++;
                i1++;
                i2++;
                i3++;
                if (i0 == old_size) i0 = 0;
                if (i1 == old_size) i1 = 0;
                if (i2 == old_size) i2 = 0;
                if (i3 == old_size) i3 = 0;

                theta_p1 = _points_polar[i1].Y();
                theta_p2 = _points_polar[i2].Y();

                if (theta_p2 < theta_p1) {
                    theta_p2 += CST_2_PI;
                }
                test = ((theta_p1 < theta_test && theta_test <= theta_p2) ||
                        (theta_test < M_PI && theta_test < theta_p1
                        && (theta_p1 < theta_test + CST_2_PI && theta_test + CST_2_PI <= theta_p2)));
            }

            const Point2d& p0 = _points_cart[i0];
            const Point2d& p1 = _points_cart[i1];
            const Point2d& p2 = _points_cart[i2];
            const Point2d& p3 = _points_cart[i3];

            // Init centripetal Catmull-Rom interpolator
            double dt0 = pow(meshit::Dist2(p0, p1), 0.25);
            double dt1 = pow(meshit::Dist2(p1, p2), 0.25);
            double dt2 = pow(meshit::Dist2(p2, p3), 0.25);

            // safety check for repeated points
            if (dt1 < 1e-4) dt1 = 1.0;
            if (dt0 < 1e-4) dt0 = dt1;
            if (dt2 < 1e-4) dt2 = dt1;

            CubicPoly px, py;
            px.InitNonuniformCatmullRom(p0.X(), p1.X(), p2.X(), p3.X(), dt0, dt1, dt2);
            py.InitNonuniformCatmullRom(p0.Y(), p1.Y(), p2.Y(), p3.Y(), dt0, dt1, dt2);

            double d_theta_num = new_angles[i] - _points_polar[i1].Y();
            double d_theta_den = _points_polar[i2].Y() - _points_polar[i1].Y();
            if (d_theta_den < 0) {
                d_theta_den += CST_2_PI;
            }
            if (d_theta_num < 0) {
                d_theta_num += CST_2_PI;
            }
            double t = d_theta_num / d_theta_den;
            double x = px.eval(t);
            double y = py.eval(t);

            newline.cart().push_back(Point2d(x, y));
        }
        std::swap(_points_cart, newline._points_cart);
        to_polar();
    }
};

inline std::ostream& operator<<(std::ostream& stream, const BoundaryLine& rhs)
{
    stream.precision(5);
    for (size_t i = 0; i < rhs.size(); i++) {
        stream.width(4);
        stream << i << std::scientific
                << " ; x = " << rhs.cart(i).X() << " y = " << rhs.cart(i).Y()
                << " ; r = " << rhs.polar(i).X() << " t = " << rhs.polar(i).Y()
                << std::endl;
    }
    return stream;
}

int main(int argc, char ** argv)
{
    if (argc < 2) {
        LOG_INFO("Usage : " << argv[0] << " NAME_OUTER [NAME_INNER]...");
        return 1;
    }

    std::string name_outer = argv[1];
    BoundaryLine bl_outer;

    meshit::SetLogLevel(INFO_LOG_LEVEL);
    LOG_INFO(" == MeshIt Palpo == ");
    try {
        LOG_INFO("Outer boundary name : " << name_outer);
        bl_outer.read_boundary_file(name_outer);
    }
    catch (const std::runtime_error& e) {
        LOG_ERROR(e.what());
        return 1;
    }
    BoundaryLine sub_outer = bl_outer;
    sub_outer.resample(100);

    // creates geometry structure
    int bc_num = 1;
    meshit::SplineGeometry2d geom;
    geom.AddLine(sub_outer.cart(), 1e99, false, ++bc_num);

    // add holes
    std::vector<BoundaryLine> sub_inner;
    int nb_inner = argc - 2;

    for (int i = 0; i < nb_inner; i++) {

        std::string name_inner = argv[nb_inner + i + 1];
        try {
            LOG_INFO("Inner boundary name : " << name_inner);
            sub_inner.push_back(BoundaryLine(name_inner));
        }
        catch (const std::runtime_error& e) {
            LOG_ERROR(e.what());
            return 1;
        }
        sub_inner[i].resample(50);
        sub_inner[i].reverse();
        geom.AddLine(sub_inner[i].cart(), 1e99, true, ++bc_num);
    }

    geom.FakeData();
    geom.SetGrading(0.2);

    meshit::Mesh mesh;
    meshit::MeshingParameters mp;
    mp.optsteps2d = 3;

    LOG_INFO("start meshing");
    mesh.BuildFromSpline2D(geom, mp);
    LOG_INFO("meshing done");

    //    mesh.FindOpenElements(0);
    //    mesh.FindOpenElements(1);
//        meshit::MeshQuality2d(mesh);
//        meshit::CheckSurfaceMesh(mesh);
//        meshit::CheckSurfaceMesh2(mesh);
    //    mesh.CheckConsistentBoundary();
    //    mesh.PrintMemInfo(std::cout);
    mesh.Export("comprehension.msh", "Gmsh2 Format");

    return 0;
}
