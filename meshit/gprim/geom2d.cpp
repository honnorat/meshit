#define _USE_MATH_DEFINES 1

#include "geom2d.hpp"
#include "geom3d.hpp"

namespace meshit
{
    Point2d::Point2d(const Point3d& p3)
        : px{p3.X()}, py{p3.Y()} { }


    inline const Box2d& Box2d::operator+=(const Box2d& b)
    {
        pmin.X() = std::min(pmin.X(), b.pmin.X());
        pmin.Y() = std::min(pmin.Y(), b.pmin.Y());
        pmax.X() = std::max(pmax.X(), b.pmax.X());
        pmax.Y() = std::max(pmax.Y(), b.pmax.Y());
        return *this;
    }


    std::ostream& operator<<(std::ostream& s, const Point2d& p)
    {
        return s << "(" << p.px << ", " << p.py << ")";
    }

    std::ostream& operator<<(std::ostream& s, const Vec2d& v)
    {
        return s << "(" << v.vx << ", " << v.vy << ")";
    }

}  // namespace meshit
