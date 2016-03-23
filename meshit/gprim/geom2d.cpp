#define _USE_MATH_DEFINES 1

#include "geom2d.hpp"

namespace meshit
{
    std::ostream& operator<<(std::ostream& s, const Point2d& p)
    {
        return s << "(" << p.px << ", " << p.py << ")";
    }

    std::ostream& operator<<(std::ostream& s, const Vec2d& v)
    {
        return s << "(" << v.vx << ", " << v.vy << ")";
    }

}  // namespace meshit
