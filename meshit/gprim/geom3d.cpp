#include "geom3d.hpp"

namespace meshit
{
    std::ostream& operator<<(std::ostream& s, const Point3d& p)
    {
        return s << "(" << p.x[0] << ", " << p.x[1] << ", " << p.x[2] << ")";
    }

    std::ostream& operator<<(std::ostream& s, const Vec3d& v)
    {
        return s << "(" << v.x[0] << ", " << v.x[1] << ", " << v.x[2] << ")";
    }

    double Angle(const Vec3d& v1, const Vec3d& v2)
    {
        double co = (v1 * v2) / (v1.Length() * v2.Length());
        if (co > 1) co = 1;
        if (co < -1) co = -1;
        return acos(co);
    }

    void Vec3d::GetNormal(Vec3d& n) const
    {
        if (fabs(X()) > fabs(Z())) {
            n.X() = -Y();
            n.Y() = X();
            n.Z() = 0;
        } else {
            n.X() = 0;
            n.Y() = Z();
            n.Z() = -Y();
        }
        double len = n.Length();
        if (len == 0) {
            n.X() = 1;
            n.Y() = n.Z() = 0;
        }
        else
            n /= len;
    }

    Box3d::Box3d(
        double aminx, double amaxx,
        double aminy, double amaxy,
        double aminz, double amaxz)
    {
        minx[0] = aminx;
        maxx[0] = amaxx;
        minx[1] = aminy;
        maxx[1] = amaxy;
        minx[2] = aminz;
        maxx[2] = amaxz;
    }

    Box3d::Box3d(const Box3d& b2)
    {
        for (int i = 0; i < 3; i++) {
            minx[i] = b2.minx[i];
            maxx[i] = b2.maxx[i];
        }
    }

    Box3d::Box3d(const Box<3>& b2)
    {
        for (int i = 0; i < 3; i++) {
            minx[i] = b2.PMin()[i];
            maxx[i] = b2.PMax()[i];
        }
    }

    Box3d::Box3d(const Point3d& p1, const Point3d& p2)
    {
        minx[0] = std::min(p1.X(), p2.X());
        minx[1] = std::min(p1.Y(), p2.Y());
        minx[2] = std::min(p1.Z(), p2.Z());
        maxx[0] = std::max(p1.X(), p2.X());
        maxx[1] = std::max(p1.Y(), p2.Y());
        maxx[2] = std::max(p1.Z(), p2.Z());
    }

    const Box3d& Box3d::operator+=(const Box3d& b)
    {
        minx[0] = std::min(minx[0], b.minx[0]);
        minx[1] = std::min(minx[1], b.minx[1]);
        minx[2] = std::min(minx[2], b.minx[2]);
        maxx[0] = std::max(maxx[0], b.maxx[0]);
        maxx[1] = std::max(maxx[1], b.maxx[1]);
        maxx[2] = std::max(maxx[2], b.maxx[2]);

        return *this;
    }

    void Transpose(Vec3d& v1, Vec3d& v2, Vec3d& v3)
    {
        std::swap(v1.Y(), v2.X());
        std::swap(v1.Z(), v3.X());
        std::swap(v2.Z(), v3.Y());
    }

    int SolveLinearSystem(
        const Vec3d& col1, const Vec3d& col2, const Vec3d& col3,
        const Vec3d& rhs, Vec3d& sol)
    {
        // changed by MW
        double matrix[3][3];
        double locrhs[3];
        int retval = 0;

        for (int i = 0; i < 3; i++) {
            matrix[i][0] = col1.X(i + 1);
            matrix[i][1] = col2.X(i + 1);
            matrix[i][2] = col3.X(i + 1);
            locrhs[i] = rhs.X(i + 1);
        }

        for (int i = 0; i < 2; i++) {
            int pivot = i;
            double maxv = fabs(matrix[i][i]);
            for (int j = i + 1; j < 3; j++)
                if (fabs(matrix[j][i]) > maxv) {
                    maxv = fabs(matrix[j][i]);
                    pivot = j;
                }

            if (fabs(maxv) > 1e-40) {
                if (pivot != i) {
                    std::swap(matrix[i][0], matrix[pivot][0]);
                    std::swap(matrix[i][1], matrix[pivot][1]);
                    std::swap(matrix[i][2], matrix[pivot][2]);
                    std::swap(locrhs[i], locrhs[pivot]);
                }
                for (int j = i + 1; j < 3; j++) {
                    double fac = matrix[j][i] / matrix[i][i];

                    for (int k = i + 1; k < 3; k++)
                        matrix[j][k] -= fac * matrix[i][k];
                    locrhs[j] -= fac * locrhs[i];
                }
            }
            else
                retval = 1;
        }

        if (fabs(matrix[2][2]) < 1e-40)
            retval = 1;

        if (retval != 0)
            return retval;


        for (int i = 2; i >= 0; i--) {
            double sum = locrhs[i];
            for (int j = 2; j > i; j--)
                sum -= matrix[i][j] * sol.X(j + 1);

            sol.X(i + 1) = sum / matrix[i][i];
        }

        return 0;
    }

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

}  // namespace meshit
