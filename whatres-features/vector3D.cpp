#include "vector3D.h"
#include <cmath>

bool Vector3D::IsFinite3(const Vector3D& p)
{
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

Vector3D Vector3D::Mean(const std::vector<Vector3D>& pts)
{
    double x=0, y=0, z=0;
    if(pts.empty()) return Vector3D(0,0,0);

    for(const auto& p : pts)
    {
        x += p.x;
        y += p.y;
        z += p.z;
    }

    return Vector3D(x/pts.size(), y/pts.size(), z/pts.size());
}

double Vector3D::Dot(const Vector3D& a, const Vector3D& b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vector3D Vector3D::PerpProject(const Vector3D& p, const Vector3D& axisUnit)
{
    double t = Dot(p, axisUnit);
    return Vector3D(p.x - t*axisUnit.x,
                    p.y - t*axisUnit.y,
                    p.z - t*axisUnit.z);
}