#ifndef VECTOR3D_H
#define VECTOR3D_H

#include "misc_math.h"
#include <vector>
#include <iostream>
#include <cmath>

#ifndef EPS
#define EPS 1.0E-7
#endif

#define COORD_TOLERANCE 0.005

class Vector3D
{
public:
    double x, y, z;

    Vector3D() : x(0), y(0), z(0) { }
    Vector3D(double xx, double yy, double zz) : x(xx), y(yy), z(zz) { }
    Vector3D(const Vector3D& rhs) : x(rhs.x), y(rhs.y), z(rhs.z) { }

    Vector3D(const std::vector<double>& v) : x(v[0]), y(v[1]), z(v[2]) { }

    operator std::vector<double>()
    {
        std::vector<double> ret;
        ret.reserve(3);
        ret.push_back(x);
        ret.push_back(y);
        ret.push_back(z);
        return ret;
    }

    static bool IsFinite3(const Vector3D& p);
    static Vector3D Mean(const std::vector<Vector3D>& pts);
    static double Dot(const Vector3D& a, const Vector3D& b);
    static Vector3D PerpProject(const Vector3D& p, const Vector3D& axisUnit);

    double Get(unsigned int i) const
    {
        if(i == 0) return x;
        return ((i < 2) ? y : z);
    }

    Vector3D operator+=(const Vector3D& rhs)
    {
        x += rhs.x; y += rhs.y; z += rhs.z;
        return *this;
    }

    Vector3D operator-=(const Vector3D& rhs)
    {
        x -= rhs.x; y -= rhs.y; z -= rhs.z;
        return *this;
    }

    Vector3D operator*=(double d)
    {
        x *= d; y *= d; z *= d;
        return *this;
    }

    Vector3D operator/=(double d)
    {
        x /= d; y /= d; z /= d;
        return *this;
    }

    void Clear()
    {
        x = 0.0; y = 0.0; z = 0.0;
    }

    Vector3D& operator=(const Vector3D& rhs)
    {
        x = rhs.x; y = rhs.y; z = rhs.z;
        return *this;
    }

    Vector3D& operator=(const std::vector<double>& rhs)
    {
        x = rhs[0]; y = rhs[1]; z = rhs[2];
        return *this;
    }

    Vector3D Normalize()
    {
        double l = std::sqrt(x*x + y*y + z*z);
        if(l > EPS)
        {
            x /= l; y /= l; z /= l;
        }
        return *this;
    }

    void print_xyz()
    {
        std::cout << x << " " << y << " " << z << std::endl;
    }
};

inline bool IsEqual(const Vector3D& v1, const Vector3D& v2)
{
    if(std::fabs(v1.x-v2.x) > COORD_TOLERANCE ||
       std::fabs(v2.y-v1.y) > COORD_TOLERANCE ||
       std::fabs(v2.z-v1.z) > COORD_TOLERANCE)
        return false;
    return true;
}

inline bool operator==(const Vector3D& v1, const Vector3D& v2)
{
    if(v1.x==v2.x && v1.y==v2.y && v2.z==v2.z) return true;
    return false;
}

inline std::ostream& operator<<(std::ostream& o, const Vector3D& rhs)
{
    o << "x=" << rhs.x << " y=" << rhs.y << " z=" << rhs.z << std::endl;
    return o;
}

template<class Type>
inline bool operator==(const std::vector<Type>& lhs, const std::vector<Type>& rhs)
{
    if(lhs.size() != rhs.size()) return false;

    typename std::vector<Type>::const_iterator it1 = lhs.begin();
    typename std::vector<Type>::const_iterator it2 = rhs.begin();
    typename std::vector<Type>::const_iterator itend = lhs.end();

    for(; it1 != itend; ++it1, ++it2)
        if(!(*it1 == *it2)) return false;

    return true;
}

inline Vector3D operator*(double d, const Vector3D& r)
{
    return Vector3D(d*r.x, d*r.y, d*r.z);
}

inline Vector3D operator+(const Vector3D& r1, const Vector3D& r2)
{
    return Vector3D(r1.x+r2.x, r1.y+r2.y, r1.z+r2.z);
}

inline Vector3D operator-(const Vector3D& r1, const Vector3D& r2)
{
    return Vector3D(r1.x-r2.x, r1.y-r2.y, r1.z-r2.z);
}

inline Vector3D operator-(const Vector3D& r)
{
    return Vector3D(-r.x, -r.y, -r.z);
}

inline Vector3D operator%(const Vector3D& r1, const Vector3D& r2)
{
    return Vector3D(r1.y*r2.z-r1.z*r2.y,
                    r1.z*r2.x-r1.x*r2.z,
                    r1.x*r2.y-r1.y*r2.x);
}

inline double operator*(const Vector3D& r1, const Vector3D& r2)
{
    return r1.x*r2.x + r1.y*r2.y + r1.z*r2.z;
}

inline Vector3D operator/(const Vector3D& r1, const double d)
{
    return Vector3D(r1.x/d, r1.y/d, r1.z/d);
}

inline double Dihedral(const Vector3D &i,
                       const Vector3D &j,
                       const Vector3D &k,
                       const Vector3D &l)
{
    double dih;

    Vector3D jk = k - j;

    Vector3D c = (i - j) % jk;
    c.Normalize();

    Vector3D d = (l - k) % jk;
    d.Normalize();

    double scal = c * d;

    if (fabs(scal + 1.0) < EPS)
        dih = 180.0;
    else if (fabs(scal - 1.0) < EPS)
        dih = 0.0;
    else
        dih = 57.2958 * acos(scal);

    double chiral = (c % d) * jk;
    if (chiral < 0)
        dih = -dih;

    return dih;
}

#endif