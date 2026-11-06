// SequenceAnalysis, Copyright Switch Lab 2026
// Joost Schymkowitz

#include <unordered_map>
#include <vector>
#include <cmath>
#include "vector3D.h"

struct SpatialGrid
{
    double cell;
    std::unordered_map<long long, std::vector<Vector3D>> bins;

    explicit SpatialGrid(double cellSize) : cell(cellSize) {}

    static long long Key(int ix, int iy, int iz)
    {
        return ((long long)ix << 42) ^ ((long long)iy << 21) ^ (long long)iz;
    }

    void Insert(const Vector3D& p)
    {
        int ix = (int)std::floor(p.x / cell);
        int iy = (int)std::floor(p.y / cell);
        int iz = (int)std::floor(p.z / cell);
        bins[Key(ix, iy, iz)].push_back(p);
    }

    void InsertAll(const std::vector<Vector3D>& pts)
    {
        for(const auto& p : pts)
            Insert(p);
    }

    // Count the number of inserted points within `dist` of `p`.
    // Note: if `p` itself was inserted its contribution is included; subtract 1 if needed.
    int CountWithin(const Vector3D& p, double dist) const
    {
        const double dist2 = dist * dist;
        int count = 0;
        int ix = (int)std::floor(p.x / cell);
        int iy = (int)std::floor(p.y / cell);
        int iz = (int)std::floor(p.z / cell);
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++)
                {
                    auto it = bins.find(Key(ix+dx, iy+dy, iz+dz));
                    if (it == bins.end()) continue;
                    for (const auto& q : it->second)
                    {
                        double ddx = p.x - q.x, ddy = p.y - q.y, ddz = p.z - q.z;
                        if (ddx*ddx + ddy*ddy + ddz*ddz <= dist2) count++;
                    }
                }
        return count;
    }

    bool AnyWithin(const Vector3D& p, double dist) const
    {
        const double dist2 = dist * dist;

        int ix = (int)std::floor(p.x / cell);
        int iy = (int)std::floor(p.y / cell);
        int iz = (int)std::floor(p.z / cell);

        for(int dx = -1; dx <= 1; dx++)
        {
            for(int dy = -1; dy <= 1; dy++)
            {
                for(int dz = -1; dz <= 1; dz++)
                {
                    auto it = bins.find(Key(ix + dx, iy + dy, iz + dz));
                    if(it == bins.end()) continue;

                    for(const auto& q : it->second)
                    {
                        double ddx = p.x - q.x;
                        double ddy = p.y - q.y;
                        double ddz = p.z - q.z;
                        double r2 = ddx*ddx + ddy*ddy + ddz*ddz;

                        if(r2 <= dist2)
                            return true;
                    }
                }
            }
        }

        return false;
    }
};

