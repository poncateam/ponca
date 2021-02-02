#pragma once

namespace Ponca {

struct Index3D
{
    int i;
    int j;
    int k;

    inline bool operator == (const Index3D& other) const
    {
        return i == other.i && j == other.j && k == other.k;
    }
};

}   
