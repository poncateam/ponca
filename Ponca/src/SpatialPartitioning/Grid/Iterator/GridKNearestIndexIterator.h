#pragma once

#include <PCA/SpacePartitioning/Query/KNearestQuery.h>

namespace Ponca {

class GridKNearestIndexIterator
{
public:
    GridKNearestIndexIterator();
    GridKNearestIndexIterator(limited_priority_queue<IndexSquaredDistance>::iterator iterator);

public:
    bool operator !=(const GridKNearestIndexIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    limited_priority_queue<IndexSquaredDistance>::iterator m_iterator;
};

}   
