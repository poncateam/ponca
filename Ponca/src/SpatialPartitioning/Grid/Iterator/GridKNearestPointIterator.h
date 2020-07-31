#pragma once

#include <PCA/SpacePartitioning/Query/KNearestQuery.h>

namespace pca {

class GridKNearestPointIterator
{
public:
    GridKNearestPointIterator();
    GridKNearestPointIterator(limited_priority_queue<IndexSquaredDistance>::iterator iterator);

public:
    bool operator !=(const GridKNearestPointIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    limited_priority_queue<IndexSquaredDistance>::iterator m_iterator;
};

} // namespace pca
