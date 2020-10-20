#pragma once

namespace pca {

class GridNearestPointIterator
{
public:
    GridNearestPointIterator();
    GridNearestPointIterator(int index);

public:
    bool operator !=(const GridNearestPointIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    int m_index;
};

} // namespace pca
