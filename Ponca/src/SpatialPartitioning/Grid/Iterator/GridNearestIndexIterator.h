#pragma once

namespace pca {

class GridNearestIndexIterator
{
public:
    GridNearestIndexIterator();
    GridNearestIndexIterator(int index);

public:
    bool operator !=(const GridNearestIndexIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    int m_index;
};

} // namespace pca
