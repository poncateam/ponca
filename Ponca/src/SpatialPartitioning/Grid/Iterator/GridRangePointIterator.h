#pragma once

namespace Ponca {

class GridRangePointQuery;

class GridRangePointIterator
{
protected:
    friend class GridRangePointQuery;

public:
    GridRangePointIterator();
    GridRangePointIterator(GridRangePointQuery* query);
    GridRangePointIterator(GridRangePointQuery* query, int index);

public:
    bool operator !=(const GridRangePointIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    GridRangePointQuery* m_query;
    int m_index;

    int m_i_start;
    int m_j_start;
    int m_k_start;
    int m_i_end;
    int m_j_end;
    int m_k_end;

    int m_i;
    int m_j;
    int m_k;

    int m_it;
    int m_idx_cell;
};

} // namespace pca
