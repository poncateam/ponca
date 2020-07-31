#pragma once

namespace pca {

class GridRangeIndexQuery;

class GridRangeIndexIterator
{
protected:
    friend class GridRangeIndexQuery;

public:
    GridRangeIndexIterator();
    GridRangeIndexIterator(GridRangeIndexQuery* query);
    GridRangeIndexIterator(GridRangeIndexQuery* query, int index);

public:
    bool operator !=(const GridRangeIndexIterator& other) const;
    void operator ++();
    int  operator * () const;

protected:
    GridRangeIndexQuery* m_query;
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
