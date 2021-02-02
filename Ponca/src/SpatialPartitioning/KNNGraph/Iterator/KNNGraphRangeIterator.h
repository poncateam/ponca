#pragma once

namespace Ponca {

class KNNGraphRangeQuery;

class KNNGraphRangeIterator
{
protected:
    friend class KNNGraphRangeQuery;

public:
    KNNGraphRangeIterator();
    KNNGraphRangeIterator(KNNGraphRangeQuery* query);
    KNNGraphRangeIterator(KNNGraphRangeQuery* query, int index);

public:
    bool operator != (const KNNGraphRangeIterator& other) const;
    void operator ++ ();
    int  operator *  () const;

protected:
    KNNGraphRangeQuery* m_query;
    int m_index;
};

}   
