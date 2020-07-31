#pragma once

#include <PCA/SpacePartitioning/Query/RangeIndexQuery.h>
#include <PCA/SpacePartitioning/KNNGraph/Iterator/KNNGraphRangeIterator.h>

#include <vector>
#include <stack>

namespace pca {

class KNNGraph;

class KNNGraphRangeQuery : public RangeIndexQuery
{
protected:
    friend class KNNGraphRangeIterator;

public:
    KNNGraphRangeQuery();
    KNNGraphRangeQuery(const KNNGraph* graph);
    KNNGraphRangeQuery(const KNNGraph* graph, Scalar radius);
    KNNGraphRangeQuery(const KNNGraph* graph, Scalar radius, int index);

public:
    KNNGraphRangeIterator begin();
    KNNGraphRangeIterator end();

protected:
    void initialize(KNNGraphRangeIterator& iterator);
    void advance(KNNGraphRangeIterator& iterator);

protected:
    const KNNGraph*   m_graph;
    std::vector<bool> m_flag;
    std::stack<int>   m_stack;
};

} // namespace pca
