#pragma once

#include <vector>

namespace pca {

class KNNGraph;

class KNNGraphQuery
{
public:
    using iterator = std::vector<int>::const_iterator;

public:
    KNNGraphQuery();
    KNNGraphQuery(const KNNGraph* graph, int index);

    iterator begin() const;
    iterator end() const;

protected:
    const KNNGraph* m_graph;
    int m_index;
};

} // namespace pca
