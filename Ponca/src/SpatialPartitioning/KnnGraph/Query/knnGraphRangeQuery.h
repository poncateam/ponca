/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../query.h"
#include "../Iterator/KnnGraphRangeIterator.h"

#include <vector>
#include <stack>

namespace Ponca {

class KnnGraph;

class KnnGraphRangeQuery : public RangeIndexQuery
{
protected:
    friend class KnnGraphRangeIterator;

public:
    KnnGraphRangeQuery();
    KnnGraphRangeQuery(const KnnGraph* graph);
    KnnGraphRangeQuery(const KnnGraph* graph, Scalar radius);
    KnnGraphRangeQuery(const KnnGraph* graph, Scalar radius, int index);

public:
    KnnGraphRangeIterator begin();
    KnnGraphRangeIterator end();

protected:
    void initialize(KnnGraphRangeIterator& iterator);
    void advance(KnnGraphRangeIterator& iterator);

protected:
    const KnnGraph*   m_graph;
    std::vector<bool> m_flag;
    std::stack<int>   m_stack;
};

} // namespace Ponca
