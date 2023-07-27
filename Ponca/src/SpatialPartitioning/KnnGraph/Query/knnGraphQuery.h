/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <vector>

namespace Ponca {

class KnnGraph;

class KnnGraphQuery
{
public:
    using iterator = std::vector<int>::const_iterator;

public:
    KnnGraphQuery();
    KnnGraphQuery(const KnnGraph* graph, int index);

    iterator begin() const;
    iterator end() const;

protected:
    const KnnGraph* m_graph;
    int m_index;
};

} // namespace Ponca
