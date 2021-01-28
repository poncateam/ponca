#pragma once

//#include "../query.h" 
#include "../iterator.h"

#include "../kdTree.h"

#include "../../query.h"
//#include "../../defines.h"

//#include <PCA/SpacePartitioning/Query/NearestPointQuery.h>
//#include <PCA/SpacePartitioning/KdTree/Query/KdTreeQuery.h>
//#include <PCA/SpacePartitioning/KdTree/Iterator/KdTreeNearestPointIterator.h>

//class Ponca::KdTree;
//struct Ponca::KdTreeQuery;
//template <typename _VectorType>
//struct Ponca::NearestPointQuery;
//struct Ponca::KdTreeNearestPointIterator;
namespace pca {

template <typename _VectorType>
class KdTreeNearestPointQuery : public NearestPointQuery<_VectorType>, public KdTreeQuery
{
    using VectorType = typename NearestPointQuery<_VectorType>::VectorType;
public:
    KdTreeNearestPointQuery();
    KdTreeNearestPointQuery(const KdTree* kdtree);
    KdTreeNearestPointQuery(const KdTree* kdtree, const VectorType& point);

public:
    //KdTreeNearestPointIterator begin();
    //KdTreeNearestPointIterator end();

protected:
    void search();
};

} // namespace pca
//#include "./KdTreeNearestPointQuery.cpp"
