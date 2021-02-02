#pragma once

//#include "../query.h" 
#include "../iterator.h"

#include "../kdTree.h"

#include "../../query.h"
//#include "../../defines.h"

//#include <Ponca/src/SpatialPartitioning/Query/NearestPointQuery.h>
//#include <Ponca/src/SpatialPartitioning/KdTree/Query/KdTreeQuery.h>
//#include <PCA/SpacePartitioning/KdTree/Iterator/KdTreeNearestPointIterator.h>

//class Ponca::KdTree;
//template <typename _VectorType>
//struct Ponca::NearestPointQuery;
//struct Ponca::KdTreeNearestPointIterator;
namespace Ponca {

struct KdTreeQuery;
template <typename _VectorType>
class KdTreeNearestPointQuery : public NearestPointQuery<_VectorType>, public KdTreeQuery
{
using VectorType = typename NearestPointQuery<_VectorType>::VectorType;
public:
    KdTreeNearestPointQuery();
    KdTreeNearestPointQuery(const KdTree* kdtree);
    KdTreeNearestPointQuery(const KdTree* kdtree, const _VectorType& point);

public:
    //KdTreeNearestPointIterator begin();
    //KdTreeNearestPointIterator end();

protected:
    void search();
};

}   
#include "./KdTreeNearestPointQuery.hpp"
