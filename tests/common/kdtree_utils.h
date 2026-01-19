#include <chrono>
#include <Eigen/Dense>
#include <cassert>

#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>
#include <Ponca/src/SpatialPartitioning/KnnGraph/knnGraph.h>
#include <chrono>
#include <random>
#include <iterator>
#include <algorithm>

#include "../common/has_duplicate.h"

#define EXPECT_EQ(a,b) assert(a==b)
#define EXPECT_TRUE(a) assert(a==true)


using namespace Eigen;
using namespace std;

using uint = unsigned int;

using Vec3 = Eigen::Matrix<float, 3, 1>;
using Vector3Array = std::vector<Vec3>;


template<typename type, int size>
class TestPoint {
public:
	enum { Dim = size };
	typedef type Scalar;
	typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;

	PONCA_MULTIARCH inline TestPoint(const VectorType& pos = VectorType::Zero())
		: _pos(pos) {}
	PONCA_MULTIARCH inline const VectorType& pos() const { return _pos; }
	PONCA_MULTIARCH inline       VectorType& pos() { return _pos; }
private:
	VectorType _pos;
};

class MyPoint {
public:
	enum { Dim = 3 };
	typedef float Scalar;
	typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;

	PONCA_MULTIARCH inline MyPoint(const VectorType& pos = VectorType::Zero())
		: _pos(pos) {}
	PONCA_MULTIARCH inline const VectorType& pos() const { return _pos; }
	PONCA_MULTIARCH inline       VectorType& pos() { return _pos; }
private:
	VectorType _pos;
};

template<typename DataPoint, typename VectorContainer, typename QueryInput, typename NeighborsIndexRange>
bool checkRangeNeighbors(const VectorContainer& points, const std::vector<int>& sampling, QueryInput& queryInput, typename DataPoint::Scalar r, NeighborsIndexRange& neighbors) {
	using Scalar     = typename DataPoint::Scalar;
	using VectorType = typename DataPoint::VectorType;

	if (hasDuplicate(neighbors))
		return false;

	VectorType point;
	constexpr bool isIndexQuery {std::is_same_v<QueryInput, int>}; // Determines if the query input is the eval point or it's index
	if constexpr (isIndexQuery) { // queryInput == index of the eval point
		auto it = std::find(neighbors.begin(), neighbors.end(), queryInput);
		if (it != neighbors.end()) return false;

		point = points[queryInput].pos();
	} else { // queryInput == eval point
		point = queryInput;
	}

	for (int idx : neighbors) {
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
			return false;
		Scalar dist = (points[idx].pos() - point).norm();
		if (r < dist)
			return false;
	}

	for (int idx : sampling) {
		if constexpr (isIndexQuery)
			if (idx == queryInput) continue;

		Scalar dist = (points[idx].pos() - point).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && r < dist)
			return false;
		if (!is_neighbor && dist < r)
			return false;
	}
	return true;
}

template<typename DataPoint, typename VectorContainer, typename QueryInput>
bool checkKNearestNeighbors(const VectorContainer& points, QueryInput& queryInput, const int k, const std::vector<int>& neighbors)
{
	using Scalar     = typename DataPoint::Scalar;
	using VectorType = typename DataPoint::VectorType;

	VectorType point;
	constexpr bool isIndexQuery {std::is_same_v<QueryInput, int>}; // Determines if the query input is the eval point or it's index
	if constexpr (isIndexQuery) { // queryInput == index of the eval point
		auto it = std::find(neighbors.begin(), neighbors.end(), queryInput);
		if (it != neighbors.end()) return false;

		point = points[queryInput].pos();
	} else { // queryInput == eval point
		point = queryInput;
	}

	if (int(points.size()) > k && int(neighbors.size()) != k)
		return false;

	if (hasDuplicate(neighbors))
		return false;

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx].pos() - point).norm());

	for (int idx = 0; idx<int(points.size()); ++idx)
	{
		if constexpr (isIndexQuery)
			if (idx == queryInput) continue;

		Scalar dist = (points[idx].pos() - point).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
			return false;
		if (!is_neighbor && dist < max_dist)
			return false;
	}
	return true;
}

template<typename DataPoint, typename VectorContainer, typename QueryInput>
bool checkKNearestNeighbors(const VectorContainer& points, const std::vector<int>& sampling, QueryInput& queryInput, const int k, const std::vector<int>& neighbors)
{
	using Scalar     = typename DataPoint::Scalar;
	using VectorType = typename DataPoint::VectorType;

	if (int(points.size()) > k && int(neighbors.size()) != k)
		return false;

	VectorType point;
	constexpr bool isIndexQuery {std::is_same_v<QueryInput, int>}; // Determines if the query input is the eval point or it's index
	if constexpr (isIndexQuery) { // queryInput == index of the eval point
		auto it = std::find(neighbors.begin(), neighbors.end(), queryInput);
		if (it != neighbors.end()) return false;

		point = points[queryInput].pos();
	} else { // queryInput == eval point
		point = queryInput;
	}

	if (hasDuplicate(neighbors))
		return false;

	for (int idx : neighbors) {
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
			return false;
	}

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx].pos() - point).norm());

	for (int idx : sampling) {
		if constexpr (isIndexQuery)
			if (idx == queryInput) continue;

		Scalar dist = (points[idx].pos() - point).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
			return false;
		if (!is_neighbor && dist < max_dist)
			return false;
	}
	return true;
}

template<typename DataPoint, typename VectorContainer, typename QueryInput>
bool checkNearestNeighbor(const VectorContainer& points, QueryInput& queryInput, int nearest)
{
    return checkKNearestNeighbors<DataPoint>(points, queryInput, 1, { nearest });
}

template<typename DataPoint, typename VectorContainer, typename QueryInput>
bool checkNearestNeighbor(const VectorContainer& points, const std::vector<int>& sampling, QueryInput& queryInput, int nearest)
{
    return checkKNearestNeighbors<DataPoint, VectorContainer>(points, sampling, queryInput, 1, { nearest });
}

template<typename DataPoint>
void generateData(std::vector<DataPoint>& points) {
	using VectorType      = typename DataPoint::VectorType;
	std::generate(points.begin(), points.end(), []() { return DataPoint(VectorType::Random()); });
}

// For subsampling
template<typename DataPoint>
std::unique_ptr<Ponca::KdTreeSparse<DataPoint>> buildSubsampledKdTree(std::vector<DataPoint>& points, std::vector<int>& sampling) {
	std::vector<int> indices(points.size());
	std::iota(indices.begin(), indices.end(), 0);
	sampling.resize(points.size() / 2);

	int seed = 0;
	std::sample(indices.begin(), indices.end(), sampling.begin(), points.size() / 2, std::mt19937(seed));

	/// [KdTree assign sparse]
	return std::make_unique<Ponca::KdTreeSparse<DataPoint>>(points, sampling);
	/// [KdTree assign sparse]
}

template<typename DataPoint>
std::unique_ptr<Ponca::KdTreeDense<DataPoint>> buildKdTreeDense(std::vector<DataPoint>& points, std::vector<int>& sampling) {
	sampling.resize(points.size());
	std::iota(sampling.begin(), sampling.end(), 0);

	/// [KdTree assign dense]
	return std::make_unique<Ponca::KdTreeDense<DataPoint>>(points);
	/// [KdTree assign dense]
}


template< bool doIndexQuery,
	typename DataPoint, typename PointContainer,
	typename QueryFunctor, typename CheckQueryFunctor,
	typename... QueryInputTypes>
std::chrono::milliseconds
testQuery(
	PointContainer &points,
	QueryFunctor callQuery,
	CheckQueryFunctor checkQuery,
	const int retry_number = 1,
	QueryInputTypes &&... outs
) {
	using VectorType = typename DataPoint::VectorType;
	const auto time = std::chrono::system_clock::now();

#ifdef NDEBUG
#pragma omp parallel for
#endif
	for (int i = 0; i < points.size(); ++i) {
		auto queryInput = [i]{
			// Call query with index input
			if constexpr (doIndexQuery)
				return i;
			// Call query with position input
			else
				return VectorType(VectorType::Random()); // values between [-1:1]
		}(); // Either an index or a position

		for (int j = 0; j < retry_number; ++j) {
			std::vector<int> resQuery;

			for (int point_idx : callQuery(queryInput, std::forward<QueryInputTypes>(outs)...))
				resQuery.push_back(point_idx);

			VERIFY((checkQuery(queryInput, resQuery)));
		}
	}
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - time);
}

template<bool doIndexQuery,
	typename DataPoint, typename PointContainer,
	typename MutableQueryFunctor, typename RegularQueryFunctor, typename CheckQueryFunctor,
	typename... QueryInputTypes>
std::chrono::milliseconds
testQuery(
	PointContainer& points,
	MutableQueryFunctor callMutableQuery,
	RegularQueryFunctor callRegularQuery,
	CheckQueryFunctor checkQuery,
	const int retry_number = 1,
	QueryInputTypes&&... outs
) {
	using VectorType = typename DataPoint::VectorType;
	const auto time = std::chrono::system_clock::now();

#ifdef NDEBUG
#pragma omp parallel for
#endif
	for (int i = 0; i < points.size(); ++i) {

		auto queryInput = [i]{
			// Do index query input
			if constexpr (doIndexQuery)
				return i;
			// Do position query input
			else
				return VectorType(VectorType::Random()); // values between [-1:1]
		}();
		auto mutableQuery = callMutableQuery();

		for (int j = 0; j < retry_number; ++j) {
			std::vector<int> resRegularQuery, resMutableQuery;

			for (int point_idx : callRegularQuery(queryInput, std::forward<QueryInputTypes>(outs)...))
				resRegularQuery.push_back(point_idx);
			for (int point_idx : mutableQuery(queryInput, std::forward<QueryInputTypes>(outs)...))
				resMutableQuery.push_back(point_idx);

			VERIFY((checkQuery(queryInput, resRegularQuery)));
			VERIFY((checkQuery(queryInput, resMutableQuery)));
		}
	}
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - time);
}
