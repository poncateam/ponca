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

template<typename Scalar, typename VectorContainer, typename NeighborsIndexRange>
bool check_range_neighbors(const VectorContainer& points, const std::vector<int>& sampling, int index, Scalar r, NeighborsIndexRange& neighbors)
{
	if (has_duplicate(neighbors))
		return false;

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
		return false;

	for (int idx : neighbors) {
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
			return false;
		Scalar dist = (points[idx].pos() - points[index].pos()).norm();
		if (r < dist)
			return false;
	}

	for (int idx : sampling) {
		if (idx == index) continue;

		Scalar dist = (points[idx].pos() - points[index].pos()).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && r < dist)
			return false;
		if (!is_neighbor && dist < r)
			return false;
	}
	return true;
}

template<typename Scalar, typename VectorContainer, typename VectorType, typename NeighborsIndexRange>
bool check_range_neighbors(const VectorContainer& points, const std::vector<int>& sampling, const VectorType& point, Scalar r, NeighborsIndexRange& neighbors)
{
	if (has_duplicate(neighbors))
		return false;

	for (int idx : neighbors) {
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
			return false;
		Scalar dist = (points[idx].pos() - point).norm();
		if (r < dist)
			return false;
	}

	for (int idx : sampling) {
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

template<typename Scalar, typename VectorContainer>
bool check_k_nearest_neighbors(const VectorContainer& points, int index, const int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) > k && int(neighbors.size()) != k)
		return false;

	if (has_duplicate(neighbors))
		return false;

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
		return false;
	

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx].pos() - points[index].pos()).norm());

	for (int idx = 0; idx<int(points.size()); ++idx)
	{
		if (idx == index) continue;

		Scalar dist = (points[idx].pos() - points[index].pos()).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
			return false;
		if (!is_neighbor && dist < max_dist)
			return false;
	}
	return true;
}

template<typename Scalar, typename VectorContainer>
bool check_k_nearest_neighbors(const VectorContainer& points, const std::vector<int>& sampling, int index, const int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) > k && int(neighbors.size()) != k)
		return false;

	if (has_duplicate(neighbors))
		return false;

	for (int idx : neighbors) {
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
			return false;
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
		return false;

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - points[index]).norm());

	for (int idx : sampling) {
		if (idx == index) continue;

		Scalar dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
			return false;
		if (!is_neighbor && dist < max_dist)
			return false;
	}
	return true;
}

template<typename Scalar, typename VectorContainer, typename VectorType>
bool check_k_nearest_neighbors(const VectorContainer& points, const std::vector<int>& sampling, const VectorType& point, const int k, const std::vector<int>& neighbors)
{
	if (int(sampling.size()) >= k && int(neighbors.size()) != k)
		return false;

	if (has_duplicate(neighbors))
		return false;

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
			return false;
	}

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - point).norm());

	for (int idx : sampling) {
		Scalar dist = (points[idx] - point).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
			return false;
		if (!is_neighbor && dist < max_dist)
			return false;
	}
	return true;
}


template<typename Scalar, typename VectorContainer, typename VectorType>
bool check_k_nearest_neighbors(const VectorContainer& points, const VectorType& point, int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) >= k && int(neighbors.size()) != k)
		return false;

	if (has_duplicate(neighbors))
		return false;

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx].pos() - point).norm());

	for (int idx = 0; idx<int(points.size()); ++idx)
	{
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


template<typename Scalar, typename VectorContainer>
bool check_nearest_neighbor(const VectorContainer& points, int index, int nearest)
{
    return check_k_nearest_neighbors<Scalar>(points, index, 1, { nearest });
}
template<typename Scalar, typename VectorContainer, typename VectorType>
bool check_nearest_neighbor(const VectorContainer& points, const VectorType& point, int nearest)
{
    return check_k_nearest_neighbors<Scalar>(points, point, 1, { nearest });
}

// template<typename Scalar, typename VectorType, typename VectorContainer>
// bool check_nearest_neighbor(const VectorContainer& points, const std::vector<int>& sampling, const VectorType& point, int nearest)
// {
//     return check_k_nearest_neighbors<Scalar, VectorContainer>(points, sampling, point, 1, { nearest });
// }
//
// template<typename Scalar, typename VectorType, typename VectorContainer>
// bool check_nearest_neighbor(const VectorContainer& points, const std::vector<int>& sampling, int index, int nearest)
// {
//     return check_k_nearest_neighbors<Scalar, VectorType, VectorContainer>(points, sampling, index, 1, { nearest });
// }


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

template<bool doIndexQuery, typename DataPoint, typename PointContainer, typename RegularQueryFunctor, typename CheckQueryFunctor>
void testQuery(
	PointContainer& points,
	RegularQueryFunctor callRegularQuery,
	CheckQueryFunctor checkQuery,
	const int reserve_number = -1
) {
	using VectorType      = typename DataPoint::VectorType;

#pragma omp parallel for
	for (int i = 0; i < points.size(); ++i) {
		auto queryInput = [&]{
			// Do index query input
			if constexpr (doIndexQuery)
				return i;
			// Do position query input
			else
				return VectorType(VectorType::Random()); // values between [-1:1]
		}();
		std::vector<int> resQuery;
		if (reserve_number > 0) {
			resQuery.reserve(reserve_number);
		}
		for (int j : callRegularQuery(queryInput))
			resQuery.push_back(j);

		VERIFY((checkQuery(queryInput, resQuery)));
	}
}

template<bool doIndexQuery, typename DataPoint, typename PointContainer, typename MutableQueryFunctor, typename RegularQueryFunctor, typename CheckQueryFunctor>
void testQuery(
	PointContainer& points,
	MutableQueryFunctor callMutableQuery,
	RegularQueryFunctor callRegularQuery,
	CheckQueryFunctor checkQuery,
	const int reserve_number = -1
) {
	using VectorType      = typename DataPoint::VectorType;

#pragma omp parallel for
	for (int i = 0; i < points.size(); ++i) {
		auto queryInput = [&]{
			// Do index query input
			if constexpr (doIndexQuery)
				return i;
			// Do position query input
			else
				return VectorType(VectorType::Random()); // values between [-1:1]
		}();
		std::vector<int> resRegularQuery, resMutableQuery;
		if (reserve_number > 0) {
			resRegularQuery.reserve(reserve_number);
			resMutableQuery.reserve(reserve_number);
		}
		for (int j : callRegularQuery(queryInput))
			resRegularQuery.push_back(j);
		for (int j : callMutableQuery(queryInput))
			resMutableQuery.push_back(j);
		VERIFY((checkQuery(queryInput, resRegularQuery)));
		VERIFY((checkQuery(queryInput, resMutableQuery)));
	}
}