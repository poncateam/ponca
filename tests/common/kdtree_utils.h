#include <chrono>
#include <Eigen/Dense>
#include <cassert>

#include <chrono>
#include <random>
#include <iterator>
#include <algorithm>


#define EXPECT_EQ(a,b) assert(a==b)
#define EXPECT_TRUE(a) assert(a==true)


using namespace Eigen;
using namespace std;

using uint = unsigned int;

using Vec3 = Eigen::Matrix<float, 3, 1>;
using Vector3Array = std::vector<Vec3>;


template <typename Scalar, int size>
class _Point {
public:
	typedef Scalar Scalar;
	typedef Eigen::Matrix<Scalar, size, 1> VectorType;
	typedef std::vector<VectorType> Vector;
	typedef Eigen::AlignedBox<Scalar, size> Aabb;
};

template<typename Scalar, typename Vector>
bool check_range_neighbors(const Vector& points, const std::vector<int>& sampling, int index, Scalar r, const std::vector<int>& neighbors)
{
	if (has_duplicate(neighbors))
	{
		return false;
	}

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
		{
			return false;
		}
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
	{
		return false;
	}

	for (int i = 0; i<int(neighbors.size()); ++i)
	{
		Scalar dist = (points[neighbors[i]] - points[index]).norm();
		if (r < dist)
		{
			return false;
		}
	}

	for (int i = 0; i<int(sampling.size()); ++i)
	{
		int idx = sampling[i];
		if (idx == index) continue;

		Scalar dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && r < dist)
		{
			return false;
		}
		if (!is_neighbor && dist < r)
		{
			return false;
		}
	}
	return true;
}

template<typename Scalar, typename Vector>
bool check_k_nearest_neighbors(const Vector& points, int index, int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) > k && int(neighbors.size()) != k)
	{
		return false;
	}

	if (has_duplicate(neighbors))
	{
		return false;
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
	{
		return false;
	}

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - points[index]).norm());

	for (int idx = 0; idx<int(points.size()); ++idx)
	{
		if (idx == index) continue;

		Scalar dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
		{
			return false;
		}
		if (!is_neighbor && dist < max_dist)
		{
			return false;
		}
	}
	return true;
}

template<typename Scalar, typename Vector>
bool check_nearest_neighbors(const Vector& points, int index, int nearest)
{
	return check_k_nearest_neighbors<Scalar, Vector>(points, index, 1, { nearest });
}

template<typename Scalar, typename Vector>
bool check_k_nearest_neighbors(const Vector& points, const std::vector<int>& sampling, int index, int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) > k && int(neighbors.size()) != k)
	{
		return false;
	}

	if (has_duplicate(neighbors))
	{
		return false;
	}

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
		{
			return false;
		}
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
	{
		return false;
	}

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - points[index]).norm());

	for (int i = 0; i<int(sampling.size()); ++i)
	{
		int idx = sampling[i];
		if (idx == index) continue;

		Scalar dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
		{
			return false;
		}
		if (!is_neighbor && dist < max_dist)
		{
			return false;
		}
	}
	return true;
}

template<typename Scalar, typename VectorType, typename Vector>
bool check_k_nearest_neighbors(const Vector& points, const std::vector<int>& sampling, const VectorType& point, int k, const std::vector<int>& neighbors)
{
	if (int(sampling.size()) >= k && int(neighbors.size()) != k)
	{
		return false;
	}

	if (has_duplicate(neighbors))
	{
		return false;
	}

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
		{
			return false;
		}
	}

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - point).norm());

	for (int i = 0; i<int(sampling.size()); ++i)
	{
		int idx = sampling[i];
		Scalar dist = (points[idx] - point).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
		{
			return false;
		}
		if (!is_neighbor && dist < max_dist)
		{
			return false;
		}
	}
	return true;
}


template<typename Scalar, typename VectorType, typename Vector>
bool check_nearest_neighbors(const Vector& points, const std::vector<int>& sampling, const VectorType& point, int nearest)
{
	return check_k_nearest_neighbors<Scalar, Vector>(points, sampling, point, 1, { nearest });
}


template<typename Scalar, typename VectorType, typename Vector>
bool check_nearest_neighbors(const Vector& points, const std::vector<int>& sampling, int index, int nearest)
{
	return check_k_nearest_neighbors<Scalar, VectorType, Vector>(points, sampling, index, 1, { nearest });
}

template<typename Scalar, typename VectorType, typename Vector>
bool check_range_neighbors(const Vector& points, const std::vector<int>& sampling, const VectorType& point, Scalar r, const std::vector<int>& neighbors)
{
	if (has_duplicate(neighbors))
	{
		return false;
	}

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
		{
			return false;
		}
	}

	for (int i = 0; i<int(neighbors.size()); ++i)
	{
		Scalar dist = (points[neighbors[i]] - point).norm();
		if (r < dist)
		{
			return false;
		}
	}

	for (int i = 0; i<int(sampling.size()); ++i)
	{
		int idx = sampling[i];
		Scalar dist = (points[idx] - point).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && r < dist)
		{
			return false;
		}
		if (!is_neighbor && dist < r)
		{
			return false;
		}
	}
	return true;
}
