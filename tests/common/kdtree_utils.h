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

template<typename Scalar, typename VectorContainer>
bool check_range_neighbors(const VectorContainer& points, const std::vector<int>& sampling, int index, Scalar r, const std::vector<int>& neighbors)
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

template<typename Scalar, typename VectorType, typename VectorContainer>
bool check_range_neighbors(const VectorContainer& points, const std::vector<int>& sampling, const VectorType& point, Scalar r, const std::vector<int>& neighbors)
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

template<typename Scalar, typename VectorContainer>
bool check_nearest_neighbor(const VectorContainer& points, int index, int nearest)
{
	return check_k_nearest_neighbors<Scalar, VectorContainer>(points, index, 1, { nearest });
}
template<typename Scalar, typename VectorType, typename VectorContainer>
bool check_nearest_neighbor(const VectorContainer& points, const VectorType& point, int nearest)
{
	return check_k_nearest_neighbors<Scalar, VectorType, VectorContainer>(points, point, 1, { nearest });
}

template<typename Scalar, typename VectorType, typename VectorContainer>
bool check_nearest_neighbor(const VectorContainer& points, const std::vector<int>& sampling, const VectorType& point, int nearest)
{
	return check_k_nearest_neighbors<Scalar, VectorContainer>(points, sampling, point, 1, { nearest });
}

template<typename Scalar, typename VectorType, typename VectorContainer>
bool check_nearest_neighbor(const VectorContainer& points, const std::vector<int>& sampling, int index, int nearest)
{
	return check_k_nearest_neighbors<Scalar, VectorType, VectorContainer>(points, sampling, index, 1, { nearest });
}

template<typename Scalar, typename VectorContainer>
bool check_k_nearest_neighbors(const VectorContainer& points, int index, int k, const std::vector<int>& neighbors)
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

template<typename Scalar, typename VectorContainer>
bool check_k_nearest_neighbors(const VectorContainer& points, const std::vector<int>& sampling, int index, int k, const std::vector<int>& neighbors)
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

template<typename Scalar, typename VectorType, typename VectorContainer>
bool check_k_nearest_neighbors(const VectorContainer& points, const std::vector<int>& sampling, const VectorType& point, int k, const std::vector<int>& neighbors)
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


template<typename Scalar, typename VectorType, typename VectorContainer>
bool check_k_nearest_neighbors(const VectorContainer& points, const VectorType& point, int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) >= k && int(neighbors.size()) != k)
	{
		return false;
	}

	if (has_duplicate(neighbors))
	{
		return false;
	}

	Scalar max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - point).norm());

	for (int idx = 0; idx<int(points.size()); ++idx)
	{
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
