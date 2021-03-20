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


template <typename Scalar, typename Vector3>
class _Point {
public:
	typedef Scalar Scalar;
	typedef Vector3 VectorType;
};


bool check_range_neighbors(const Vector3Array& points, const std::vector<int>& sampling, int index, float r, const std::vector<int>& neighbors)
{
	if (has_duplicate(neighbors))
	{
		////PCP_DEBUG_ERROR;
		return false;
	}

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
		{
			////PCP_DEBUG_ERROR;
			return false;
		}
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
	{
		////PCP_DEBUG_ERROR;
		return false;
	}

	for (int i = 0; i<int(neighbors.size()); ++i)
	{
		float dist = (points[neighbors[i]] - points[index]).norm();
		if (r < dist)
		{
			////PCP_DEBUG_ERROR;
			return false;
		}
	}

	for (int i = 0; i<int(sampling.size()); ++i)
	{
		int idx = sampling[i];
		if (idx == index) continue;

		float dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && r < dist)
		{
			////PCP_DEBUG_ERROR;
			return false;
		}
		if (!is_neighbor && dist < r)
		{
			////PCP_DEBUG_ERROR;
			return false;
		}
	}
	return true;
}

bool check_k_nearest_neighbors(const Vector3Array& points, int index, int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) > k && int(neighbors.size()) != k)
	{
		//////PCP_DEBUG_ERROR;
		return false;
	}

	if (has_duplicate(neighbors))
	{
		//////PCP_DEBUG_ERROR;
		return false;
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
	{
		//////PCP_DEBUG_ERROR;
		return false;
	}

	float max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - points[index]).norm());

	for (int idx = 0; idx<int(points.size()); ++idx)
	{
		if (idx == index) continue;

		float dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
		{
			//////PCP_DEBUG_ERROR;
			return false;
		}
		if (!is_neighbor && dist < max_dist)
		{
			//////PCP_DEBUG_ERROR;
			return false;
		}
	}
	return true;
}

bool check_nearest_neighbors(const Vector3Array& points, int index, int nearest)
{
	return check_k_nearest_neighbors(points, index, 1, { nearest });
}

bool check_k_nearest_neighbors(const Vector3Array& points, const std::vector<int>& sampling, int index, int k, const std::vector<int>& neighbors)
{
	if (int(points.size()) > k && int(neighbors.size()) != k)
	{
		//////PCP_DEBUG_ERROR;
		return false;
	}

	if (has_duplicate(neighbors))
	{
		//////PCP_DEBUG_ERROR;
		return false;
	}

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
		{
			//////PCP_DEBUG_ERROR;
			return false;
		}
	}

	auto it = std::find(neighbors.begin(), neighbors.end(), index);
	if (it != neighbors.end())
	{
		//////PCP_DEBUG_ERROR;
		return false;
	}

	float max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - points[index]).norm());

	for (int i = 0; i<int(sampling.size()); ++i)
	{
		int idx = sampling[i];
		if (idx == index) continue;

		float dist = (points[idx] - points[index]).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
		{
			//////PCP_DEBUG_ERROR;
			return false;
		}
		if (!is_neighbor && dist < max_dist)
		{
			//////PCP_DEBUG_ERROR;
			return false;
		}
	}
	return true;
}

bool check_k_nearest_neighbors(const Vector3Array& points, const std::vector<int>& sampling, const Vec3& point, int k, const std::vector<int>& neighbors)
{
	if (int(sampling.size()) >= k && int(neighbors.size()) != k)
	{
		//PCP_DEBUG_ERROR;
		return false;
	}

	if (has_duplicate(neighbors))
	{
		//PCP_DEBUG_ERROR;
		return false;
	}

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
	}

	float max_dist = 0;
	for (int idx : neighbors)
		max_dist = std::max(max_dist, (points[idx] - point).norm());

	for (int i = 0; i<int(sampling.size()); ++i)
	{
		int idx = sampling[i];
		float dist = (points[idx] - point).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && max_dist < dist)
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
		if (!is_neighbor && dist < max_dist)
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
	}
	return true;
}


bool check_nearest_neighbors(const Vector3Array& points, const std::vector<int>& sampling, const Vec3& point, int nearest)
{
	return check_k_nearest_neighbors(points, sampling, point, 1, { nearest });
}



bool check_nearest_neighbors(const Vector3Array& points, const std::vector<int>& sampling, int index, int nearest)
{
	return check_k_nearest_neighbors(points, sampling, index, 1, { nearest });
}

bool check_range_neighbors(const Vector3Array& points, const std::vector<int>& sampling, const Vec3& point, float r, const std::vector<int>& neighbors)
{
	if (has_duplicate(neighbors))
	{
		//PCP_DEBUG_ERROR;
		return false;
	}

	for (int idx : neighbors)
	{
		if (std::find(sampling.begin(), sampling.end(), idx) == sampling.end())
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
	}

	for (int i = 0; i<int(neighbors.size()); ++i)
	{
		float dist = (points[neighbors[i]] - point).norm();
		if (r < dist)
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
	}

	for (int i = 0; i<int(sampling.size()); ++i)
	{
		int idx = sampling[i];
		float dist = (points[idx] - point).norm();
		auto it = std::find(neighbors.begin(), neighbors.end(), idx);
		bool is_neighbor = it != neighbors.end();

		if (is_neighbor && r < dist)
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
		if (!is_neighbor && dist < r)
		{
			//PCP_DEBUG_ERROR;
			return false;
		}
	}
	return true;
}
