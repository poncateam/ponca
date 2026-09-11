#include <iostream>
#include <fstream>
#include <iomanip>

#include <Ponca/Ponca>
#include "utils.h"

// Main defines
using Point   = Ponca::PointPositionNormal<double, 3>;
using Factory = Ponca::Factory<Point, Ponca::DistWeightFilter<Point, Ponca::SmoothWeightKernel<typename Point::Scalar>>, Ponca::FitSpaceDer>;

// 
using Scalar = typename Point::Scalar;
using VectorType = typename Point::VectorType;

void GenerateComputeTestCase(
    std::ofstream& file, 
    const std::vector<Point>& points, 
    const std::vector<VectorType> analysisLocations, 
    const std::vector<Scalar>     analysisScale
)
{   
    // Write input
    WriteArray<Point, true>(file, "pointcloud", points); 
    file << ",";

    WriteArray<Point>(file, "analysisLocation", analysisLocations);
    file << ",";

    WriteArray<Point>(file, "analysisScale", analysisScale);
    file << ",";

    file << "\"runs\":[";
    {
        bool first = true;
        // Potential
        Factory::Filter<Ponca::ImplicitPrimitiveProvider>().foreach([&](auto& x){
            std::vector<Scalar> potentials;
            for (unsigned int i = 0; i < analysisLocations.size(); ++i)
            { 
                // Copy to avoid sharing between invocations
                auto fit = x.object;
                fit.init();
                fit.setNeighborFilter({ analysisLocations[i], analysisScale[i] });
                fit.compute(points);

                potentials.push_back(fit.potential());
            }
            
            if (first) first = false;
            else file << ",";

            WriteResult<Point>(file, x.name, "POTENTIAL", std::vector<Scalar>{}, potentials);
        });
    }
    file << "]";
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "At least one output file must be provided" << std::endl;
        return 1;
    }

    srand((unsigned int) 20081976);
    
    std::ofstream file(argv[1]);
    file << std::setprecision(20);

    constexpr unsigned int N = 16; // Cloud size
    constexpr unsigned int L =  2;  // Analysis location count

    file << "{";
    {
        const VectorType center = VectorType::Random();
        const Scalar radius     = std::abs(VectorType::Random()[0]); // Quickest way to get a random float 

        std::vector<Point> cloud(N);
        for (unsigned int i = 0; i < N; ++i)
            cloud[i] = Ponca::getPointOnSphere<Point>(radius, center);

        // Assumes cloud are in [-1, 1], so the vector is also within the boundary
        std::vector<VectorType> analysisLocations(L);  
        std::vector<Scalar>     analysisScalars(L, (Scalar)std::numeric_limits<Scalar>::max());
        for (unsigned int i = 0; i < L; ++i)
            analysisLocations[i] = VectorType::Random(); 

        GenerateComputeTestCase(file, cloud, analysisLocations, analysisScalars);
    }
    file << "}";
    return 0;
}