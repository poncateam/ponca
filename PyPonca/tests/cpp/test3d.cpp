#include <iostream>
#include <fstream>
#include <iomanip>

#include <Ponca/Ponca>
#include "utils.h"

// Main defines
using Point   = Ponca::PointPositionNormal<double, 3>;
using Factory = Ponca::Factory<Point, Ponca::DistWeightFilter<Point, Ponca::SmoothWeightKernel<typename Point::Scalar>>,
                               Ponca::FitSpaceDer>;

//
using Scalar     = typename Point::Scalar;
using VectorType = typename Point::VectorType;

void GenerateComputeTestCase(std::ofstream& file, const std::vector<Point>& points,
                             const std::vector<VectorType> analysisLocations, const std::vector<Scalar> analysisScale)
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
        std::vector<Scalar> none{};
        Factory::Filter<Ponca::ImplicitPrimitiveProvider>().foreach ([&](auto& x) {
            std::vector<Scalar> potentials;
            std::vector<VectorType> gradients;
            for (unsigned int i = 0; i < analysisLocations.size(); ++i)
            {
                // Copy to avoid sharing between invocations
                auto fit = x.object;
                fit.init();
                fit.setNeighborFilter({analysisLocations[i], analysisScale[i]});
                fit.compute(points);

                potentials.push_back(fit.potential());
                // With the cast to be consistent with pyponca
                gradients.push_back(fit.implicitPrimitive().primitiveGradient());
            }

            if (first)
                first = false;
            else
                file << ",";

            WriteResult<Point>(file, x.name, "POTENTIAL", x.object.isSigned(), none, potentials);
            file << ",";
            WriteResult<Point>(file, x.name, "PRIMITIVE_GRADIENT", x.object.isSigned(), none, gradients);
        });


        Factory::Filter<Ponca::ProjectionOperatorProvider>().foreach ([&](auto& x) {
            std::vector<VectorType> results;
            for (unsigned int i = 0; i < analysisLocations.size(); ++i)
            {
                // Copy to avoid sharing between invocations
                auto fit = x.object;
                fit.init();
                fit.setNeighborFilter({analysisLocations[i], analysisScale[i]});
                fit.compute(points);
                
                results.push_back(fit.project(analysisLocations[i]));
            }

            if (first)
                first = false;
            else
                file << ",";

            WriteResult<Point>(file, x.name, "PROJECTION", false, analysisLocations, results);
        });

        Factory::Filter<Ponca::PrincipalCurvaturesProvider>().foreach ([&](auto& x) {
            std::vector<std::array<Scalar, 3 + 2 * Point::Dim>> data;
            for (unsigned int i = 0; i < analysisLocations.size(); ++i)
            {
                // Copy to avoid sharing between invocations
                auto fit = x.object;
                fit.init();
                fit.setNeighborFilter({analysisLocations[i], analysisScale[i]});
                fit.compute(points);
                
                data.emplace_back();
                data.back()[0] = fit.kmin();
                data.back()[1] = fit.kmax();
                data.back()[2] = fit.GaussianCurvature();
                
                auto mindir = fit.kminDirection();
                auto maxdir = fit.kmaxDirection();

                for (size_t i = 0; i < Point::Dim; ++i)
                    data.back()[3 + i] = mindir[i];
                for (size_t i = 0; i < Point::Dim; ++i)
                    data.back()[3 + Point::Dim + i] = maxdir[i];
            }

            if (first)
                first = false;
            else
                file << ",";

            WriteResult<Point>(file, x.name, "PRINCIPAL_CURVATURE", false, none, data);
        });

        Factory::Filter<Ponca::AlgebraicSphereProvider>().foreach ([&](auto& x) {
            std::vector<std::array<Scalar, 2 + Point::Dim>> data;
            for (unsigned int i = 0; i < analysisLocations.size(); ++i)
            {
                // Copy to avoid sharing between invocations
                auto fit = x.object;
                fit.init();
                fit.setNeighborFilter({analysisLocations[i], analysisScale[i]});
                fit.compute(points);
                
                data.emplace_back();
                data.back()[0] = fit.radius();
                data.back()[1] = fit.prattNorm();

                auto center = fit.center();
                for (size_t i = 0; i < Point::Dim; ++i)
                    data.back()[2 + i] = center[i];
            }

            if (first)
                first = false;
            else
                file << ",";

            WriteResult<Point>(file, x.name, "ALGEBRAIC_SPHERE", false, none, data);
        });

        Factory::Filter<Ponca::GLSParamProvider>().foreach ([&](auto& x) {
            std::vector<std::array<Scalar, 2 + Point::Dim>> data;
            std::vector<std::array<Scalar, 2 + Point::Dim>> data_normalized;
            
            for (unsigned int i = 0; i < analysisLocations.size(); ++i)
            {
                // Copy to avoid sharing between invocations
                auto fit = x.object;
                fit.init();
                fit.setNeighborFilter({analysisLocations[i], analysisScale[i]});
                fit.compute(points);
                
                data.emplace_back();
                data.back()[0] = fit.tau();
                data.back()[1] = fit.kappa();
                
                data_normalized.emplace_back();
                data_normalized.back()[0] = fit.tau_normalized();
                data_normalized.back()[1] = fit.kappa_normalized();

                auto eta = fit.eta();
                auto eta_normalized = fit.eta_normalized();
                for (size_t i = 0; i < Point::Dim; ++i)
                {
                    data.back()[2 + i] = eta(i);
                    data_normalized.back()[2 + i] = eta_normalized(i);
                }
            }

            if (first)
                first = false;
            else
                file << ",";

            WriteResult<Point>(file, x.name, "GLS", false, none, data);
            file << ",";
            WriteResult<Point>(file, x.name, "GLS_NORMALIZED", false, none, data_normalized);
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

    srand((unsigned int)20081976);

    std::ofstream file(argv[1]);
    file << std::setprecision(20);

    constexpr unsigned int N = 128; // Cloud size
    constexpr unsigned int L = 64;  // Analysis location count

    file << "{";
    {
        const VectorType center = VectorType::Random();
        const Scalar radius     = std::abs(VectorType::Random()[0]); // Quickest way to get a random float

        std::vector<Point> cloud(N);
        for (unsigned int i = 0; i < N; ++i)
            cloud[i] = Ponca::getPointOnSphere<Point>(radius, center);

        // Assumes cloud are in [-1, 1], so the vector is also within the boundary
        std::vector<VectorType> analysisLocations(L);
        std::vector<Scalar> analysisScalars(L, 4 * radius / 3); // Average distance of two points sampled on a sphere
        for (unsigned int i = 0; i < L; ++i)
            analysisLocations[i] = VectorType::Random();

        GenerateComputeTestCase(file, cloud, analysisLocations, analysisScalars);
    }
    file << "}";
    return 0;
}
