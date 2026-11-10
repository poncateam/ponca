#include <Ponca/Ponca>
#include "../common/testing.h"
#include "../common/testUtils.h"
#include "../split_test_helper.h"

template <typename Fit>
void testInterpolation()
{
    using DataPoint  = typename Fit::DataPoint;
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;

    int nbPoints = Eigen::internal::random<int>(100, 1000);

    Scalar radius = Eigen::internal::random<Scalar>(Scalar(0.1), Scalar(10.));

    Scalar analysisScale = Scalar(10.) * std::sqrt(Scalar(4. * M_PI) * radius * radius / nbPoints);
    Scalar centerScale   = Eigen::internal::random<Scalar>(1, 10);
    VectorType center    = VectorType::Random() * centerScale;

    Scalar epsilon = testEpsilon<Scalar>();

    std::vector<DataPoint> vectorPoints(nbPoints);
    for (unsigned int i = 0; i < vectorPoints.size(); ++i)
        vectorPoints[i] = Ponca::getPointOnSphere<DataPoint>(radius, center, false, false);

    Fit fit;
    fit.init();
    fit.compute(vectorPoints.begin(), vectorPoints.end());

    for (size_t i = 0; i < vectorPoints.size(); ++i)
    {
        const auto& p = vectorPoints[i];
        VERIFY(std::abs(fit.potential(p.pos())) < epsilon);
        VERIFY((fit.primitiveGradient(p.pos()) - p.normal()).norm() < epsilon);
    }
}

template <typename Scalar, int Dim>
void callSubTests()
{
    using Point = Ponca::PointPositionNormal<Scalar, Dim>;
    testInterpolation<Ponca::HRBF<Point>>();
}

int main(const int argc, char** argv)
{
    if (!init_testing(argc, argv))
        return EXIT_FAILURE;

    std::cout << "Test for the HRBF fit method" << std::endl;

    callSubTests<float, 3>();
    callSubTests<double, 3>();
}
