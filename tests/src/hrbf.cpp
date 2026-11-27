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

    // Probability of ill-conditionned system is very very very low
    VERIFY(fit.compute(vectorPoints.begin(), vectorPoints.end()) == Ponca::FIT_RESULT::STABLE);

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

    for (int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST((testInterpolation<Ponca::HRBF<Point>>()));
    }
}

using namespace std;
using namespace Ponca;

int main(const int argc, char** argv)
{
    if (!init_testing(argc, argv))
        return EXIT_FAILURE;

    std::cout << "Test for the HRBF fit method" << std::endl;
    CALL_SUBTEST_1((callSubTests<float, 3>()));
    std::cout << " (ok), float" << std::flush;
    CALL_SUBTEST_2((callSubTests<double, 3>()));
    std::cout << " (ok), double" << std::flush;
    CALL_SUBTEST_3((callSubTests<long double, 3>()));
    std::cout << " (ok), long double" << std::flush;
    CALL_SUBTEST_4((callSubTests<double, 2>()));
    std::cout << " (ok), double 2D" << std::flush;

    return 0;
}
