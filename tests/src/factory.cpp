#include <Ponca/Ponca>
#include <Ponca/src/Fitting/Factory/factory.h>

#include "../common/testing.h"
#include "../common/testUtils.h"

int main(int argc, char** argv)
{
    if (!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    using namespace Ponca;
    using Point = PointPosition<double, 2>;
    using NF    = DistWeightFilter<Point, SmoothWeightKernel<double>>;
    using F     = Factory<Point, NF, 0>;

    F::foreach ([](const auto& x) { std::cout << x.name << std::endl; });

    std::cout << "---" << std::endl;
    auto spheres = F::Filter<NotDerivativesProvider>();
    spheres.foreach ([](const auto& x) { std::cout << x.name << std::endl; });

    // auto m = spheres.GetMethod<Method::APSS>();
    // std::cout << m.name << std::endl;
}
