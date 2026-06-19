#include <Ponca/Ponca>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace Ponca;
using namespace std;

using Point      = PointPositionNormal<double, 3>;
using Scalar     = typename Point::Scalar;
using VectorType = typename Point::VectorType;
using NF         = DistWeightFilter<Point, SmoothWeightKernel<Scalar>>;
using MyFactory  = Factory<Point, NF, Ponca::FitSpaceDer>;

int main(int argc, char** argv)
{
    // init input data
    constexpr int n       = 10000;
    constexpr Scalar tmax = Scalar(100.0);

    // set evaluation point and scale
    vector<Point> vecs(n);
    VectorType _p = VectorType::Random();
    std::generate(vecs.begin(), vecs.end(), getRandomPoint<Point>);

    //! [FactoryBasicForeach]
    // Run some code for each element within the factory
    // When using the Factory, state *ARE NOT persistent*
    MyFactory::foreach ([&](auto& x) {
        std::cout << x.name << std::endl;
        x.object.setNeighborFilter({_p, tmax});
        x.object.init();
        x.object.compute(vecs.begin(), vecs.end());
    });
    //! [FactoryBasicForeach]

    //! [FactoryFilterPersistent]
    // We can also filter (this is a conjunction)
    // Here:
    //  * Object that do not provide derivatives
    //  * Builds and implicit representation of the primitive
    //  * Can project any point to this primitive
    auto list = MyFactory::Filter<NotDerivativesProvider, ImplicitPrimitiveProvider, ProjectionOperatorProvider>();
    list.foreach ([&](auto& x) {
        x.object.setNeighborFilter({_p, tmax});
        x.object.init();
        x.object.compute(vecs.begin(), vecs.end());

        std::cout << "Potential for method '" << x.name << "' is: " << x.object.potential() << std::endl;
    });

    // For filtered list, state are persistent if *PASSED BY REFERENCE*
    // (but further filter would lose the current state !):
    list.foreach ([&](auto& x) {
        std::cout << "Projection of:\n" << _p << "\nby '" << x.name << "' is:\n" << x.object.project(_p) << std::endl;
    });
    //! [FactoryFilterPersistent]

    //! [FactoryClassicalMethods]
    // Some classical methods can also be obtained by name.
    auto apss = MyFactory::GetMethod<Method::APSS>();
    apss.setNeighborFilter({_p, tmax});
    apss.init();
    apss.compute(vecs.begin(), vecs.end());
    //! [FactoryClassicalMethods]
}
