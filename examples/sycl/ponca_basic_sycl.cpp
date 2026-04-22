#include <sycl/sycl.hpp>

#include <Ponca/Ponca>

class mytask;

/// Convenience function that try to create a GPU device, or a CPU if not available
template <bool verbose = true>
sycl::device initDevice()
{
    sycl::device d;
    try
    {
        d = sycl::device(sycl::gpu_selector_v);
    }
    catch (sycl::exception const& e)
    {
        d = sycl::device(sycl::cpu_selector_v);
    }

    if (verbose)
        std::cout << "Running on: " << d.get_info<sycl::info::device::name>() << "\n";
    return d;
}

using namespace Ponca;
using Scalar     = float;
using MyPoint    = PointPositionNormal<Scalar, 3>;
using WeightFunc = DistWeightFunc<MyPoint, SmoothWeightKernel<Scalar>>;
using Fit        = Basket<MyPoint, WeightFunc, OrientedSphereFit>;

int main(int, char**)
{

    auto queue = sycl::queue(initDevice());

    queue.submit([&](sycl::handler& cgh) {
        auto os = sycl::stream{128, 128, cgh};
        cgh.single_task<mytask>([=]() {
            os << "Hello World! (on device)\n";
            Fit f;
            f.init();
        });
    });

    return 0;
}
