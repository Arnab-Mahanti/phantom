#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <ranges>
#include <string>
#include <valarray>
#include <matplot/matplot.h>
#include <highfive/H5Easy.hpp>

namespace range = std::ranges;
namespace h5 = HighFive;
int main()
{
    auto grid_size = 10;
    auto grid_refinement = 5;
    auto modes = 3;
    auto nos = 10;
    auto T_min = 1000.0;
    auto T_max = 2000.0;
    std::string filename = "MAMLAS_test_gaussian_3_2d.hdf5";

    std::random_device entropy_source;
    std::mt19937 mt(entropy_source());
    auto dist = std::uniform_real_distribution<double>();

    auto T_data = std::vector(nos,
                              matplot::vector_2d(grid_refinement * grid_size,
                                                 matplot::vector_1d(grid_refinement * grid_size, 0.0f)));

    auto x = matplot::linspace(0.0, 1.0, grid_refinement * grid_size);
    auto [X, Y] = matplot::meshgrid(x);

    for (auto &&T : T_data)
    {
        for (size_t n = 0; n < modes; n++)
        {
            auto peak_x = dist(mt);
            auto peak_y = dist(mt);
            auto peak_variance = dist(mt);
            auto peak_magnitude = dist(mt);
            auto t = matplot::transform(X, Y, [&](auto x, auto y)
                                        { return peak_magnitude * std::exp(-((x - peak_x) * (x - peak_x) + (y - peak_y) * (y - peak_y)) / (2 * peak_variance * peak_variance)); });
            T = matplot::transform(T, t, [](auto a, auto b)
                                   { return a + b; });
        }
        T = matplot::transform(T, [](auto val)
                               { return val / 3; });
    }

    auto file = h5::File(filename, h5::File::Create | h5::File::Overwrite);
    file.createDataSet("data",T_data);
    file.createDataSet("T_min",T_min);
    file.createDataSet("T_max",T_max);
    matplot::surf(X, Y, T_data[0]);
    matplot::show();
}