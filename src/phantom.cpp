#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <ranges>
#include <string>
#include <valarray>
#include <matplot/matplot.h>

namespace range = std::ranges;
int main()
{
    auto grid_size = 10;
    auto grid_refinement = 5;
    auto modes = 3;
    auto nos = 1;
    auto T_min = 1000.0;
    auto T_max = 2000.0;
    std::string filename = "MAMLAS_test_gaussian_3_2d.hdf5";

    std::random_device entropy_source;
    std::mt19937 mt(entropy_source());
    auto dist = std::uniform_real_distribution<double>();

    auto T_data = std::valarray(std::valarray<double>(0.0, grid_refinement * grid_size), nos);
    std::valarray<double> x(grid_refinement * grid_size);
    std::iota(std::begin(x), std::end(x), 0.0);

    range::transform(x, std::begin(x), [&x](auto val)
                     { return val / (x.size() - 1); });

    range::transform(T_data, std::begin(T_data), [&](auto val)
                     {
            for (size_t i = 0; i < modes; i++)
            {
                auto peak_x = dist(mt);
                auto peak_y = dist(mt);
                auto peak_variance = dist(mt);
                val += dist(mt) * std::exp(-((x - peak_x) * (x - peak_x) + (x - peak_y) * (x - peak_y)) / (2 * peak_variance * peak_variance));
            }
    val /= 3;
    return val; });

   

    // range::copy(T_data[0], std::ostream_iterator<double>(std::cout, "\n"));
    // std::cout << "\n";
}