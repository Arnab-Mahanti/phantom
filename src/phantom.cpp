#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <string>
#include <matplot/matplot.h>
#include <highfive/H5Easy.hpp>

namespace h5 = HighFive;
int main()
{
    // Parameters
    const auto grid_size = 10;
    const auto grid_refinement = 5;
    const auto modes = 3;
    const auto nos = 5000;
    const auto T_min = 1000.0;
    const auto T_max = 2000.0;
    const std::string filename = "MAMLAS_test_gaussian_3_2d.hdf5";

    // Random engine and distributions
    std::random_device entropy_source;
    std::mt19937 mt(entropy_source());
    auto dist13 = std::uniform_real_distribution<double>(0.1, 0.3);
    auto dist19 = std::uniform_real_distribution<double>(0.1, 0.9);
    auto dist310 = std::uniform_real_distribution<double>(0.3, 1.0);
    auto dist = std::uniform_real_distribution<double>(0.0, 1.0);

    // Temperature phantoms
    auto T_data = std::vector(nos,
                              matplot::vector_2d(grid_refinement * grid_size,
                                                 matplot::vector_1d(grid_refinement * grid_size, 0.0f)));

    // Meshgrid
    auto x = matplot::linspace(0.0, 1.0, grid_refinement * grid_size);
    auto [X, Y] = matplot::meshgrid(x);
    auto max_T_data = 0.0;

    // Generating phantoms
    for (auto &&T : T_data)
    {
        // Generate initial phantom with absolute positions
        auto peak_x_0 = dist19(mt);
        auto peak_y_0 = dist19(mt);
        auto peak_variance = dist13(mt);
        auto peak_magnitude = dist310(mt);
        auto t = matplot::transform(X, Y, [&](auto x, auto y)
                                    { return peak_magnitude * std::exp(-((x - peak_x_0) * (x - peak_x_0) + (y - peak_y_0) * (y - peak_y_0)) / (2 * peak_variance * peak_variance)); });
        T = matplot::transform(T, t, [](auto a, auto b)
                               { return a + b; });

        // Add N-1 phantoms with relative position
        auto rel_dist_x = std::uniform_real_distribution<double>(-peak_x_0, 1 - peak_x_0);
        auto rel_dist_y = std::uniform_real_distribution<double>(-peak_y_0, 1 - peak_y_0);
        for (size_t n = 1; n < modes; n++)
        {
            auto peak_x = peak_x_0 + rel_dist_x(mt);
            auto peak_y = peak_y_0 + rel_dist_y(mt);
            auto peak_variance = dist13(mt);
            auto peak_magnitude = dist(mt);
            auto t = matplot::transform(X, Y, [&](auto x, auto y)
                                        { return peak_magnitude * std::exp(-((x - peak_x) * (x - peak_x) + (y - peak_y) * (y - peak_y)) / (2 * peak_variance * peak_variance)); });
            T = matplot::transform(T, t, [](auto a, auto b)
                                   { return a + b; });
        }
        // T = matplot::transform(T, [&](auto val)
        //                        { return val / modes; });

        // Store maximum of all T_data
        max_T_data = matplot::max(max_T_data, matplot::max(T));
    }

    // Scale all T_data wrt max_T_data
    for (auto &&T : T_data)
    {
        T = matplot::transform(T, [&max_T_data](auto val)
                               { return val / max_T_data; });
    }

    std::cout << "Maximum of T data: " << max_T_data << std::endl;
    // Store the phantoms
    auto file = h5::File(filename, h5::File::Create | h5::File::Overwrite);
    file.createDataSet("data", T_data);
    file.createDataSet("T_min", T_min);
    file.createDataSet("T_max", T_max);
    matplot::surf(X, Y, T_data[std::uniform_int_distribution(0, nos)(mt)]);
    matplot::show();
}
