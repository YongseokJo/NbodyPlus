/**
 * High-performance potential energy computation using pybind11.
 *
 * Compile with:
 *   c++ -O3 -Wall -shared -std=c++17 -fPIC $(python3 -m pybind11 --includes) \
 *       -fopenmp potential_energy.cpp -o potential_energy$(python3-config --extension-suffix)
 *
 * Or use the provided build script: python setup_potential.py build_ext --inplace
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <omp.h>

namespace py = pybind11;

/**
 * Compute gravitational potential energy: -G * sum(mi * mj / rij)
 * In code units where G = 1.
 *
 * Uses OpenMP parallelization for O(n^2) pairwise computation.
 */
double compute_potential_energy(
    py::array_t<double, py::array::c_contiguous> pos,
    py::array_t<double, py::array::c_contiguous> mass
) {
    auto pos_buf = pos.unchecked<2>();
    auto mass_buf = mass.unchecked<1>();

    const ssize_t n = pos_buf.shape(0);
    if (n != mass_buf.shape(0)) {
        throw std::runtime_error("Position and mass arrays must have same length");
    }
    if (pos_buf.shape(1) != 3) {
        throw std::runtime_error("Position array must have shape (N, 3)");
    }

    double potential = 0.0;

    #pragma omp parallel reduction(+:potential)
    {
        #pragma omp for schedule(dynamic, 64)
        for (ssize_t i = 0; i < n - 1; ++i) {
            double local_sum = 0.0;
            const double xi = pos_buf(i, 0);
            const double yi = pos_buf(i, 1);
            const double zi = pos_buf(i, 2);
            const double mi = mass_buf(i);

            for (ssize_t j = i + 1; j < n; ++j) {
                const double dx = xi - pos_buf(j, 0);
                const double dy = yi - pos_buf(j, 1);
                const double dz = zi - pos_buf(j, 2);
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);

                if (r > 0.0) {
                    local_sum += mi * mass_buf(j) / r;
                }
            }
            potential += local_sum;
        }
    }

    return -potential;
}

PYBIND11_MODULE(potential_energy, m) {
    m.doc() = "High-performance potential energy computation";
    m.def("compute_potential_energy", &compute_potential_energy,
          py::arg("pos"), py::arg("mass"),
          R"doc(
Compute gravitational potential energy.

Parameters:
    pos: numpy array of shape (N, 3) with particle positions
    mass: numpy array of shape (N,) with particle masses

Returns:
    Potential energy (negative value, in code units where G=1)
)doc");
}
