# Simple TC Solver

[![License](https://img.shields.io/github/license/NaokiHori/SimpleTCSolver)](https://opensource.org/license/MIT)
[![Last Commit](https://img.shields.io/github/last-commit/NaokiHori/SimpleTCSolver/main)](https://github.com/NaokiHori/SimpleTCSolver/commits/main)
[![CI](https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/ci.yml/badge.svg)](https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/ci.yml)
[![Documentation](https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/doc.yml/badge.svg)](https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/doc.yml)

![Thumbnail](https://github.com/NaokiHori/SimpleTCSolver/blob/main/docs/source/thumbnail.png)

## Overview

This library numerically solves the incompressible Navier-Stokes equations in both two-dimensional and three-dimensional planar and curved channels using the finite-difference method. While the primary objective is to simulate Taylor-Couette flows, the solver can also be applied to planar flow simulations. 

For details on the governing equations and numerical methods ensuring conservative and consistent schemes, refer to the [documentation](https://naokihori.github.io/SimpleTCSolver).

## Features

- Energy-consistent treatment of advective, pressure-gradient, and diffusive terms, ensuring proper conservation properties.
- [MPI parallelization](https://github.com/NaokiHori/SimpleDecomp).
- Efficient FFT-based direct Poisson solver.
- Explicit and implicit treatments of diffusive terms in all spatial directions.
- Scalar transport simulation.

## Dependencies

- [C compiler](https://gcc.gnu.org)
- [MPI](https://www.open-mpi.org)
- [FFTW3](https://www.fftw.org)
- [Python3](https://www.python.org) (only required for initializing flow fields in `NPY` format)

## Quick Start

1. **Set up your workspace**

   ```console
   mkdir -p /path/to/your/directory
   cd /path/to/your/directory
   ```

2. **Clone the repository**

   ```console
   git clone --recurse-submodules https://github.com/NaokiHori/SimpleTCSolver
   cd SimpleTCSolver
   ```

3. **Set initial conditions**

   `Python3` is used to generate the initial flow fields as `NPY` files.

   ```console
   cd initial_condition
   make output
   bash main.sh
   cd ..
   ```

4. **Build the solver**

   ```console
   make output
   make all
   ```

## Example Results

### Instantaneous Velocity Field

![Velocity Field](https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/snapshot.png)

### Maximum Divergence

![Divergence](https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/divergence.png)

### Normalized Energy Injection and Dissipation

![Energy Balance](https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/balance_main.png)

The black-dashed line represents reference results from Ostilla et al., *J. Fluid Mech. (719), 2013*.

The numerical scheme is designed to ensure energy injection and dissipation remain perfectly balanced (up to rounding error) in steady states:

![Energy Dissipation](https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/balance_dif.png)

## 2D Version

Since Taylor-Couette flows are inherently three-dimensional, the default solver is designed for 3D simulations. However, a 2D version that extracts radial-azimuthal motions is available in the `2d` branch.

## Planar Flows

This solver can also be used for planar flow simulations (e.g., normal channel flows). To enable this mode, set the `is_curved` flag in the flow initializer to `false`. For more details, refer to the [documentation](https://naokihori.github.io/SimpleTCSolver).

## Acknowledgement

I would like to thank [Dr. Kazuyasu Sugiyama](https://researchmap.jp/50466786) for insightful discussions at *Flow for Future - PoF25* and the *37th CFD Symposium*.

