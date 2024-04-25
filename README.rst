################
Simple TC Solver
################

|License| |LastCommit|

|CI| |DOC|

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleTCSolver
.. _License: https://opensource.org/license/MIT

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleTCSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleTCSolver/commits/main

.. |CI| image:: https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/ci.yml/badge.svg
.. _CI: https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/ci.yml

.. |DOC| image:: https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/doc.yml/badge.svg
.. _DOC: https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/doc.yml

.. image:: https://github.com/NaokiHori/SimpleTCSolver/blob/main/docs/source/thumbnail.png
   :width: 100%

********
Overview
********

This library numerically solves the incompressible Navier-Stokes equations in two-dimensional and three-dimensional planar and curved channels using a finite-difference method.
Although the main objective is to simulate Taylor-Couette flows, this project also serves as a simulator for planar domains.
The governing equations and their numerical descriptions to achieve conservative and consistent schemes are written in `the documentation <https://naokihori.github.io/SimpleTCSolver>`_.

********
Features
********

* An energy-consistent treatment of advective, pressure-gradient, and diffusive terms, correctly replicating properties of the conservation laws.
* `MPI parallelisation <https://github.com/NaokiHori/SimpleDecomp>`_.
* Efficient FFT-based direct Poisson solver.
* Explicit / implicit treatments of diffusive terms in all spatial directions.
* Scalar transport.

**********
Dependency
**********

* `C compiler <https://gcc.gnu.org>`_
* `MPI <https://www.open-mpi.org>`_
* `FFTW3 <https://www.fftw.org>`_
* `Python3 <https://www.python.org>`_

`Python` is only needed to initialise flow fields as the `NPY` format.

***********
Quick start
***********

#. Prepare workplace

   .. code-block:: console

      mkdir -p /path/to/your/directory
      cd       /path/to/your/directory

#. Get source

   .. code-block:: console

      git clone --recurse-submodules https://github.com/NaokiHori/SimpleTCSolver
      cd SimpleTCSolver

#. Set initial condition

   Here ``Python3`` is used to initialise the flow fields as ``NPY`` files.

   .. code-block:: console

      cd initial_condition
      make output
      bash main.sh
      cd ..

#. Build NS solver

   .. code-block:: console

      make output
      make all

Typical results are as follows.

An instantaneous velocity field:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/snapshot.png
   :width: 60%

Maximum divergence:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/divergence.png
   :width: 60%

Normalised energy injection and dissipation:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/balance_main.png
   :width: 60%

The black-dashed line is the literature result to compare with (Ostilla et al., J. Fluid Mech. (719), 2013).

The numerical scheme is designed such that the energy injection and dissipation perfectly (up to rounding error) balances when the flow fields are in steady states:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/balance_dif.png
   :width: 60%

**********
2D version
**********

Since Taylor-Couette flows are essentially three-dimensional, the three-dimensional version is set as a default.
However, there is a two-dimensional version which extracts the radial-azimuthal motions for completeness, which is available at ``2d`` branch.

************
Planar flows
************

This solver is designed to be used as a solver for planar flows (normal channel flows).
Change ``is_curved`` flag defined in the flow initialiser to ``false``.
See `the documentation <https://naokihori.github.io/SimpleTCSolver>`_ for more details.

***************
Acknowledgement
***************

I would like to acknowledge `Dr. Kazuyasu Sugiyama <https://researchmap.jp/50466786>`_ for fruitful discussions at *Flow for Future - PoF25* and *37th CFD Symposium*.

