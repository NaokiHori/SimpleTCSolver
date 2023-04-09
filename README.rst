################
Simple TC Solver
################

|License|_ |CI|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleTCSolver
.. _License: https://opensource.org/license/MIT

.. |CI| image:: https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/ci.yml/badge.svg
.. _CI: https://github.com/NaokiHori/SimpleTCSolver/actions/workflows/ci.yml

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleTCSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleTCSolver/commits/main

.. image:: https://github.com/NaokiHori/SimpleTCSolver/blob/main/docs/source/thumbnail.png
   :width: 100%

********
Overview
********

This library numerically solves the incompressible Navier-Stokes equations in three-dimensional cylindrical domains using a finite-difference method.
In particular, it aims to simulate the motion of fluid between two independently rotating coaxial cylinders, namely Taylor-Couette flows.
This is based on `SimpleNSSolver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_, a massively parallelised Navier-Stokes solver in Cartesian domains.
Please refer to `the documentation <https://naokihori.github.io/SimpleTCSolver>`_ for the governing equations and their discretisations to achieve a conservative and consistent scheme.

********
Features
********

* A second-order-accurate and energy-conserving finite difference scheme, achieving stable integration.
* `MPI parallelisation <https://github.com/NaokiHori/SimpleDecomp>`_
* `Efficient FFT-based direct Poisson solver <https://naokihori.github.io/SimpleNavierStokesSolver/implementation/fluid/compute_potential/main.html>`_
* Explicit/implicit treatments of diffusive terms in all (radial, azimuthal, and axial) directions are easily switchable.

*******
Caveats
*******

* ``x`` and ``y`` in the code

  This library is based on the Cartesian version and I try to keep the interface as similar as possible.
  So ``x`` and ``y`` are used in the code to describe the radial and the azimuthal directions, respectively.
  Also note that the indices ``i``, ``j``, and ``k`` correspond to radial, azimuthal, and axial loops, respectively.
  To see the correspondence between the implementations and the equations, please refer to `the documentation <https://naokihori.github.io/SimpleTCSolver>`_.

* Boundary conditions

  I only consider *idealised* Taylor-Couette flows whose azimuthal and axial directions are both periodic.
  Unfortunately it is non-trivial to change them.

* Singularity at the pole

  The Navier-Stokes equations are singular at the pole ``r = 0``, which is out of focus in this library.

* Fixed parameters

  The radius of the inner cylinder is fixed to 1.
  Also the azimuthal velocity of the inner cylinder is set to 1, while the outer cylinder is assumed to be at rest.
  They can be modified easily (please check the last part of ``src/fluid/boundary_conditions/uy.c``).

**********
Dependency
**********

* `C compiler <https://gcc.gnu.org>`_
* `MPI <https://www.open-mpi.org>`_
* `FFTW3 <https://www.fftw.org>`_

***********
Quick start
***********

First, prepare the initial condition:

.. code-block:: console

   $ cd initial_conditions
   $ make output   # create directories for artifacts
   $ make all      # compile
   $ vi exec.sh    # change parameters if needed
   $ bash exec.sh  # create initial condition
   $ cd ..         # go back to root

Build and execute main solver:

.. code-block:: console

   $ make output   # create directories for artifacts
   $ make all      # compile
   $ bash exec.sh  # change parameters if needed

For more details, please check the README of `SimpleNSSolver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_.

Typical results are as follows.

Velocity fields:

.. image:: https://naokihori.github.io/SimpleTCSolver/_images/snapshot.png
   :width: 60%

Maximum divergence:

.. image:: https://naokihori.github.io/SimpleTCSolver/_images/divergence.png
   :width: 60%

Normalised torque:

.. image:: https://naokihori.github.io/SimpleTCSolver/_images/nusselt.png
   :width: 60%

The black-dashed line is the literature result (Ostilla et al., J. Fluid Mech. (719), 2013).

*************
Documentation
*************

Documentation can be found `here <https://naokihori.github.io/SimpleTCSolver>`_.
Many explanations are omitted in order to avoid duplication.
Please refer to the `documentation of SimpleNSSolver (Cartesian version) <https://naokihori.github.io/SimpleNavierStokesSolver>`_ for more details.

