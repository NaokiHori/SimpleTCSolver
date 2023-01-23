################
Simple TC Solver
################

UNDER ACTIVE DEVELOPMENT.

More rigorous CI will be added soon, and the documentation will be enriched in the near future in order to clearly tell why and how things are treated.

Also implementations are partially incomplete, in particular 1. the implicit treatment of the diffusive terms in the azimuthal direction and 2. the computation of the dissipation rate are under construction (Cartesian version is directly used).

|License|_ |CI|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleTCSolver
.. _License: https://opensource.org/licenses/MIT

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
This is built on top of `SimpleNSSolver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_.
Please refer to `the documentation <https://naokihori.github.io/SimpleTCSolver>`_ for the governing equations and their discretisations.

**********
Dependency
**********

* `C compiler <https://gcc.gnu.org>`_
* `MPI <https://www.open-mpi.org>`_
* `FFTW3 <https://www.fftw.org>`_

***********
Quick start
***********

Please check the README of `SimpleNSSolver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_.

.. code-block::

   $ make clean
   $ make output # create directories for artifacts
   $ make all    # compile
   $ sh exec.sh  # change parameters if needed

Simulations would dump the whole flow field and some log files, whose visualisation leads

.. image:: https://naokihori.github.io/SimpleTCSolver/_images/result.png
   :width: 100%

Here the averaged azimuthal velocity profile is visualised (left), which shows a disturbed flow field (Taylor vortices).
Also the normalised torques measured on the inner and outer cylinder walls are plotted (right), which should agree to the literature result (e.g., Ostilla et al., J. Fluid Mech. (719), 2013).

Please refer to `the workflow <https://github.com/NaokiHori/SimpleTCSolver/blob/main/.github/workflows/ci.yml>`_ to see how this figure is obtained.

*************
Documentation
*************

Documentation can be found `here <https://naokihori.github.io/SimpleTCSolver>`_ (still on half way).
Since the interface is quite similar to the Cartesian version, descriptions are omitted for simplicity.
Please refer to the `documentation of SimpleNSSolver <https://naokihori.github.io/SimpleNavierStokesSolver>`_ for other details.

