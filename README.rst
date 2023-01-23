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
   $ make output # just in case, to create directories for artifacts
   $ make all    # compile
   $ sh exec.sh  # change parameters if needed

.. image:: https://github.com/NaokiHori/SimpleTCSolver/blob/artifacts/result.png
   :width: 100%

*************
Documentation
*************

Documentation can be found `here <https://naokihori.github.io/SimpleTCSolver>`_ (still on half way).
Since the interface is quite similar to the Cartesian version, descriptions are omitted for simplicity.
Please refer to the `documentation of SimpleNSSolver <https://naokihori.github.io/SimpleNavierStokesSolver>`_ for other details.

