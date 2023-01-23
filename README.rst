################
Simple TC Solver
################

UNDER ACTIVE DEVELOPMENT.
CI will be added soon, and the documentation will be enriched in the near future.

|License|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleTCSolver
.. _License: https://opensource.org/licenses/MIT

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

*************
Documentation
*************

Documentation can be found `here <https://naokihori.github.io/SimpleTCSolver>`_ (still on half way).
Since the interface is quite similar to the Cartesian version, descriptions are omitted for simplicity.
Please refer to the `documentation of SimpleNSSolver <https://naokihori.github.io/SimpleNavierStokesSolver>`_ for other details.

