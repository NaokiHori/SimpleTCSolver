#################
initial_condition
#################

********
Overview
********

This directory contains source files to initialise the flow fields.
Note that this initialiser is not parallelised to ease the configuration of the flow fields.

*************
Configuration
*************

All parameters are specified in ``exec.sh``:

* ``dirname``

   where to output the data

* ``lx``, ``ly``, ``lz``

   domain sizes

* ``glisize``, ``gljsize``, ``glksize``

   resolutions

* ``uniformx``

   use uniform grid spacing in the wall-normal (x) direction or not

Other conditions, such as the logging schedule or the maximum simulation time, are not contained here since they are not relevant to the initial conditions.
Those additional parameters should be specified when run-time (see the main ``exec.sh`` in one-directory up).

*****
Usage
*****

.. code-block:: console

   $ make output
   $ make all
   $ vi exec.sh
   $ bash exec.sh

giving several ``NPY`` files under the specified directory.

These files will be loaded by the main simulator.

