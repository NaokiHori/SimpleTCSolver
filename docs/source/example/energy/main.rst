######
Energy
######

.. include:: /reference/reference.txt

Since the current project adopts the energy-conserving scheme, the total discrete kinetic energy should be conserved.
However, since the explicit Runge-Kutta scheme is used, the discrete kinetic energy is dissipated, which is confirmed here.

.. mydetails:: Configuration

   .. literalinclude:: data/exec.sh
      :language: bash

Total discrete kinetic energy as a function of time:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/energy/energy1.png
   :width: 600

Discrete kinetic energy decay as a function of time step size:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/energy/energy2.png
   :width: 600

The third-order convergence is expected (e.g. |MORINISHI1998|).

