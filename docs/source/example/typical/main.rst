############
Typical case
############

.. include:: /reference/reference.txt

To check whether the code works properly, a low Reynolds number case is simulated to see the Taylor rolls.

.. mydetails:: Configuration

   .. literalinclude:: data/exec.sh
      :language: bash

The azimuthal velocity field, radial-axial velocity field, and scalar field averaged in the azimuthal direction are visualised:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/snapshot.png
   :width: 600

The maximum divergence as a function of time:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/divergence.png
   :width: 600

The normalised energy injection and dissipation (divided by the laminar reference value) inside the domain is plotted as a function of time.
Note that these normalised energy budgets coincide with the definition of the Nusselt number of the angular velocity flux (|ECKHARDT2007|), and the black-dashed line is displayed as a counterpart (|OSTILLA2013|).

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/balance_main.png
   :width: 600

Recall that the schemes adopted in this project are designed so that the energy injection and the dissipation balances statistically, which can be confirmed by the following plot showing the difference between these two quantities as well as an analogous relation with respect to the passive scalar:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleTCSolver/artifacts/artifacts/typical/balance_dif.png
   :width: 600

Initially they are different, which is because the flow is not in a steady state.
Eventually the flow reaches a stationary state, at which the deviations vanish.

