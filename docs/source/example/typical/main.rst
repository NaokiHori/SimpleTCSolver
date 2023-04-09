############
Taylor rolls
############

.. include:: /references.txt

To check whether the code works properly, a low Reynolds number case is simulated to see the Taylor rolls.

.. mydetails:: Configuration

   .. literalinclude:: data/exec.sh
      :language: bash

The azimuthal velocity field and the radial-axial velocity field (both averaged in the azimuthal direction) are visualised:

.. image:: data/snapshot.png
   :width: 600

The maximum divergence as a function of time:

.. image:: data/divergence.png
   :width: 600

The normalised torque measured on the inner and the outer cylinders are plotted as a function of time.
The black-dashed line is the result from a literature (|OSTILLA2013|).

.. image:: data/nusselt.png
   :width: 600

