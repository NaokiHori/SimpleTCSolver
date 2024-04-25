
.. _discrete_incompressibility:

############################
Incompressibility constraint
############################

The incompressibility constraint is defined at each cell center:

.. math::

   \ddiv.

The (left-hand-side of) incompressibility constraint is computed to find the right-hand-side term of the Poisson equation:

.. myliteralinclude:: /../../src/fluid/compute_potential.c
   :language: c
   :tag: local divergence

and to monitor the maximum divergence of the flow field as a logging process:

.. myliteralinclude:: /../../src/logging/divergence.c
   :language: c
   :tag: local divergence

