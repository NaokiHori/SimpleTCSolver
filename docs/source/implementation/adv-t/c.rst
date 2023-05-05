###############
Coriolis effect
###############

.. math::

   -
   \frac{1}{J}
   \ave{
      \frac{J}{\vr}
      \ave{\ut}{\vt}
      \ave{\ur}{\vr}
   }{\vt}
   =
   -
   \frac{1}{\vr}
   \ave{
      \ave{\ut}{\vt}
      \ave{\ur}{\vr}
   }{\vt}.

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: Coriolis contribution

