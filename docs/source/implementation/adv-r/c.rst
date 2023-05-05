##################
Centrifugal effect
##################

.. math::

   +
   \frac{1}{J}
   \ave{
      \frac{J}{\vr}
      \ave{\ut}{\vt}
      \ave{\ut}{\vt}
   }{\vr}
   =
   +
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \dr \dt \dz
      \ave{\ut}{\vt}
      \ave{\ut}{\vt}
   }{\vr}
   =
   +
   \frac{1}{\vr \dr}
   \ave{
      \dr
      \ave{\ut}{\vt}
      \ave{\ut}{\vt}
   }{\vr}.

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: centrifugal contribution

