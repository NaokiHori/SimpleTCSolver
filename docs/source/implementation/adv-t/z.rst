#####################################################
Azimuthal momentum transported in the axial direction
#####################################################

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vz}} \uz
      }{\vt}
      \dif{\ut}{\vz}
   }{\vz}
   =
   -
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \ave{
         \vr \dr \dt \uz
      }{\vt}
      \dif{\ut}{\vz}
   }{\vz}
   =
   -
   \frac{1}{\dz}
   \ave{
      \ave{
         \uz
      }{\vt}
      \dif{\ut}{\vz}
   }{\vz},

which is

.. math::

   l \vat{\ut}{\tic, \tjc, \tkmm}
   +
   c \vat{\ut}{\tic, \tjc, \tkc }
   +
   u \vat{\ut}{\tic, \tjc, \tkpp},

where

.. math::

   l
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\dz}
   \vat{
      \ave{
         \uz
      }{\vt}
   }{\tic, \tjc, \tkm},

.. math::

   u
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\dz}
   \vat{
      \ave{
         \uz
      }{\vt}
   }{\tic, \tjc, \tkp},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: uy is transported by uz

