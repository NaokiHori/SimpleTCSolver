#################################################
Axial momentum transported in the axial direction
#################################################

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vz}} \uz
      }{\vz}
      \dif{\uz}{\vz}
   }{\vz}
   =
   -
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \ave{
         \vr \dr \dt \uz
      }{\vz}
      \dif{\uz}{\vz}
   }{\vz}
   =
   -
   \frac{1}{\dz}
   \ave{
      \ave{
         \uz
      }{\vz}
      \dif{\uz}{\vz}
   }{\vz},

which is

.. math::

   l \vat{\uz}{\zic, \zjc, \zkmm}
   +
   c \vat{\uz}{\zic, \zjc, \zkc }
   +
   u \vat{\uz}{\zic, \zjc, \zkpp},

where

.. math::

   l
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\Delta \vz}
   \vat{
      \ave{\uz}{\vz}
   }{\zic, \zjc, \zkm},

.. math::

   u
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\Delta \vz}
   \vat{
      \ave{\uz}{\vz}
   }{\zic, \zjc, \zkp},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/compute_rhs/uz.c
   :language: c
   :tag: uz is transported by uz

