#####################################################
Axial momentum transported in the azimuthal direction
#####################################################

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vt}} \ut
      }{\vz}
      \dif{\uz}{\vt}
   }{\vt}
   =
   -
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \ave{
         \dr \dz \ut
      }{\vz}
      \dif{\uz}{\vt}
   }{\vt}
   =
   -
   \frac{1}{\vr \dt}
   \ave{
      \ave{
         \ut
      }{\vz}
      \dif{\uz}{\vt}
   }{\vt},

which is

.. math::

   l \vat{\uz}{\zic, \zjmm, \zkc}
   +
   c \vat{\uz}{\zic, \zjc , \zkc}
   +
   u \vat{\uz}{\zic, \zjpp, \zkc},

where

.. math::

   l
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\vr_{\zic}}
   \frac{1}{\Delta \vt}
   \vat{
      \ave{\ut}{\vz}
   }{\zic, \zjm, \zkc},

.. math::

   u
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\vr_{\zic}}
   \frac{1}{\Delta \vt}
   \vat{
      \ave{\ut}{\vz}
   }{\zic, \zjp, \zkc},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/compute_rhs/uz.c
   :language: c
   :tag: uz is transported by uy

