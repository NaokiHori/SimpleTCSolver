##################################################
Axial momentum transported in the radial direction
##################################################

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vr}} \ur
      }{\vz}
      \dif{\uz}{\vr}
   }{\vr}
   =
   -
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \ave{
         \vr \dt \dz \ur
      }{\vz}
      \dif{\uz}{\vr}
   }{\vr}
   =
   -
   \frac{1}{\vr \dr}
   \ave{
      \ave{
         \vr \ur
      }{\vz}
      \dif{\uz}{\vr}
   }{\vr},

which is

.. math::

   l \vat{\ut}{\zimm, \zjc, \zkc}
   +
   c \vat{\ut}{\zic , \zjc, \zkc}
   +
   u \vat{\ut}{\zipp, \zjc, \zkc},

where

.. math::

   l
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\vr_{\zic}}
   \frac{1}{\dr_{\zic}}
   \vat{
      \ave{
         \vr \ur
      }{\vz}
   }{\zim, \zjc, \zkc},

.. math::

   u
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\vr_{\zic}}
   \frac{1}{\dr_{\zic}}
   \vat{
      \ave{
         \vr \ur
      }{\vz}
   }{\zip, \zjc, \zkc},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/compute_rhs/uz.c
   :language: c
   :tag: uz is transported by ux

