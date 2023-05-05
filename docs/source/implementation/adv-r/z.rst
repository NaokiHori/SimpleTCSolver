##################################################
Radial momentum transported in the axial direction
##################################################

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vz}} \uz
      }{\vr}
      \dif{\ur}{\vz}
   }{\vz}
   =
   -
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \ave{
         \vr \dr \dt \uz
      }{\vr}
      \dif{\ur}{\vz}
   }{\vz}
   =
   -
   \frac{1}{\vr \dr \dz}
   \ave{
      \ave{
         \vr \dr \uz
      }{\vr}
      \dif{\ur}{\vz}
   }{\vz},

which is

.. math::

   l \vat{\ur}{\ric, \rjc, \rkmm}
   +
   c \vat{\ur}{\ric, \rjc, \rkc }
   +
   u \vat{\ur}{\ric, \rjc, \rkpp},

where

.. math::

   l
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\vr_{\ric}}
   \frac{1}{\Delta \vr_{\ric}}
   \frac{1}{\Delta \vz}
   \vat{
      \ave{\vr \dr \uz}{\vr}
   }{\ric, \rjc, \rkm},

.. math::

   u
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\vr_{\ric}}
   \frac{1}{\Delta \vr_{\ric}}
   \frac{1}{\Delta \vz}
   \vat{
      \ave{\vr \dr \uz}{\vr}
   }{\ric, \rjc, \rkp},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: ux is transported by uz

