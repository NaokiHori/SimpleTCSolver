###################################################
Radial momentum transported in the radial direction
###################################################

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vr}} \ur
      }{\vr}
      \dif{\ur}{\vr}
   }{\vr}
   =
   -
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \ave{
         \vr \dt \dz \ur
      }{\vr}
      \dif{\ur}{\vr}
   }{\vr}
   =
   -
   \frac{1}{\vr \dr}
   \ave{
      \ave{
         \vr \ur
      }{\vr}
      \dif{\ur}{\vr}
   }{\vr},

which is

.. math::

   l \vat{\ur}{\rimm, \rjc, \rkc}
   +
   c \vat{\ur}{\ric,  \rjc, \rkc}
   +
   u \vat{\ur}{\ripp, \rjc, \rkc},

where

.. math::

   l
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\vr_{\ric}}
   \frac{1}{\Delta \vr_{\ric}}
   \vat{
      \ave{\vr \ur}{\vr}
   }{\rim, \rjc, \rkc},

.. math::

   u
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\vr_{\ric}}
   \frac{1}{\Delta \vr_{\ric}}
   \vat{
      \ave{\vr \ur}{\vr}
   }{\rip, \rjc, \rkc},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: ux is transported by ux

