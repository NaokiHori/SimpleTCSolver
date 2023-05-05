######################################################
Radial momentum transported in the azimuthal direction
######################################################

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vt}} \ut
      }{\vr}
      \dif{\ur}{\vt}
   }{\vt}
   =
   -
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \ave{
        \dr \dz \ut
      }{\vr}
      \dif{\ur}{\vt}
   }{\vt}
   =
   -
   \frac{1}{\vr \dr \dt}
   \ave{
      \ave{
        \dr \ut
      }{\vr}
      \dif{\ur}{\vt}
   }{\vt},

which is

.. math::

   l \vat{\ur}{\ric, \rjmm, \rkc}
   +
   c \vat{\ur}{\ric, \rjc , \rkc}
   +
   u \vat{\ur}{\ric, \rjpp, \rkc},

where

.. math::

   l
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\vr_{\ric}}
   \frac{1}{\Delta \vr_{\ric}}
   \frac{1}{\Delta \vt}
   \vat{
      \ave{\dr \ut}{\vr}
   }{\ric, \rjm, \rkc},

.. math::

   u
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\vr_{\ric}}
   \frac{1}{\Delta \vr_{\ric}}
   \frac{1}{\Delta \vt}
   \vat{
      \ave{\dr \ut}{\vr}
   }{\ric, \rjp, \rkc},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: ux is transported by uy

