######################################################
Azimuthal momentum transported in the radial direction
######################################################

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vr}} \ur
      }{\vt}
      \dif{\ut}{\vr}
   }{\vr}
   =
   -
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \ave{
         \vr \dt \dz \ur
      }{\vt}
      \dif{\ut}{\vr}
   }{\vr}
   =
   -
   \frac{1}{\vr \dr}
   \ave{
      \ave{
         \vr \ur
      }{\vt}
      \dif{\ut}{\vr}
   }{\vr},

which is

.. math::

   l \vat{\ut}{\timm, \tjc, \tkc}
   +
   c \vat{\ut}{\tic , \tjc, \tkc}
   +
   u \vat{\ut}{\tipp, \tjc, \tkc},

where

.. math::

   l
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\vr_{\tic}}
   \frac{1}{\dr_{\tic}}
   \vat{
      \ave{
         \vr \ur
      }{\vt}
   }{\tim, \tjc, \tkc},

.. math::

   u
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\vr_{\tic}}
   \frac{1}{\dr_{\tic}}
   \vat{
      \ave{
         \vr \ur
      }{\vt}
   }{\tip, \tjc, \tkc},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: uy is transported by ux

