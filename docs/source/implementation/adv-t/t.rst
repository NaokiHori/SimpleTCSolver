#########################################################
Azimuthal momentum transported in the azimuthal direction
#########################################################

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vt}} \ut
      }{\vt}
      \dif{\ut}{\vt}
   }{\vt}
   =
   -
   \frac{1}{\vr \dr \dt \dz}
   \ave{
      \ave{
         \dr \dz \ut
      }{\vt}
      \dif{\ut}{\vt}
   }{\vt}
   =
   -
   \frac{1}{\vr \dt}
   \ave{
      \ave{
         \ut
      }{\vt}
      \dif{\ut}{\vt}
   }{\vt},

which is

.. math::

   l \vat{\ut}{\tic, \tjmm, \tkc}
   +
   c \vat{\ut}{\tic, \tjc , \tkc}
   +
   u \vat{\ut}{\tic, \tjpp, \tkc},

where

.. math::

   l
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\vr_{\tic}}
   \frac{1}{\dt}
   \vat{
      \ave{
         \ut
      }{\vt}
   }{\tic, \tjm, \tkc},

.. math::

   u
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\vr_{\tic}}
   \frac{1}{\dt}
   \vat{
      \ave{
         \ut
      }{\vt}
   }{\tic, \tjp, \tkc},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: uy is transported by uy

