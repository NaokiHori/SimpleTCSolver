#########################
Radial momentum advection
#########################

.. include:: /references.txt

I consider variables at :math:`\left( \ric, \rjc, \rkc \right)` if not specified.

Scheme proposed by |MORINISHI2004| is adopted to discretise

.. math::

     \ur L_{\vr \vr}
   + \ut L_{\vr \vt}
   + \uz L_{\vr \vz}
   = \ur \der{\ur}{\vr}
   + \ut \left( \frac{1}{\vr} \der{\ur}{\vt} - \frac{\ut}{\vr} \right)
   + \uz \der{\ur}{\vz}

in order to conserve the discrete kinetic energy.

Let

.. math::

   h_{\vr}
   \equiv
   \Delta \vr,

.. math::

   h_{\vt}
   \equiv
   \vr \Delta \vt,

.. math::

   h_{\vz}
   \equiv
   \Delta \vz,

and

.. math::
   J
   \equiv
   h_{\vr}
   h_{\vt}
   h_{\vz}
   =
   \vr \Delta \vr \Delta \vt \Delta \vz,

Radial momentum transported in the radial direction:

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vr}} \ur
      }{\vr}
      \dif{\ur}{\vr}
   }{\vr}:

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: ux is transported by ux

Radial momentum transported in the azimuthal direction:

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vt}} \ut
      }{\vr}
      \dif{\ur}{\vt}
   }{\vt}:

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: ux is transported by uy

Radial momentum transported in the axial direction:

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vz}} \uz
      }{\vr}
      \dif{\ur}{\vz}
   }{\vz}:

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: ux is transported by uz

Centrifugal effect:

.. math::

   +
   \frac{1}{J}
   \ave{
      \frac{J}{\vr}
      \ave{\ut}{\vt}
      \ave{\ut}{\vt}
   }{\vr}:

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: centrifugal contribution

