############################
Azimuthal momentum advection
############################

.. include:: /references.txt

I consider variables at :math:`\left( \tic, \tjc, \tkc \right)` if not specified.

Scheme proposed by |MORINISHI2004| is adopted to discretise

.. math::

     \ur L_{\vt \vr}
   + \ut L_{\vt \vt}
   + \uz L_{\vt \vz}
   = \ur \der{\ut}{\vr}
   + \ut \left( \frac{1}{\vr} \der{\ut}{\vt} + \frac{\ur}{\vr} \right)
   + \uz \der{\ut}{\vz}

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

Azimuthal momentum transported in the radial direction:

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vr}} \ur
      }{\vt}
      \dif{\ut}{\vr}
   }{\vr}:

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: uy is transported by ux

Azimuthal momentum transported in the azimuthal direction:

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vt}} \ut
      }{\vt}
      \dif{\ut}{\vt}
   }{\vt}:

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: uy is transported by uy

Azimuthal momentum transported in the axial direction:

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vz}} \uz
      }{\vt}
      \dif{\ut}{\vz}
   }{\vz}:

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: uy is transported by uz

Coriolis effect:

.. math::

   -
   \frac{1}{J}
   \ave{
      \frac{J}{\vr}
      \ave{\ut}{\vt}
      \ave{\ur}{\vr}
   }{\vt}:

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: Coriolis contribution

