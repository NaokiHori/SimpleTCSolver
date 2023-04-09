########################
Axial momentum advection
########################

.. include:: /references.txt

I consider variables at :math:`\left( \zic, \zjc, \zkc \right)` if not specified.

Scheme proposed by |MORINISHI2004| is adopted to discretise

.. math::

     \ur L_{\vz \vr}
   + \ut L_{\vz \vt}
   + \uz L_{\vz \vz}
   = \ur \der{\uz}{\vr}
   + \ut \frac{1}{\vr} \der{\uz}{\vt}
   + \uz \der{\uz}{\vz}

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

Axial momentum transported in the radial direction:

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vr}} \ur
      }{\vz}
      \dif{\uz}{\vr}
   }{\vr}

.. myliteralinclude:: /../../src/fluid/compute_rhs/uz.c
   :language: c
   :tag: uz is transported by ux

Axial momentum transported in the azimuthal direction:

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vt}} \ut
      }{\vz}
      \dif{\uz}{\vt}
   }{\vt}

.. myliteralinclude:: /../../src/fluid/compute_rhs/uz.c
   :language: c
   :tag: uz is transported by uy

Axial momentum transported in the axial direction:

.. math::

   -
   \frac{1}{J}
   \ave{
      \ave{
         \frac{J}{h_{\vz}} \uz
      }{\vz}
      \dif{\uz}{\vz}
   }{\vz}

.. myliteralinclude:: /../../src/fluid/compute_rhs/uz.c
   :language: c
   :tag: uz is transported by uz

