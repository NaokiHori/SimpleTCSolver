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
   \vr \Delta \vr \Delta \vt \Delta \vz.

.. toctree::
   :maxdepth: 1

   r
   t
   z
   c

