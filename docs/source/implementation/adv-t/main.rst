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
   \vr \Delta \vr \Delta \vt \Delta \vz.

.. toctree::
   :maxdepth: 1

   r
   t
   z
   c

