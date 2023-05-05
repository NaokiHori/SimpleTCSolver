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
   \vr \Delta \vr \Delta \vt \Delta \vz.

.. toctree::
   :maxdepth: 1

   r
   t
   z

