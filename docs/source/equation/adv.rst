###############
Advective terms
###############

********
Momentum
********

.. math::

   &
   u_j L_{ij}
   = \\
   &
   \left(
      \er \ur
      +
      \et \ut
      +
      \ez \uz
   \right)
   \cdot
   \begin{pmatrix}
        + \left( \er \otimes \er \right) L_{\vr \vr}
      & + \left( \et \otimes \er \right) L_{\vr \vt}
      & + \left( \ez \otimes \er \right) L_{\vr \vz} \\
        + \left( \er \otimes \et \right) L_{\vt \vr}
      & + \left( \et \otimes \et \right) L_{\vt \vt}
      & + \left( \ez \otimes \et \right) L_{\vt \vz} \\
        + \left( \er \otimes \ez \right) L_{\vz \vr}
      & + \left( \et \otimes \ez \right) L_{\vz \vt}
      & + \left( \ez \otimes \ez \right) L_{\vz \vz}
   \end{pmatrix}.

Using

.. math::

   \underline{e}_i \cdot \left( \underline{e}_j \otimes \underline{e}_k \right)
   =
   \left( \underline{e}_i \cdot \underline{e}_j \right) \underline{e}_k
   =
   \delta_{ij} \underline{e}_k,

I have the following.

======
Radial
======

.. math::

   \ur L_{\vr \vr}
   +
   \ut L_{\vr \vt}
   +
   \uz L_{\vr \vz}
   =
   \ur \left( \lrr \right)
   +
   \ut \left( \lrt \right)
   +
   \uz \left( \lrz \right).

=========
Azimuthal
=========

.. math::

   \ur L_{\vt \vr}
   +
   \ut L_{\vt \vt}
   +
   \uz L_{\vt \vz}
   =
   \ur \left( \ltr \right)
   +
   \ut \left( \ltt \right)
   +
   \uz \left( \ltz \right).

=====
Axial
=====

.. math::

   \ur L_{\vz \vr}
   +
   \ut L_{\vz \vt}
   +
   \uz L_{\vz \vz}
   =
   \ur \left( \lzr \right)
   +
   \ut \left( \lzt \right)
   +
   \uz \left( \lzz \right).

******
Energy
******

By taking the inner product of the velocity vector and the momentum equation, I have the following.

======
Radial
======

.. math::

   \ur \der{\kr}{\vr}
   +
   \ut \frac{1}{\vr} \der{\kr}{\vt}
   -
   \frac{\ur \ut^2}{\vr}
   +
   \uz \der{\kr}{\vz},

where :math:`\kr \equiv \ur^2 / 2`.

=========
Azimuthal
=========

.. math::

   \ur \der{\kt}{\vr}
   +
   \ut \frac{1}{\vr} \der{\kt}{\vt}
   +
   \frac{\ur \ut^2}{\vr}
   +
   \uz \der{\kt}{\vz},

where :math:`\kt \equiv \ut^2 / 2`.

=====
Axial
=====

.. math::

   \ur \der{\kz}{\vr}
   +
   \ut \frac{1}{\vr} \der{\kz}{\vt}
   +
   \uz \der{\kz}{\vz},

where :math:`\kz \equiv \uz^2 / 2`.

