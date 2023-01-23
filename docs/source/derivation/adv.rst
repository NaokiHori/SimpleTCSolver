###############
Advective terms
###############

I consider

.. math::

   u_j L_{ij},

namely

.. math::

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

:math:`27` terms exist in total, and :math:`9` terms survive using

.. math::

   \underline{e}_i \cdot \left( \underline{e}_j \otimes \underline{e}_k \right)
   =
   \left( \underline{e}_i \cdot \underline{e}_j \right) \underline{e}_k
   =
   \delta_{ij} \underline{e}_k.

In conclusion, I obtain

.. math::

   \ur L_{\vr \vr}
   +
   \ut L_{\vr \vt}
   +
   \uz L_{\vr \vz}

in the radial direction,

.. math::

   \ur L_{\vt \vr}
   +
   \ut L_{\vt \vt}
   +
   \uz L_{\vt \vz}

in the azimuthal direction, and

.. math::

   \ur L_{\vz \vr}
   +
   \ut L_{\vz \vt}
   +
   \uz L_{\vz \vz}

in the axial direction.

Explicitly, I have

.. math::

   \ur \left( \lrr \right)
   +
   \ut \left( \lrt \right)
   +
   \uz \left( \lrz \right),

.. math::

   \ur \left( \ltr \right)
   +
   \ut \left( \ltt \right)
   +
   \uz \left( \ltz \right),

.. math::

   \ur \left( \lzr \right)
   +
   \ut \left( \lzt \right)
   +
   \uz \left( \lzz \right).

