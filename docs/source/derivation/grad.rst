#################
Velocity gradient
#################

****************************
Incompressibility constraint
****************************

.. math::

   \der{u_i}{x_i}
   & =
   \left(
      \er               \der{}{\vr}
    + \et \frac{1}{\vr} \der{}{\vt}
    + \ez               \der{}{\vz}
   \right)
   \cdot
   \left(
      \er \ur
    + \et \ut
    + \ez \uz
   \right) \\
   & =
   \frac{1}{\vr} \der{\vr \ur}{\vr}
   +
   \frac{1}{\vr} \der{\ut}{\vt}
   +
   \der{\uz}{\vz}
   =
   0.

************************
Velocity-gradient tensor
************************

I consider

.. math::

   \der{u_i}{x_j}
   =
   \left(
      \er               \der{}{\vr}
    + \et \frac{1}{\vr} \der{}{\vt}
    + \ez               \der{}{\vz}
   \right)
   \otimes
   \left(
      \er \ur
    + \et \ut
    + \ez \uz
   \right)

and define a tensor :math:`L_{ij}` as

.. math::

   \begin{pmatrix}
      L_{\vr \vr} & L_{\vr \vt} & L_{\vr \vz} \\
      L_{\vt \vr} & L_{\vt \vt} & L_{\vt \vz} \\
      L_{\vz \vr} & L_{\vz \vt} & L_{\vz \vz}
   \end{pmatrix}
   =
   \begin{pmatrix}
      \lrr & \lrt & \lrz \\
      \ltr & \ltt & \ltz \\
      \lzr & \lzt & \lzz
   \end{pmatrix}

Note that :math:`L_{ij}` is the scalar components in front of :math:`\underline{e}_j \otimes \underline{e}_i`.

**********
Derivation
**********

===============
:math:`j = \vr`
===============

.. math::
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \left( \er \der{}{\vr} \right) \otimes \left( \er \ur \right)
   & =
   \left( \er \otimes \er \right) \der{\ur}{\vr} \\
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \left( \er \der{}{\vr} \right) \otimes \left( \et \ut \right)
   & =
   \left( \er \otimes \et \right) \der{\ut}{\vr} \\
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \left( \er \der{}{\vr} \right) \otimes \left( \ez \uz \right)
   & =
   \left( \er \otimes \ez \right) \der{\uz}{\vr}

===============
:math:`j = \vt`
===============

.. math::
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \left( \et \frac{1}{\vr} \der{}{\vt} \right) \otimes \left( \er \ur \right)
   & =
   \left( \et \otimes \et \right) \frac{\ur}{\vr}
   +
   \left( \et \otimes \er \right) \frac{1}{\vr} \der{\ur}{\vt} \\
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \left( \et \frac{1}{\vr} \der{}{\vt} \right) \otimes \left( \et \ut \right)
   & =
   \left( \et \otimes \er \right) \left( - \frac{\ut}{\vr} \right)
   +
   \left( \et \otimes \et \right) \frac{1}{\vr} \der{\ut}{\vt} \\
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \left( \et \frac{1}{\vr} \der{}{\vt} \right) \otimes \left( \ez \uz \right)
   & =
   \left( \et \otimes \ez \right) \frac{1}{\vr} \der{\uz}{\vt}

===============
:math:`j = \vz`
===============

.. math::
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \left( \ez \der{}{\vz} \right) \otimes \left( \er \ur \right)
   & =
   \left( \ez \otimes \er \right) \der{\ur}{\vz} \\
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \left( \ez \der{}{\vz} \right) \otimes \left( \et \ut \right)
   & =
   \left( \ez \otimes \et \right) \der{\ut}{\vz} \\
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   \left( \ez \der{}{\vz} \right) \otimes \left( \ez \uz \right)
   & =
   \left( \ez \otimes \ez \right) \der{\uz}{\vz}

