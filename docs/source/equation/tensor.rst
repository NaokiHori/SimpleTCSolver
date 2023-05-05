###############
Tensor calculus
###############

*******************
Derivative of basis
*******************

.. math::

   \der{}{\vt} \er
   =
   \et

.. math::

   \der{}{\vt} \et
   =
   -
   \er

********************
Divergence of vector
********************

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

******************
Gradient of vector
******************

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

Each component yields

.. math::

   \left( \er \der{}{\vr} \right) \otimes \left( \er \ur \right)
   =
   \left( \er \otimes \er \right) \der{\ur}{\vr}

.. math::

   \left( \er \der{}{\vr} \right) \otimes \left( \et \ut \right)
   =
   \left( \er \otimes \et \right) \der{\ut}{\vr}

.. math::

   \left( \er \der{}{\vr} \right) \otimes \left( \ez \uz \right)
   =
   \left( \er \otimes \ez \right) \der{\uz}{\vr}

.. math::

   \left( \et \frac{1}{\vr} \der{}{\vt} \right) \otimes \left( \er \ur \right)
   =
   \left( \et \otimes \et \right) \frac{\ur}{\vr}
   +
   \left( \et \otimes \er \right) \frac{1}{\vr} \der{\ur}{\vt}

.. math::

   \left( \et \frac{1}{\vr} \der{}{\vt} \right) \otimes \left( \et \ut \right)
   =
   \left( \et \otimes \er \right) \left( - \frac{\ut}{\vr} \right)
   +
   \left( \et \otimes \et \right) \frac{1}{\vr} \der{\ut}{\vt}

.. math::

   \left( \et \frac{1}{\vr} \der{}{\vt} \right) \otimes \left( \ez \uz \right)
   =
   \left( \et \otimes \ez \right) \frac{1}{\vr} \der{\uz}{\vt}

.. math::

   \left( \ez \der{}{\vz} \right) \otimes \left( \er \ur \right)
   =
   \left( \ez \otimes \er \right) \der{\ur}{\vz}

.. math::

   \left( \ez \der{}{\vz} \right) \otimes \left( \et \ut \right)
   =
   \left( \ez \otimes \et \right) \der{\ut}{\vz}

.. math::

   \left( \ez \der{}{\vz} \right) \otimes \left( \ez \uz \right)
   =
   \left( \ez \otimes \ez \right) \der{\uz}{\vz}

In summary, I define a tensor :math:`L_{ij}`:

.. math::

   \begin{pmatrix}
      L_{\vr \vr} & L_{\vr \vt} & L_{\vr \vz} \\
      L_{\vt \vr} & L_{\vt \vt} & L_{\vt \vz} \\
      L_{\vz \vr} & L_{\vz \vt} & L_{\vz \vz} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      \lrr & \lrt & \lrz \\
      \ltr & \ltt & \ltz \\
      \lzr & \lzt & \lzz \\
   \end{pmatrix}

.. note::

   :math:`L_{ij}` is the scalar components in front of :math:`\underline{e}_j \otimes \underline{e}_i`.

********************
Divergence of tensor
********************

.. math::

   \partial_j T_{ij}

.. math::

   \left(
      \er \der{}{\vr}
      +
      \et \frac{1}{\vr} \der{}{\vt}
      +
      \ez \der{}{\vz}
   \right)
   \cdot
   \begin{pmatrix}
        + \left( \er \otimes \er \right) T_{\vr \vr}
      & + \left( \et \otimes \er \right) T_{\vr \vt}
      & + \left( \ez \otimes \er \right) T_{\vr \vz} \\
        + \left( \er \otimes \et \right) T_{\vt \vr}
      & + \left( \et \otimes \et \right) T_{\vt \vt}
      & + \left( \ez \otimes \et \right) T_{\vt \vz} \\
        + \left( \er \otimes \ez \right) T_{\vz \vr}
      & + \left( \et \otimes \ez \right) T_{\vz \vt}
      & + \left( \ez \otimes \ez \right) T_{\vz \vz} \\
   \end{pmatrix}

Using

.. math::

   \underline{e}_i
   \cdot
   \left(
      \underline{e}_j
      \otimes
      \underline{e}_k
   \right)
   \equiv
   \left(
      \underline{e}_i
      \cdot
      \underline{e}_j
   \right)
   \underline{e}_k
   =
   \delta_{ij}
   \underline{e}_k,

I have the following.

======
Radial
======

.. math::

   \der{T_{\vr \vr}}{\vr}
   +
   \frac{1}{\vr}
   \der{T_{\vr \vt}}{\vt}
   +
   \der{T_{\vr \vz}}{\vz}
   +
   \frac{T_{\vr \vr}}{\vr}
   -
   \frac{T_{\vt \vt}}{\vr}
   =
   \frac{1}{\vr}
   \der{\left( \vr T_{\vr \vr} \right)}{\vr}
   +
   \frac{1}{\vr}
   \der{T_{\vr \vt}}{\vt}
   +
   \der{T_{\vr \vz}}{\vz}
   -
   \frac{T_{\vt \vt}}{\vr}

=========
Azimuthal
=========

.. math::

   \der{T_{\vt \vr}}{\vr}
   +
   \frac{1}{\vr}
   \der{T_{\vt \vt}}{\vt}
   +
   \der{T_{\vt \vz}}{\vz}
   +
   \frac{T_{\vr \vt}}{\vr}
   +
   \frac{T_{\vt \vr}}{\vr}
   =
   \frac{1}{\vr}
   \der{\left( \vr T_{\vt \vr} \right)}{\vr}
   +
   \frac{1}{\vr}
   \der{T_{\vt \vt}}{\vt}
   +
   \der{T_{\vt \vz}}{\vz}
   +
   \frac{T_{\vr \vt}}{\vr}

=====
Axial
=====

.. math::

   \der{T_{\vz \vr}}{\vr}
   +
   \frac{1}{\vr}
   \der{T_{\vz \vt}}{\vt}
   +
   \der{T_{\vz \vz}}{\vz}
   +
   \frac{T_{\vz \vr}}{\vr}
   =
   \frac{1}{\vr}
   \der{\left( \vr T_{\vz \vr} \right)}{\vr}
   +
   \frac{1}{\vr}
   \der{T_{\vz \vt}}{\vt}
   +
   \der{T_{\vz \vz}}{\vz}

