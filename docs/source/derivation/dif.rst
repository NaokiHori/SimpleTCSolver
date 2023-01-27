###############
Diffusive terms
###############

*************************************************
Prerequisite: divergence of a second-order tensor
*************************************************

I consider the divergence of a second-order tensor :math:`T_{ij}` in cylindrical coordinates

.. math::

   \partial_j T_{ij},

namely

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
      & + \left( \ez \otimes \ez \right) T_{\vz \vz}
   \end{pmatrix}.

:math:`27` terms exist in total.
Some of them vanish while others do not.

Since the radial derivative does not affect the basis, terms with :math:`\er \otimes \underline{e}_i` (first column) survive:

.. math::

   \er \der{T_{\vr \vr}}{\vr}
   +
   \et \der{T_{\vt \vr}}{\vr}
   +
   \ez \der{T_{\vz \vr}}{\vr}.

Also, since the axial derivative does not affect the basis, terms with :math:`\ez \otimes \underline{e}_i` (third column) survive:

.. math::

   \er \der{T_{\vr \vz}}{\vz}
   +
   \et \der{T_{\vt \vz}}{\vz}
   +
   \ez \der{T_{\vz \vz}}{\vz}.

Azimuthal derivatives make things complicated.
First, similar to the other components, I have

.. math::

   \er \frac{1}{\vr} \der{T_{\vr \vt}}{\vt}
   +
   \et \frac{1}{\vr} \der{T_{\vt \vt}}{\vt}
   +
   \ez \frac{1}{\vr} \der{T_{\vz \vt}}{\vt}

from the second column.

Second, since the azimuthal derivative of the radial basis yields the azimuthal basis, terms with :math:`\er \otimes \underline{e}_i` (first column) survive:

.. math::

   \er \frac{T_{\vr \vr}}{\vr}
   +
   \et \frac{T_{\vt \vr}}{\vr}
   +
   \ez \frac{T_{\vz \vr}}{\vr}

as well as :math:`\et \otimes \er T_{\vr \vt}`:

.. math::

   \et \frac{T_{\vr \vt}}{\vr}.

Finally, since the azimuthal derivative of the azimuthal basis yields the radial basis with negative sign, :math:`\et \otimes \et T_{\vt \vt}` survives:

.. math::

   - \er \frac{T_{\vt \vt}}{\vr}.

In summary, I have

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

in the radial direction,

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

in the azimuthal direction, and

.. math::

   \der{T_{\vz \vr}}{\vr}
   +
   \frac{1}{\vr}
   \der{T_{\vz \vt}}{\vt}
   +
   \der{T_{\vz \vz}}{\vz}
   +
   \frac{T_{\vz \vr}}{\vr}

in the axial direction.

***************
Diffusive terms
***************

By assuming Newtonian fluid, I can model the stress tensor :math:`T_{ij}` using :math:`L_{ij}`:

.. math::

   2 \mu \left( \frac{1}{2} L_{ij} + \frac{1}{2} L_{ji} \right),

giving

.. math::

   \begin{pmatrix}
      T_{\vr \vr} & T_{\vr \vt} & T_{\vr \vz} \\
      T_{\vt \vr} & T_{\vt \vt} & T_{\vt \vz} \\
      T_{\vz \vr} & T_{\vz \vt} & T_{\vz \vz}
   \end{pmatrix}
   =
   \begin{pmatrix}
      \trr & \trt & \trz \\
      \ttr & \ttt & \ttz \\
      \tzr & \tzt & \tzz
   \end{pmatrix},

where I use an identity:

.. math::

   \der{q}{r}
   -
   \frac{q}{r}
   =
   r \der{}{r} \left( \frac{q}{r} \right).

Note that this relation is not always valid from a numerical point of view (after being discretised).

Also, because of the symmetric nature of :math:`T_{ij}`, I have

.. math::

   \frac{T_{\vr \vr}}{\vr}
   -
   \frac{T_{\vt \vt}}{\vr}
   =
   \mu
   \frac{2}{\vr}
   \left\{
      \vr \der{}{\vr} \left(
         \frac{\ur}{\vr}
      \right)
      -
      \frac{1}{\vr}
      \der{\ut}{\vt}
   \right\},

and

.. math::

   \frac{T_{\vr \vt}}{\vr}
   +
   \frac{T_{\vt \vr}}{\vr}
   =
   \mu
   \frac{2}{\vr}
   \left\{
      \vr \der{}{\vr} \left(
         \frac{\ut}{\vr}
      \right)
      +
      \frac{1}{\vr}
      \der{\ur}{\vt}
   \right\}.

Finally I obtain the following relations, the radial component:

.. math::

   \der{}{\vr} \left\{
      \trr
   \right\}
   & +
   \frac{1}{\vr}
   \der{}{\vt} \left\{
      \trt
   \right\}
   +
   \der{}{\vz} \left\{
      \trz
   \right\} \\
   & +
   \mu \frac{2}{\vr}
   \left\{
      \vr \der{}{\vr} \left(
         \frac{\ur}{\vr}
      \right)
      -
      \frac{1}{\vr}
      \der{\ut}{\vt}
   \right\},

the azimuthal component:

.. math::

   \der{}{\vr} \left\{
      \ttr
   \right\}
   & +
   \frac{1}{\vr}
   \der{}{\vt} \left\{
      \ttt
   \right\}
   +
   \der{}{\vz} \left\{
      \ttz
   \right\} \\
   & +
   \mu \frac{2}{\vr}
   \left\{
      \vr \der{}{\vr} \left(
         \frac{\ut}{\vr}
      \right)
      +
      \frac{1}{\vr}
      \der{\ur}{\vt}
   \right\},

and the axial component:

.. math::

   \der{}{\vr} \left\{
      \tzr
   \right\}
   & +
   \frac{1}{\vr}
   \der{}{\vt} \left\{
      \tzt
   \right\}
   +
   \der{}{\vz} \left\{
      \tzz
   \right\} \\
   & +
   \frac{1}{\vr} \left\{
      \tzr
   \right\}.

When the viscosity is constant

.. math::

   \mu
   =
   \mu_0
   =
   const.,

which is of my interest in this project, these relations can be greatly simplified.

Roughly speaking, since

.. math::

   \der{}{x_j} \der{u_j}{x_i}
   =
   \der{}{x_i} \der{u_j}{x_j}
   =
   0,

the anti-symmetric part of the stress tensor vanishes.

Although this is apparent in Cartesian domain, interchanging :math:`\der{}{x_j}` and :math:`\der{}{x_i}` is not always valid numerically, which will be discussed in the :ref:`implementation <implementation>`.

