###############
Diffusive terms
###############

********
Momentum
********

By assuming Newtonian fluid, I model the stress tensor :math:`T_{ij}` using :math:`L_{ij}`:

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

.. note::

   These relations are not always valid numerically.

Finally I obtain the following relations.

======
Radial
======

.. math::

   & +
   \der{}{\vr} \left\{
      \trr
   \right\} \\
   & +
   \frac{1}{\vr}
   \der{}{\vt} \left\{
      \trt
   \right\} \\
   & +
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
   \right\}.

When the viscosity is constant:

.. math::

   \mu \frac{1}{\vr} \der{}{\vr} \left(
      \vr \der{\ur}{\vr}
   \right)
   +
   \mu \frac{1}{\vr} \der{}{\vt} \left(
      \frac{1}{\vr} \der{\ur}{\vt}
   \right)
   +
   \mu \der{}{\vz} \left(
      \der{\ur}{\vz}
   \right)
   -
   \mu \frac{\ur}{\vr^2}
   -
   \mu \frac{2}{\vr^2} \der{\ut}{\vt}.

=========
Azimuthal
=========

.. math::

   & +
   \der{}{\vr} \left\{
      \ttr
   \right\} \\
   & +
   \frac{1}{\vr}
   \der{}{\vt} \left\{
      \ttt
   \right\} \\
   & +
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
   \right\}.

When the viscosity is constant:

.. math::

   \mu \frac{1}{\vr} \der{}{\vr} \left(
      \vr \der{\ut}{\vr}
   \right)
   +
   \mu \frac{1}{\vr} \der{}{\vt} \left(
      \frac{1}{\vr} \der{\ut}{\vt}
   \right)
   +
   \mu \der{}{\vz} \left(
      \der{\ut}{\vz}
   \right)
   -
   \mu \frac{\ut}{\vr^2}
   +
   \mu \frac{2}{\vr^2} \der{\ur}{\vt}.

=====
Axial
=====

.. math::

   & +
   \der{}{\vr} \left\{
      \tzr
   \right\} \\
   & +
   \frac{1}{\vr}
   \der{}{\vt} \left\{
      \tzt
   \right\} \\
   & +
   \der{}{\vz} \left\{
      \tzz
   \right\} \\
   & +
   \frac{1}{\vr} \left\{
      \tzr
   \right\}.

When the viscosity is constant:

.. math::

   \mu \frac{1}{\vr} \der{}{\vr} \left(
      \vr \der{\uz}{\vr}
   \right)
   +
   \mu \frac{1}{\vr} \der{}{\vt} \left(
      \frac{1}{\vr} \der{\uz}{\vt}
   \right)
   +
   \mu \der{}{\vz} \left(
      \der{\uz}{\vz}
   \right).

.. note::

   Roughly speaking, the constant viscosity case is derived using

   .. math::

      \der{}{x_j} \der{u_j}{x_i}
      =
      \der{}{x_i} \der{u_j}{x_j}
      =
      0.

   Although this exchange is clear in the continuous domain, it is not always valid numerically.

******
Energy
******

.. note::

   For now I consider the constant viscosity cases.

By taking the inner product of the velocity vector and the momentum equation, I have the following.

======
Radial
======

.. math::

   &
   +
   \mu \frac{1}{\vr} \der{}{\vr} \left(
      \vr \der{\kr}{\vr}
   \right)
   +
   \mu \frac{1}{\vr} \der{}{\vt} \left(
      \frac{1}{\vr} \der{\kr}{\vt}
   \right)
   +
   \mu \der{}{\vz} \left(
      \der{\kr}{\vz}
   \right) \\
   &
   -
   \mu \left( \der{\ur}{\vr} \right)^2
   -
   \mu \left( \frac{1}{\vr} \der{\ur}{\vt} \right)^2
   -
   \mu \left( \der{\ur}{\vz} \right)^2
   -
   \mu \left( \frac{\ur}{\vr} \right)^2
   -
   \mu \left( 2 \frac{\ur}{\vr} \frac{1}{\vr} \der{\ut}{\vt} \right).

=========
Azimuthal
=========

.. math::

   &
   +
   \mu \frac{1}{\vr} \der{}{\vr} \left(
      \vr \der{\kt}{\vr}
   \right)
   +
   \mu \frac{1}{\vr} \der{}{\vt} \left(
      \frac{1}{\vr} \der{\kt}{\vt}
   \right)
   +
   \mu \der{}{\vz} \left(
      \der{\kt}{\vz}
   \right) \\
   &
   -
   \mu \left( \der{\ut}{\vr} \right)^2
   -
   \mu \left( \frac{1}{\vr} \der{\ut}{\vt} \right)^2
   -
   \mu \left( \der{\ut}{\vz} \right)^2
   -
   \mu \left( \frac{\ut}{\vr} \right)^2
   +
   \mu \left( 2 \frac{\ut}{\vr} \frac{1}{\vr} \der{\ur}{\vt} \right).

=====
Axial
=====

.. math::

   &
   +
   \mu \frac{1}{\vr} \der{}{\vr} \left(
      \vr \der{\kz}{\vr}
   \right)
   +
   \mu \frac{1}{\vr} \der{}{\vt} \left(
      \frac{1}{\vr} \der{\kz}{\vt}
   \right)
   +
   \mu \der{}{\vz} \left(
      \der{\kz}{\vz}
   \right) \\
   &
   -
   \mu \left( \der{\uz}{\vr} \right)^2
   -
   \mu \left( \frac{1}{\vr} \der{\uz}{\vt} \right)^2
   -
   \mu \left( \der{\uz}{\vz} \right)^2.

