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

.. Although it is not correct numerically, I adopt the relation

.. .. math::

..    \der{q}{\vr}
..    +
..    \frac{q}{\vr}
..    =
..    \frac{1}{\vr} \der{}{\vr} \left( \vr q \right)

.. to yield

.. .. math::

..    \frac{1}{\vr}
..    \der{}{\vr} \left( \vr T_{\vr \vr} \right)
..    +
..    \frac{1}{\vr}
..    \der{T_{\vr \vt}}{\vt}
..    +
..    \der{T_{\vr \vz}}{\vz}
..    -
..    \frac{T_{\vt \vt}}{\vr},

.. .. math::

..    \frac{1}{\vr}
..    \der{}{\vr} \left( \vr T_{\vt \vr} \right)
..    +
..    \frac{1}{\vr}
..    \der{T_{\vt \vt}}{\vt}
..    +
..    \der{T_{\vt \vz}}{\vz}
..    +
..    \frac{T_{\vr \vt}}{\vr},

.. .. math::

..    \frac{1}{\vr}
..    \der{}{\vr} \left( \vr T_{\vz \vr} \right)
..    +
..    \frac{1}{\vr}
..    \der{T_{\vz \vt}}{\vt}
..    +
..    \der{T_{\vz \vz}}{\vz}.

***************
Diffusive terms
***************

By assuming Newtonian fluid, I can model the stress tensor :math:`T_{ij}` using :math:`L_{ij}`:

.. math::

   2 \mu \left( \frac{1}{2} L_{ij} + \frac{1}{2} L_{ji} \right).

Also, because of the symmetric nature of :math:`T_{ij}`, I have

.. math::

   \frac{T_{\vr \vr}}{\vr}
   -
   \frac{T_{\vt \vt}}{\vr}
   =
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
   \frac{2}{\vr}
   \left\{
      \vr \der{}{\vr} \left(
         \frac{\ut}{\vr}
      \right)
      +
      \frac{1}{\vr}
      \der{\ur}{\vt}
   \right\},

where I use an identity:

.. math::

   \der{q}{r}
   -
   \frac{q}{r}
   =
   r \der{}{r} \left( \frac{q}{r} \right).

Finally I obtain the following relations, the radial component:

.. math::

   2 \der{}{\vr} \lrr
   & +
   \frac{1}{\vr}
   \der{}{\vt} \left\{
      % \lrt + \ltr
      \frac{1}{\vr} \der{\ur}{\vt}
      +
      \vr \der{}{\vr} \left( \frac{\ut}{\vr} \right)
   \right\}
   +
   \der{}{\vz} \left(
      \lrz
      +
      \lzr
   \right) \\
   & +
   \frac{2}{\vr}
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

   \der{}{\vr} \left(
      % \ltr + \lrt
      \vr \der{}{\vr} \left( \frac{\ut}{\vr} \right)
      +
      \frac{1}{\vr} \der{\ur}{\vt}
   \right)
   & +
   \frac{2}{\vr}
   \der{}{\vt} \left(
      \ltt
   \right)
   +
   \der{}{\vz} \left(
      \ltz
      +
      \lzt
   \right) \\
   & +
   \frac{2}{\vr}
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

   \der{}{\vr} \left(
      \lzr
      +
      \lrz
   \right)
   & +
   \frac{1}{\vr}
   \der{}{\vt} \left(
      \lzt
      +
      \ltz
   \right)
   +
   2 \der{}{\vz} \left(
      \lzz
   \right) \\
   & +
   \frac{1}{\vr} \left(
      \lzr
      +
      \lrz
   \right).

