#######################
Contravariant component
#######################

***********
Description
***********

The divergence of the velocity vector reads

.. math::

   \left(
      \sum_j
      \vec{E}^j
      \pder{}{X^j}
   \right)
   \cdot
   \left(
      \sum_i
      \vec{E}_i
      U^i
   \right)
   =
   \sum_{ij}
   \vec{E}^j
   \cdot
   \pder{\vec{E}_i}{X^j}
   U^i
   +
   \sum_{ij}
   \vec{E}^j
   \cdot
   \vec{E}_i
   \pder{U^i}{X^j}.

The first term leads to

.. math::

   \sum_{ijk}
   \left(
      \vec{E}^j
      \cdot
      \vec{E}_k
      =
      \delta_k^j
   \right)
   \Gamma_{ij}^k
   U^i
   &
   =
   \sum_{ij}
   \Gamma_{ij}^j
   U^i \\
   &
   =
   \frac{1}{J}
   \sum_i
   \pder{J}{X^i}
   U^i,

where :ref:`the relation <sum_of_christoffel>`

.. math::

   \sumofchristoffel

is used.

The second term is

.. math::

   \sum_{ij}
   \delta_i^j
   \pder{U^i}{X^j}
   &
   =
   \sum_i
   \pder{U^i}{X^i} \\
   &
   =
   \frac{1}{J}
   \sum_i
   J
   \pder{U^i}{X^i}.

As a consequence, the sum yields

.. math::

   \frac{1}{J}
   \sum_i
   \pder{}{X^i}
   \left(
      J
      U^i
   \right).

*******
Example
*******

=======================
Cylindrical coordinates
=======================

.. math::

   \frac{1}{X^1}
   \pder{X^1 U^1}{X^1}
   +
   \pder{U^2}{X^2}
   +
   \pder{U^3}{X^3}.

=======================
Rectilinear coordinates
=======================

.. math::

   \frac{1}{J}
   \pder{
      J
      U^1
   }{
      X^1
   }
   +
   \frac{1}{J}
   \pder{
      J
      U^2
   }{
      X^2
   }
   +
   \frac{1}{J}
   \pder{
      J
      U^3
   }{
      X^3
   }.

