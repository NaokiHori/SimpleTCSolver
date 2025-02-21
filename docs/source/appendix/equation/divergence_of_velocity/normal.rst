####################
Normalised component
####################

***********
Description
***********

The divergence of the velocity vector reads

.. math::

   \left(
      \sum_j
      \frac{1}{H_j}
      \vec{\hat{E}}_j
      \pder{}{X^j}
   \right)
   \cdot
   \left(
      \sum_i
      \vec{\hat{E}}_i
      \hat{U}^i
   \right)
   =
   \sum_{ij}
   \frac{1}{H_j}
   \vec{\hat{E}}_j
   \cdot
   \pder{\vec{\hat{E}}_i}{X^j}
   \hat{U}^i
   +
   \sum_{ij}
   \frac{1}{H_j}
   \vec{\hat{E}}_j
   \cdot
   \vec{\hat{E}}_i
   \pder{\hat{U}^i}{X^j},

which is

.. math::

   \sum_i
   \frac{1}{H_i}
   \pder{\hat{U}^i}{X^i}
   +
   \sum_{ij}
   \frac{1}{H_i}
   \frac{1}{H_j}
   \pder{H_j}{X^i}
   \hat{U}^i
   -
   \sum_i
   \frac{1}{H_i}
   \frac{1}{H_i}
   \pder{H_i}{X^i}
   \hat{U}^i.

Writing down explicitly would be helpful to find the resulting relation:

.. math::

   \frac{1}{J}
   \sum_i
   \pder{}{X^i}
   \left(
      \frac{J}{H_i}
      \hat{U}^i
   \right).

*******
Example
*******

In cylindrical coordinates, I have

.. math::

   \frac{1}{X^1}
   \pder{X^1 U^1}{X^1}
   +
   \frac{1}{X^1}
   \pder{U^2}{X^2}
   +
   \pder{U^3}{X^3}.

