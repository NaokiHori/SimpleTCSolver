############
Basis vector
############

***********
Description
***********

An arbitrary vector (often referred to as the radius vector) :math:`\vec{r}` is

.. _vector_comp_basis:

.. math::

   \vectorcompbasis,

where :math:`X^i` is the contravariant component, while :math:`\vec{E}_i` is the covariant basis vector defined as

.. math::

   \vec{E}_i
   \equiv
   \pder{\vec{r}}{X^i}.

Note that :math:`\vec{E}_i` is not necessarily normalised; namely

.. math::

   \vec{E}_i
   \cdot
   \vec{E}_i

is not necessarily :math:`1`.
The normalised version :math:`\vec{\hat{E}}_i` is given by

.. math::

   \vec{\hat{E}}_i
   \equiv
   \frac{\vec{E}_i}{H_i},

where :math:`H_i` is the scale factor elaborated :ref:`later <transformation>`.
Since the basis vectors are orthogonal to each other, I have

.. _orthogonal:

.. math::

   \orthogonal,

where :math:`\delta` is the Kronecker delta.

The same vector can be written using the covariant component and the contravariant basis vector as well:

.. math::

   \vec{r}
   =
   \sum_i
   X_i
   \vec{E}^i,

where the covariant and the contravariant vectors satisfy (by definition)

.. math::

   \vec{E}_i
   \cdot
   \vec{E}^j
   \equiv
   \delta_i^j.

Although the covariant vectors are mostly used in this document, the contravariant vector will appear to define :ref:`the nabla operator <nabla>`.

In summary, the three different basis vectors are related by

.. math::

   \vec{\hat{E}}_i
   =
   \frac{\vec{E}_i}{H_i}
   =
   H_i
   \vec{E}^i,

as well as the three components:

.. math::

   \hat{X}_i
   =
   H_i X^i
   =
   \frac{1}{H_i} X_i.

Note that, for the Cartesian coordinate, :math:`h_i \equiv 1` and all the three representations are identical:

.. math::

   \vec{\hat{e}}_i
   \equiv
   \vec{e}_i
   \equiv
   \vec{e}^i,

.. math::

   \hat{x}_i
   \equiv
   x^i
   \equiv
   x_i.

