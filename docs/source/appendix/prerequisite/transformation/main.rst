
.. _transformation:

##############
Transformation
##############

***********
Description
***********

======
Vector
======

Again I consider :ref:`the radius vector <vector_comp_basis>`:

.. math::

   \vectorcompbasis,

and its small displacement (which is known as the line segment):

.. math::

   d\vec{r}
   =
   \sum_i
   dx^i
   \vec{e}_i.

By using

.. math::

   dx^i
   =
   \sum_j
   \pder{x^i}{X^j}
   dX^j,

I have

.. math::

   d\vec{r}
   &
   =
   \sum_i
   \left(
      \sum_j
      \pder{x^i}{X^j}
      dX^j
   \right)
   \vec{e}_i \\
   &
   =
   \sum_{ij}
   \pder{x^i}{X^j}
   dX^j
   \vec{e}_i \\
   &
   =
   \sum_j
   dX^j
   \left(
      \sum_i
      \pder{x^i}{X^j}
      \vec{e}_i
   \right).

Since the same line segment can be written as

.. math::

   d\vec{r}
   =
   \sum_j
   dX^j
   \vec{E}_j,

I obtain

.. _from_c_to_g:

.. math::

   \fromctog.

Similarly, I have

.. _from_g_to_c:

.. math::

   \fromgtoc.

The following relations for the normalised vector are obvious from the results.

.. math::

   \vec{\hat{E}}_i
   =
   \frac{1}{H_i}
   \sum_j
   \pder{x^j}{X^i}
   \vec{e}_j,

.. math::

   \vec{e}_i
   =
   \sum_j
   H_j
   \pder{X^j}{x^i}
   \vec{\hat{E}}_j.

=========
Component
=========

Assigning :ref:`the relation <from_g_to_c>`

.. math::

   \fromgtoc

to

.. math::

   \vec{u}
   =
   \sum_i
   u^i
   \vec{e}_i,

and assigning :ref:`the other relation <from_c_to_g>`

.. math::

   \fromctog

to

.. math::

   \vec{u}
   =
   \sum_i
   U^i
   \vec{E}_i

yield

.. math::

   \vec{u}
   =
   \sum_i
   \left(
      \sum_j
      u^j
      \pder{X^i}{x^j}
   \right)
   \vec{E}_i

and

.. math::

   \vec{u}
   =
   \sum_i
   \left(
      \sum_j
      U^j
      \pder{x^i}{X^j}
   \right)
   \vec{e}_i,

respectively.

By taking the inner product, the following relations are obtained:

.. math::

   U^i
   =
   \sum_j
   \pder{X^i}{x^j}
   u^j,

.. math::

   u^i
   =
   \sum_j
   \pder{x^i}{X^j}
   U^j.

The normalised relations are

.. math::

   \hat{U}^i
   =
   H_i
   \sum_j
   \pder{X^i}{x^j}
   u^j,

.. math::

   u^i
   =
   \sum_j
   \frac{1}{H_j}
   \pder{x^i}{X^j}
   \hat{U}^j.

=====================
Transformation matrix
=====================

Taking the inner product of :ref:`the relation <from_c_to_g>`

.. math::

   \fromctog

and :math:`\vec{e}_k` gives

.. math::

   \vec{E}_i
   \cdot
   \vec{e}_k
   &
   =
   \sum_j
   \pder{x^j}{X^i}
   \vec{e}_j
   \cdot
   \vec{e}_k \\
   &
   =
   \sum_j
   \pder{x^j}{X^i}
   \delta_{jk} \\
   &
   =
   \pder{x^k}{X^i},

while the inner product between :ref:`the relation <from_g_to_c>`

.. math::

   \fromgtoc

and :math:`\vec{E}_k` yields

.. math::

   \vec{e}_i
   \cdot
   \vec{E}_k
   &
   =
   \sum_j
   \pder{X^j}{x^i}
   \vec{E}_j
   \cdot
   \vec{E}_k \\
   &
   =
   \sum_j
   \pder{X^j}{x^i}
   H_j
   H_k
   \delta_{jk} \\
   &
   =
   \pder{X^k}{x^i}
   H_k H_k,

where :ref:`the orthogonality <orthogonal>`

.. math::

   \orthogonal

is adopted.

By comparing these two relations (note that the indices are dummy and thus are interchangeable), I obtain

.. _jacobi_conv:

.. math::

   \jacobiconv,

giving the relation of the transformation matrix and its inversed one.

====================
Jacobian determinant
====================

I define the determinant of the transformation matrix :math:`J` as

.. math::

   J
   \equiv
   \prod_i
   H_i.

*******
Example
*******

.. toctree::
   :maxdepth: 1

   cylindrical
   rectilinear
   application

