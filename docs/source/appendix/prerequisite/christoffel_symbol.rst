##################
Christoffel symbol
##################

***********
Description
***********

==========
Definition
==========

I consider the Christoffel symbol of the second kind

.. _christoffel:

.. math::

   \christoffel,

which will appear when I consider :ref:`the change of basis vectors <change_of_basis_vector>`.
Here I assume :math:`x^l` is twice-differentiable to find :math:`\Gamma_{ij}^k = \Gamma_{ji}^k`.

==========================
Relation with scale factor
==========================

By using :ref:`the relation for the transformation matrices <jacobi_conv>`:

.. math::

   \jacobiconv,

I have

.. math::

   \Gamma_{ij}^k
   =
   \frac{1}{H_k H_k}
   \sum_l
   \pder{}{X^i}
   \left(
      \pder{x^l}{X^j}
   \right)
   \pder{x^l}{X^k}.

Here I recall :ref:`the orthogonality <orthogonal>`:

.. math::

   \orthogonal

to find

.. _metric_tensor:

.. math::

   \metrictensor,

because

.. math::

   \vec{E}_i
   \cdot
   \vec{E}_j
   &
   =
   \left(
      \sum_k
      \pder{x^k}{X^i}
      \vec{e}_k
   \right)
   \cdot
   \left(
      \sum_l
      \pder{x^l}{X^j}
      \vec{e}_l
   \right) \\
   &
   =
   \sum_{kl}
   \pder{x^k}{X^i}
   \pder{x^l}{X^j}
   \vec{e}_k
   \cdot
   \vec{e}_l \\
   &
   =
   \sum_{kl}
   \pder{x^k}{X^i}
   \pder{x^l}{X^j}
   \delta_{kl} \\
   &
   =
   \sum_k
   \pder{x^k}{X^i}
   \pder{x^k}{X^j}.

Furthermore, I consider its derivative with respect to :math:`X^l`:

.. math::

   \pder{}{X^l}
   \left(
      H_i
      H_j
      \delta_{ij}
   \right)
   &
   =
   \pder{}{X^l}
   \left(
      \sum_k
      \pder{x^k}{X^i}
      \pder{x^k}{X^j}
   \right) \\
   &
   =
   \sum_k
   \pder{}{X^l}
   \left(
      \pder{x^k}{X^i}
   \right)
   \pder{x^k}{X^j}
   +
   \sum_k
   \pder{}{X^l}
   \left(
      \pder{x^k}{X^j}
   \right)
   \pder{x^k}{X^i} \\
   &
   =
   H_j H_j \Gamma_{il}^j
   +
   H_i H_i \Gamma_{jl}^i,

which is equivalent to

.. math::

   \newcommand{\csymb}[3]{
      H_{#1} H_{#1} \Gamma_{#2#3}^{#1}
      =
      \pder{}{X^#3}
      \left(
         H_{#1}
         H_{#2}
         \delta_{#1#2}
      \right)
      -
      H_{#2} H_{#2} \Gamma_{#1#3}^{#2}
   }
   \csymb{i}{j}{l},

   \csymb{j}{l}{i},

   \csymb{l}{i}{j}.

Finally I obtain an explicit relation of the Christoffel symbol of the second kind and the scale factors (the indices are interchanged for later convenience):

.. _from_christoffel_to_scale_factor:

.. math::

   \fromchristoffeltoscalefactor.

================
Useful relations
================

Assuming :math:`k = j` yields

.. math::

   \Gamma_{ij}^j
   =
   \frac{1}{H_j}
   \pder{H_j}{X^i}.

Now I consider

.. math::

   \sum_j
   \Gamma_{ij}^j
   \mathcal{Q},

where :math:`\mathcal{Q}` is an arbitrary order of tensor.

This relation is

.. math::

   \sum_j
   \frac{1}{H_j}
   \pder{H_j}{X^i}
   \mathcal{Q}
   =
   \sum_j
   \frac{1}{J}
   \frac{J}{H_j}
   \pder{H_j}{X^i}
   \mathcal{Q}.

In summary, I have

.. _sum_of_christoffel:

.. math::

   \sumofchristoffel.

