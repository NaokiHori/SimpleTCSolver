###############
Covariant basis
###############

***********
Description
***********

I consider

.. math::

   \pder{\vec{E}_i}{X^j}.

By using :ref:`the relation <from_c_to_g>`

.. math::

   \fromctog,

this is equal to

.. math::

   \pder{}{X^j}
   \left(
      \sum_l
      \pder{x^l}{X^i}
      \vec{e}_l
   \right)
   =
   \sum_l
   \pder{}{X^j}
   \left(
      \pder{x^l}{X^i}
   \right)
   \vec{e}_l
   +
   \sum_l
   \pder{x^l}{X^i}
   \pder{\vec{e}_l}{X^j}.

Since the reference Cartesian basis vectors are constant, the second term disappears and I am left with

.. math::

   \sum_l
   \pder{}{X^j}
   \left(
      \pder{x^l}{X^i}
   \right)
   \vec{e}_l
   &
   =
   \sum_l
   \pder{}{X^j}
   \left(
      \pder{x^l}{X^i}
   \right)
   \left(
      \sum_k
      \pder{X^k}{x^l}
      \vec{E}_k
   \right) \\
   &
   =
   \sum_{kl}
   \pder{}{X^j}
   \left(
      \pder{x^l}{X^i}
   \right)
   \pder{X^k}{x^l}
   \vec{E}_k \\
   &
   =
   \sum_k
   \left\{
      \sum_l
      \pder{}{X^j}
      \left(
         \pder{x^l}{X^i}
      \right)
      \pder{X^k}{x^l}
   \right\}
   \vec{E}_k \\
   &
   =
   \sum_k
   \Gamma_{ij}^k
   \vec{E}_k,

where :math:`\Gamma_{ij}^k` is :ref:`the Christoffel symbol of the second kind <christoffel>`:

.. math::

   \christoffel,

and by using :ref:`the relation between the symbol and the scale factors <from_christoffel_to_scale_factor>`:

.. math::

   \fromchristoffeltoscalefactor,

the change of the basis vector is given by

.. math::

   \newcommand{\tmp}[3]{
      \pder{}{X^{#3}}
      \left(
         \delta_{#1#2}
         H_{#1}
         H_{#2}
      \right)
   }
   &
   \sum_k
   \frac{1}{2 H_k H_k}
   \left\{
      -
      \tmp{i}{j}{k}
      +
      \tmp{k}{i}{j}
      +
      \tmp{j}{k}{i}
   \right\}
   \vec{E}_k \\
   =
   &
   -
   \sum_k
   \delta_{ij}
   \frac{H_j}{2 H_k H_k}
   \pder{H_i}{X^k}
   \vec{E}_k
   -
   \sum_k
   \delta_{ij}
   \frac{H_i}{2 H_k H_k}
   \pder{H_j}{X^k}
   \vec{E}_k
   +
   \frac{1}{H_i}
   \pder{H_i}{X^j}
   \vec{E}_i
   +
   \frac{1}{H_j}
   \pder{H_j}{X^i}
   \vec{E}_j,

and in summary:

.. _dedx:

.. math::

   \dedx{i}{j}{k}.

This relation is written down explicitly below for an intuitive understanding.

.. math::

   \newcommand{\diag}[1]{
      \sum_k
      \frac{H_#1}{H_k H_k}
      \pder{H_#1}{X^k}
      \vec{E}_k
   }
   \newcommand{\nodiag}[2]{
      \frac{1}{H_{#1}}
      \pder{H_{#1}}{X^{#2}}
      \vec{E}_{#1}
      +
      \frac{1}{H_{#2}}
      \pder{H_{#2}}{X^{#1}}
      \vec{E}_{#2}
   }
   \begin{pmatrix}
      \pder{\vec{E}_1}{X^1} & \pder{\vec{E}_2}{X^1} & \pder{\vec{E}_3}{X^1} \\
      \pder{\vec{E}_1}{X^2} & \pder{\vec{E}_2}{X^2} & \pder{\vec{E}_3}{X^2} \\
      \pder{\vec{E}_1}{X^3} & \pder{\vec{E}_2}{X^3} & \pder{\vec{E}_3}{X^3} \\
   \end{pmatrix}
   =
   &
   -
   \begin{pmatrix}
      \diag{1} & 0 & 0 \\
      0 & \diag{2} & 0 \\
      0 & 0 & \diag{3} \\
   \end{pmatrix} \\
   &
   +
   \begin{pmatrix}
      0 & \nodiag{1}{2} & \nodiag{1}{3} \\
      sym. & 0 & \nodiag{2}{3} \\
      sym. & sym. & 0 \\
   \end{pmatrix}.

*******
Example
*******

The changes of the basis vectors

.. math::

   \begin{pmatrix}
      \pder{\vec{E}_1}{X^1} & \pder{\vec{E}_2}{X^1} & \pder{\vec{E}_3}{X^1} \\
      \pder{\vec{E}_1}{X^2} & \pder{\vec{E}_2}{X^2} & \pder{\vec{E}_3}{X^2} \\
      \pder{\vec{E}_1}{X^3} & \pder{\vec{E}_2}{X^3} & \pder{\vec{E}_3}{X^3} \\
   \end{pmatrix}

read as follows.

=======================
Cylindrical coordinates
=======================

.. math::

   \begin{pmatrix}
      0 & \frac{1}{X^1} \vec{E}_2 & 0 \\
      \frac{1}{X^1} \vec{E}_2 & - X^1 \vec{E}_1 & 0 \\
      0 & 0 & 0 \\
   \end{pmatrix}.

=======================
Rectilinear coordinates
=======================

.. math::

   \begin{pmatrix}
      - \frac{1}{H_1} \pder{H_1}{X^1} \vec{E}_1 & 0 & 0 \\
      0 & - \frac{1}{H_2} \pder{H_2}{X^2} \vec{E}_2 & 0 \\
      0 & 0 & - \frac{1}{H_3} \pder{H_3}{X^3} \vec{E}_3 \\
   \end{pmatrix}.

===========
Application
===========

.. math::

   \begin{pmatrix}
      - \frac{1}{H_1} \pder{H_1}{X^1} \vec{E}_1 & \frac{1}{H_2} \pder{H_2}{X^1} \vec{E}_2 & 0 \\
      \frac{1}{H_2} \pder{H_2}{X^1} \vec{E}_2 & - \frac{H_2}{H_1 H_1} \pder{H_2}{X^1} \vec{E}_1 & 0 \\
      0 & 0 & 0 \\
   \end{pmatrix}.

