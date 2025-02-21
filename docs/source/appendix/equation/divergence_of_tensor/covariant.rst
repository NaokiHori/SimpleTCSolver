###############
Covariant basis
###############

***********
Description
***********

The divergence of a second-order tensor reads

.. math::

   \left(
      \sum_k
      \vec{E}^k
      \pder{}{X^k}
   \right)
   \cdot
   \left\{
      \sum_{ij}
      \left(
         \vec{E}_j
         \otimes
         \vec{E}_i
      \right)
      S^{ij}
   \right\},

which has three components:

.. math::

   \sum_{ijk}
   \left(
      \vec{E}^k
      \cdot
      \pder{\vec{E}_j}{X^k}
   \right)
   \vec{E}_i
   S^{ij}
   =
   &
   \sum_{ijkl}
   \left(
      \vec{E}^k
      \cdot
      \vec{E}_l
      =
      \delta_l^k
   \right)
   \vec{E}_i
   \Gamma_{jk}^l
   S^{ij} \\
   =
   &
   \sum_{ijk}
   \vec{E}_i
   \Gamma_{jk}^k
   S^{ij} \\
   =
   &
   \sum_{ij}
   \vec{E}_i
   \frac{1}{J}
   \pder{J}{X^j}
   S^{ij},

.. math::

   \sum_{ijk}
   \left(
      \vec{E}^k
      \cdot
      \vec{E}_j
      =
      \delta_j^k
   \right)
   \pder{\vec{E}_i}{X^k}
   S^{ij}
   =
   &
   \sum_{ij}
   \pder{\vec{E}_i}{X^j}
   S^{ij} \\
   =
   &
   \sum_{ij}
   \left(
      \frac{1}{H_i}
      \pder{H_i}{X^j}
      \vec{E}_i
      +
      \frac{1}{H_j}
      \pder{H_j}{X^i}
      \vec{E}_j
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
   \right)
   S^{ij} \\
   =
   &
   \newcommand{\tmpa}[2]{
      \frac{1}{H_{#1}}
      \pder{H_{#1}}{X^{#2}}
      \left(
         S^{#1 #2}
         +
         S^{#2 #1}
      \right)
   }
   \newcommand{\tmpb}[2]{
      \frac{H_{#2}}{H_{#1} H_{#1}}
      \pder{H_{#2}}{X^{#1}}
      S^{#2 #2}
   }
   \sum_{ij}
   \vec{E}_i
   \tmpa{i}{j}
   -
   \sum_{ij}
   \vec{E}_i
   \tmpb{i}{j},

.. math::

   \sum_{ijk}
   \left(
      \vec{E}^k
      \cdot
      \vec{E}_j
      =
      \delta_j^k
   \right)
   \vec{E}_i
   \pder{S^{ij}}{X^k}
   =
   \sum_{ij}
   \vec{E}_i
   \pder{S^{ij}}{X^j}.

The sum of the first and the last elements yields

.. math::

   \newcommand{\tmpc}[2]{
      \frac{1}{J}
      \pder{}{X^{#2}}
      \left(
         J
         S^{#1 #2}
      \right)
   }
   \sum_{ij}
   \vec{E}_i
   \tmpc{i}{j}.

Explicitly, the sum of

.. math::

   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      \tmpc{1}{1}
      +
      \tmpc{1}{2}
      +
      \tmpc{1}{3}
      \\
      \tmpc{2}{1}
      +
      \tmpc{2}{2}
      +
      \tmpc{2}{3}
      \\
      \tmpc{3}{1}
      +
      \tmpc{3}{2}
      +
      \tmpc{3}{3}
   \end{pmatrix},

.. math::

   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      \tmpa{1}{1}
      +
      \tmpa{1}{2}
      +
      \tmpa{1}{3}
      \\
      \tmpa{2}{1}
      +
      \tmpa{2}{2}
      +
      \tmpa{2}{3}
      \\
      \tmpa{3}{1}
      +
      \tmpa{3}{2}
      +
      \tmpa{3}{3}
   \end{pmatrix},

.. math::

   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      -
      \tmpb{1}{1}
      -
      \tmpb{1}{2}
      -
      \tmpb{1}{3}
      \\
      -
      \tmpb{2}{1}
      -
      \tmpb{2}{2}
      -
      \tmpb{2}{3}
      \\
      -
      \tmpb{3}{1}
      -
      \tmpb{3}{2}
      -
      \tmpb{3}{3}
   \end{pmatrix}.

*******
Example
*******

=======================
Cylindrical coordinates
=======================

The sum of

.. math::

   \vec{E}_1
   \left(
      \frac{1}{X^1}
      \pder{X^1 S^{1 1}}{X^1}
      +
      \pder{S^{1 2}}{X^2}
      +
      \pder{S^{1 3}}{X^3}
      -
      X^1
      S^{2 2}
   \right),

.. math::

   \vec{E}_2
   \left(
      \frac{1}{X^1}
      \pder{X^1 S^{2 1}}{X^1}
      +
      \pder{S^{2 2}}{X^2}
      +
      \pder{S^{2 3}}{X^3}
      +
      \frac{
         S^{2 1}
         +
         S^{1 2}
      }{
         X^1
      }
   \right),

.. math::

   \vec{E}_3
   \left(
      \frac{1}{X^1}
      \pder{X^1 S^{3 1}}{X^1}
      +
      \pder{S^{3 2}}{X^2}
      +
      \pder{S^{3 3}}{X^3}
   \right).

=======================
Rectilinear coordinates
=======================

The sum of

.. math::

   \vec{E}_1
   \left(
      \frac{1}{J}
      \pder{}{X^1}
      \left(
         J
         S^{11}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^2}
      \left(
         J
         S^{12}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^3}
      \left(
         J
         S^{13}
      \right)
      +
      \frac{1}{H_1}
      \pder{H_1}{X^1}
      S^{11}
   \right),

.. math::

   \vec{E}_2
   \left(
      \frac{1}{J}
      \pder{}{X^1}
      \left(
         J
         S^{21}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^2}
      \left(
         J
         S^{22}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^3}
      \left(
         J
         S^{23}
      \right)
      +
      \frac{1}{H_2}
      \pder{H_2}{X^2}
      S^{22}
   \right),

.. math::

   \vec{E}_3
   \left(
      \frac{1}{J}
      \pder{}{X^1}
      \left(
         J
         S^{31}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^2}
      \left(
         J
         S^{32}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^3}
      \left(
         J
         S^{33}
      \right)
      +
      \frac{1}{H_3}
      \pder{H_3}{X^3}
      S^{33}
   \right).

===========
Application
===========

The sum of

.. math::

   \vec{E}_1
   \left(
      \frac{1}{J}
      \pder{}{X^1}
      \left(
         J
         S^{11}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^2}
      \left(
         J
         S^{12}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^3}
      \left(
         J
         S^{13}
      \right)
      +
      \frac{1}{H_1}
      \pder{H_1}{X^1}
      S^{11}
      -
      \frac{H_2}{H_1 H_1}
      \pder{H_2}{X^1}
      S^{22}
   \right),

.. math::

   \vec{E}_2
   \left(
      \frac{1}{J}
      \pder{}{X^1}
      \left(
         J
         S^{21}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^2}
      \left(
         J
         S^{22}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^3}
      \left(
         J
         S^{23}
      \right)
      +
      \frac{1}{H_2}
      \pder{H_2}{X^1}
      S^{21}
      +
      \frac{1}{H_2}
      \pder{H_2}{X^1}
      S^{12}
   \right),

.. math::

   \vec{E}_3
   \left(
      \frac{1}{J}
      \pder{}{X^1}
      \left(
         J
         S^{31}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^2}
      \left(
         J
         S^{32}
      \right)
      +
      \frac{1}{J}
      \pder{}{X^3}
      \left(
         J
         S^{33}
      \right)
   \right).

