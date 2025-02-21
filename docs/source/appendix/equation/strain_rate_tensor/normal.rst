################
Normalised basis
################

***********
Description
***********

========================
Velocity-gradient tensor
========================

.. math::

   \left(
      \sum_j
      \frac{1}{H_j}
      \vec{\hat{E}}_j
      \pder{}{X^j}
   \right)
   \otimes
   \left(
      \sum_i
      \vec{\hat{E}}_i
      \hat{U}^i
   \right)
   =
   \sum_{ij}
   \left(
      \vec{\hat{E}}_j
      \otimes
      \pder{\vec{\hat{E}}_i}{X^j}
   \right)
   \frac{1}{H_j}
   \hat{U}^i
   +
   \sum_{ij}
   \left(
      \vec{\hat{E}}_j
      \otimes
      \vec{\hat{E}}_i
   \right)
   \frac{1}{H_j}
   \pder{\hat{U}^i}{X^j}.

By using :ref:`the changes of the normalised basis vectors <dehatdx>`:

.. math::

   \dehatdx{i}{j}{k},

the first term yields

.. math::

   +
   \sum_{ij}
   \left(
      \vec{\hat{E}}_j
      \otimes
      \frac{1}{H_i}
      \pder{H_j}{X^i}
      \vec{\hat{E}}_j
   \right)
   \frac{1}{H_j}
   \hat{U}^i
   -
   \sum_{ij}
   \left(
      \vec{\hat{E}}_j
      \otimes
      \sum_k
      \delta_{ij}
      \frac{1}{H_k}
      \pder{H_i}{X^k}
      \vec{\hat{E}}_k
   \right)
   \frac{1}{H_j}
   \hat{U}^i.

By doing arithmetic and simplifying (e.g. interchanging the free indices) the relation, this leads to

.. math::

   +
   \sum_j
   \left(
      \vec{\hat{E}}_j
      \otimes
      \vec{\hat{E}}_j
   \right)
   \sum_k
   \frac{1}{H_j}
   \frac{1}{H_k}
   \pder{H_j}{X^k}
   \hat{U}^k
   -
   \sum_{ij}
   \left(
      \vec{\hat{E}}_j
      \otimes
      \vec{\hat{E}}_i
   \right)
   \frac{1}{H_i}
   \frac{1}{H_j}
   \pder{H_j}{X^i}
   \hat{U}^j.

Thus the conclusive relation is, written explicitly, as follows.

.. math::

   \newcommand{\normal}[2]{
      \frac{1}{H_#2} \pder{\hat{U}^#1}{X^#2}
   }
   \newcommand{\diag}[2]{
      +
      \frac{1}{H_#1 H_#2} \pder{H_#2}{X^#1} \hat{U}^#1
   }
   \newcommand{\nondiag}[2]{
      -
      \frac{1}{H_#1 H_#2} \pder{H_#2}{X^#1} \hat{U}^#2
   }
   &
   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
   \end{pmatrix}
   \begin{pmatrix}
      \normal{1}{1} & \normal{2}{1} & \normal{3}{1} \\
      \normal{1}{2} & \normal{2}{2} & \normal{3}{2} \\
      \normal{1}{3} & \normal{2}{3} & \normal{3}{3} \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{\hat{E}}_1
      \\
      \vec{\hat{E}}_2
      \\
      \vec{\hat{E}}_3
   \end{pmatrix} \\
   &
   + \\
   &
   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
   \end{pmatrix}
   \begin{pmatrix}
      \diag{2}{1}
      \diag{3}{1}
      &
      \nondiag{2}{1}
      &
      \nondiag{3}{1}
      \\
      \nondiag{1}{2}
      &
      \diag{1}{2}
      \diag{3}{2}
      &
      \nondiag{3}{2}
      \\
      \nondiag{1}{3}
      &
      \nondiag{2}{3}
      &
      \diag{1}{3}
      \diag{2}{3}
      \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{\hat{E}}_1
      \\
      \vec{\hat{E}}_2
      \\
      \vec{\hat{E}}_3
   \end{pmatrix}.

==================
Strain-rate tensor
==================

The strain-rate tensor is defined as the symmetric part of it:

.. math::

   \newcommand{\nondiag}[2]{
      -
      \frac{1}{2}
      \frac{1}{H_#1 H_#2}
      \pder{H_#2}{X^#1}
      \hat{U}^#2
   }
   &
   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
   \end{pmatrix}
   \begin{pmatrix}
      \normal{1}{1} & \frac{1}{2} \normal{2}{1} + \frac{1}{2} \normal{1}{2} & \frac{1}{2} \normal{3}{1} + \frac{1}{2} \normal{1}{3} \\
      sym.          & \normal{2}{2} & \frac{1}{2} \normal{3}{2} + \frac{1}{2} \normal{2}{3} \\
      sym.          & sym.          & \normal{3}{3} \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{\hat{E}}_1
      \\
      \vec{\hat{E}}_2
      \\
      \vec{\hat{E}}_3
   \end{pmatrix} \\
   &
   + \\
   &
   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
   \end{pmatrix}
   \begin{pmatrix}
      \diag{2}{1}
      \diag{3}{1}
      &
      \nondiag{2}{1}
      \nondiag{1}{2}
      &
      \nondiag{3}{1}
      \nondiag{1}{3}
      \\
      sym.
      &
      \diag{1}{2}
      \diag{3}{2}
      &
      \nondiag{3}{2}
      \nondiag{2}{3}
      \\
      sym.
      &
      sym.
      &
      \diag{1}{3}
      \diag{2}{3}
      \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{\hat{E}}_1
      \\
      \vec{\hat{E}}_2
      \\
      \vec{\hat{E}}_3
   \end{pmatrix}.

*******
Example
*******

In cylindrical coordinates, the velocity-gradient tensor is given by

.. math::

   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
   \end{pmatrix}
   \begin{pmatrix}
      \pder{\hat{U}^1}{X^1}
      &
      \pder{\hat{U}^2}{X^1}
      &
      \pder{\hat{U}^3}{X^1}
      \\
      \frac{1}{X^1} \pder{\hat{U}^1}{X^2} - \frac{\hat{U}^2}{X^1}
      &
      \frac{1}{X^1} \pder{\hat{U}^2}{X^2} + \frac{\hat{U}^1}{X^1}
      &
      \frac{1}{X^1} \pder{\hat{U}^3}{X^2}
      \\
      \pder{\hat{U}^1}{X^3}
      &
      \pder{\hat{U}^2}{X^3}
      &
      \pder{\hat{U}^3}{X^3}
      \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{\hat{E}}_1
      \\
      \vec{\hat{E}}_2
      \\
      \vec{\hat{E}}_3
   \end{pmatrix}.

The strain-rate tensor is given by

.. math::

   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
   \end{pmatrix}
   \begin{pmatrix}
      \pder{\hat{U}^1}{X^1}
      &
      \frac{1}{2}
      \pder{\hat{U}^2}{X^1}
      +
      \frac{1}{2}
      \frac{1}{X^1}
      \pder{\hat{U}^1}{X^2}
      -
      \frac{1}{2}
      \frac{\hat{U}^2}{X^1}
      &
      \frac{1}{2}
      \pder{\hat{U}^3}{X^1}
      +
      \frac{1}{2}
      \pder{\hat{U}^1}{X^3}
      \\
      sym.
      &
      \frac{1}{X^1}
      \pder{\hat{U}^2}{X^2}
      +
      \frac{\hat{U}^1}{X^1}
      &
      \frac{1}{2}
      \frac{1}{X^1}
      \pder{\hat{U}^3}{X^2}
      +
      \frac{1}{2}
      \pder{\hat{U}^2}{X^3}
      \\
      sym.
      &
      sym.
      &
      \pder{\hat{U}^3}{X^3}
      \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{\hat{E}}_1
      \\
      \vec{\hat{E}}_2
      \\
      \vec{\hat{E}}_3
   \end{pmatrix}.

