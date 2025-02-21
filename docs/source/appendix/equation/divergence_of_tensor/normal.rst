################
Normalised basis
################

***********
Description
***********

The divergence of a second-order tensor reads

.. math::

   \left(
      \sum_k
      \frac{1}{H_k}
      \vec{\hat{E}}_k
      \pder{}{X^k}
   \right)
   \cdot
   \left\{
      \sum_{ij}
      \left(
         \vec{\hat{E}}_j
         \otimes
         \vec{\hat{E}}_i
      \right)
      \hat{S}^{ij}
   \right\},

which has the following three components.

The first one is

.. math::

   \sum_{ijk}
   \left(
      \vec{\hat{E}}_k
      \cdot
      \pder{\vec{\hat{E}}_j}{X^k}
   \right)
   \frac{1}{H_k}
   \vec{\hat{E}}_i
   \hat{S}^{ij}.

By using :ref:`the relation <dehatdx>`:

.. math::

   \dehatdx{j}{k}{l},

this leads to

.. math::

   &
   \sum_{ijk}
   \left(
      \vec{\hat{E}}_k
      \cdot
      \frac{1}{H_j}
      \pder{H_k}{X^j}
      \vec{\hat{E}}_k
      -
      \vec{\hat{E}}_k
      \cdot
      \sum_l
      \delta_{jk}
      \frac{1}{H_l}
      \pder{H_j}{X^l}
      \vec{\hat{E}}_l
   \right)
   \frac{1}{H_k}
   \vec{\hat{E}}_i
   \hat{S}^{ij} \\
   =
   &
   \sum_{ijk}
   \vec{\hat{E}}_i
   \frac{1}{H_j}
   \frac{1}{H_k}
   \pder{H_k}{X^j}
   \hat{S}^{ij}
   -
   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_j}
   \frac{1}{H_j}
   \pder{H_j}{X^j}
   \hat{S}^{ij} \\
   =
   &
   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_j}
   \frac{1}{J}
   \pder{J}{X^j}
   \hat{S}^{ij}
   -
   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_j}
   \frac{1}{H_j}
   \pder{H_j}{X^j}
   \hat{S}^{ij},

where :ref:`the relation <sum_of_christoffel>`:

.. math::

   \sumofchristoffel

is adopted.

The second one is

.. math::

   &
   \sum_{ijk}
   \left(
      \vec{\hat{E}}_k
      \cdot
      \vec{\hat{E}}_j
      =
      \delta_{kj}
   \right)
   \frac{1}{H_k}
   \pder{\vec{\hat{E}}_i}{X^k}
   \hat{S}^{ij} \\
   =
   &
   \sum_{ij}
   \frac{1}{H_j}
   \pder{\vec{\hat{E}}_i}{X^j}
   \hat{S}^{ij}.

By using :ref:`the relation <dehatdx>`:

.. math::

   \dehatdx{i}{j}{k},

this leads to

.. math::

   &
   \sum_{ij}
   \frac{1}{H_j}
   \frac{1}{H_i}
   \pder{H_j}{X^i}
   \vec{\hat{E}}_j
   \hat{S}^{ij}
   -
   \sum_{ijk}
   \frac{1}{H_j}
   \delta_{ij}
   \frac{1}{H_k}
   \pder{H_i}{X^k}
   \vec{\hat{E}}_k
   \hat{S}^{ij} \\
   =
   &
   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_i}
   \frac{1}{H_j}
   \pder{H_i}{X^j}
   \hat{S}^{ji}
   -
   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_i}
   \frac{1}{H_j}
   \pder{H_j}{X^i}
   \hat{S}^{jj}.

The third one is

.. math::

   &
   \sum_{ijk}
   \left(
      \vec{\hat{E}}_k
      \cdot
      \vec{\hat{E}}_j
      =
      \delta_{kj}
   \right)
   \frac{1}{H_k}
   \vec{\hat{E}}_i
   \pder{\hat{S}^{ij}}{X^k} \\
   =
   &
   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_j}
   \pder{\hat{S}^{ij}}{X^j}.

The sum of the whole elements yields

.. math::

   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_j}
   \frac{1}{J}
   \pder{
      J
      \hat{S}^{ij}
   }{
      X^j
   }
   -
   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_j}
   \frac{1}{H_j}
   \pder{H_j}{X^j}
   \hat{S}^{ij}
   +
   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_i}
   \frac{1}{H_j}
   \pder{H_i}{X^j}
   \hat{S}^{ji}
   -
   \sum_{ij}
   \vec{\hat{E}}_i
   \frac{1}{H_i}
   \frac{1}{H_j}
   \pder{H_j}{X^i}
   \hat{S}^{jj},

where the part of the first element and the last element were unified to yield the Jacobian relation.

Furthermore, the first two terms can be unified to yield

.. math::

   \newcommand{\tmpa}[2]{
      \frac{1}{J}
      \pder{}{X^#2}
      \left(
         \frac{J}{H_#2}
         \hat{S}^{#1 #2}
      \right)
   }
   \newcommand{\tmpb}[2]{
      \frac{1}{H_#1}
      \frac{1}{H_#2}
      \pder{H_#1}{X^#2}
      \hat{S}^{#2 #1}
   }
   \newcommand{\tmpc}[2]{
      \frac{1}{H_#1}
      \frac{1}{H_#2}
      \pder{H_#2}{X^#1}
      \hat{S}^{#2 #2}
   }
   \sum_{ij}
   \vec{\hat{E}}_i
   \tmpa{i}{j}
   +
   \sum_{ij}
   \vec{\hat{E}}_i
   \tmpb{i}{j}
   -
   \sum_{ij}
   \vec{\hat{E}}_i
   \tmpc{i}{j}.

Explicitly, the sum of

.. math::

   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
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
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
   \end{pmatrix}
   \begin{pmatrix}
      \tmpb{1}{1}
      +
      \tmpb{1}{2}
      +
      \tmpb{1}{3}
      \\
      \tmpb{2}{1}
      +
      \tmpb{2}{2}
      +
      \tmpb{2}{3}
      \\
      \tmpb{3}{1}
      +
      \tmpb{3}{2}
      +
      \tmpb{3}{3}
   \end{pmatrix},

.. math::

   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
   \end{pmatrix}
   \begin{pmatrix}
      -
      \tmpc{1}{1}
      -
      \tmpc{1}{2}
      -
      \tmpc{1}{3}
      \\
      -
      \tmpc{2}{1}
      -
      \tmpc{2}{2}
      -
      \tmpc{2}{3}
      \\
      -
      \tmpc{3}{1}
      -
      \tmpc{3}{2}
      -
      \tmpc{3}{3}
   \end{pmatrix}.

*******
Example
*******

The sum of

.. math::

   \vec{\hat{E}}_1
   \left(
      \frac{1}{X^1}
      \pder{X^1 \hat{S}^{1 1}}{X^1}
      +
      \frac{1}{X^1}
      \pder{\hat{S}^{1 2}}{X^2}
      +
      \pder{\hat{S}^{1 3}}{X^3}
      -
      \frac{
         \hat{S}^{2 2}
      }{
         X^1
      }
   \right),

.. math::

   \vec{\hat{E}}_2
   \left(
      \frac{1}{X^1}
      \pder{X^1 \hat{S}^{2 1}}{X^1}
      +
      \frac{1}{X^1}
      \pder{\hat{S}^{2 2}}{X^2}
      +
      \pder{\hat{S}^{2 3}}{X^3}
      +
      \frac{
         \hat{S}^{1 2}
      }{
         X^1
      }
   \right),

.. math::

   \vec{\hat{E}}_3
   \left(
      \frac{1}{X^1}
      \pder{X^1 \hat{S}^{3 1}}{X^1}
      +
      \frac{1}{X^1}
      \pder{\hat{S}^{3 2}}{X^2}
      +
      \pder{\hat{S}^{3 3}}{X^3}
   \right).

