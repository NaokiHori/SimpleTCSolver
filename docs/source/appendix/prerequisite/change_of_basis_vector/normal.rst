################
Normalised basis
################

***********
Description
***********

By substituting

.. math::

   \vec{E}_i
   =
   H_i
   \vec{\hat{E}}_i

into :ref:`the covariant relation <dedx>`:

.. math::

   \dedx{i}{j}{k},

I obtain

.. math::

   \pder{}{X^j}
   \left(
      H_i
      \vec{\hat{E}}_i
   \right)
   &
   =
   \pder{H_i}{X^j}
   \vec{\hat{E}}_i
   +
   H_i
   \pder{\hat{E}_i}{X^j} \\
   &
   =
   -
   \sum_k
   \delta_{ij}
   \frac{H_j}{2 H_k}
   \pder{H_i}{X^k}
   \vec{\hat{E}}_k
   -
   \sum_k
   \delta_{ij}
   \frac{H_i}{2 H_k}
   \pder{H_j}{X^k}
   \vec{\hat{E}}_k
   +
   \pder{H_i}{X^j}
   \vec{\hat{E}}_i
   +
   \pder{H_j}{X^i}
   \vec{\hat{E}}_j,

or simply (:math:`j = i` is used to unify the first two terms)

.. _dehatdx:

.. math::

   \dehatdx{i}{j}{k}.

This relation is written down explicitly below for an intuitive understanding.

.. math::

   \newcommand{\diag}[2]{
      -
      \frac{1}{H_{#2}}
      \pder{H_{#1}}{X^{#2}}
      \vec{\hat{E}}_{#2}
   }
   \newcommand{\nodiag}[2]{
      \frac{1}{H_{#1}}
      \pder{H_{#2}}{X^{#1}}
      \vec{\hat{E}}_{#2}
   }
   \begin{pmatrix}
      \pder{\vec{\hat{E}}_1}{X^1} & \pder{\vec{\hat{E}}_2}{X^1} & \pder{\vec{\hat{E}}_3}{X^1} \\
      \pder{\vec{\hat{E}}_1}{X^2} & \pder{\vec{\hat{E}}_2}{X^2} & \pder{\vec{\hat{E}}_3}{X^2} \\
      \pder{\vec{\hat{E}}_1}{X^3} & \pder{\vec{\hat{E}}_2}{X^3} & \pder{\vec{\hat{E}}_3}{X^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      \diag{1}{2}
      \diag{1}{3}
      &
      \nodiag{2}{1}
      &
      \nodiag{3}{1}
      \\
      \nodiag{1}{2}
      &
      \diag{2}{1}
      \diag{2}{3}
      &
      \nodiag{3}{2}
      \\
      \nodiag{1}{3}
      &
      \nodiag{2}{3}
      &
      \diag{3}{1}
      \diag{3}{2}
   \end{pmatrix}.

*******
Example
*******

In cylindrical coordinates, I have

.. math::

   \begin{pmatrix}
      \pder{\vec{\hat{E}}_1}{X^1} & \pder{\vec{\hat{E}}_2}{X^1} & \pder{\vec{\hat{E}}_3}{X^1} \\
      \pder{\vec{\hat{E}}_1}{X^2} & \pder{\vec{\hat{E}}_2}{X^2} & \pder{\vec{\hat{E}}_3}{X^2} \\
      \pder{\vec{\hat{E}}_1}{X^3} & \pder{\vec{\hat{E}}_2}{X^3} & \pder{\vec{\hat{E}}_3}{X^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      0 & 0 & 0 \\
      \vec{\hat{E}}_2 & - \vec{\hat{E}}_1 & 0 \\
      0 & 0 & 0 \\
   \end{pmatrix}.

