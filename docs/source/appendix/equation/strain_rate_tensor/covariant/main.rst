###############
Covariant basis
###############

***********
Description
***********

========================
Velocity-gradient tensor
========================

The gradient of a velocity vector is

.. math::

   \left(
      \sum_j
      \vec{E}^j
      \pder{}{X^j}
   \right)
   \otimes
   \left(
      \sum_i
      \vec{E}_i
      U^i
   \right)
   =
   \sum_{ij}
   \left(
      \vec{E}^j
      \otimes
      \pder{\vec{E}_i}{X^j}
   \right)
   U^i
   +
   \sum_{ij}
   \left(
      \vec{E}^j
      \otimes
      \vec{E}_i
   \right)
   \pder{U^i}{X^j}.

By using :ref:`the changes of the basis vectors <dedx>`:

.. math::

   \dedx{i}{j}{k},

I find that the vector gradient is given by the sum of the following elements.

.. math::

   \newcommand{\tmpa}[2]{
      \frac{1}{H_{#2} H_{#2}}
      \frac{1}{H_{#1}}
      \pder{H_{#1}}{X^{#2}}
      U^{#1}
   }
   \sum_{ij}
   \left(
      \vec{E}_j
      \otimes
      \vec{E}_i
   \right)
   \tmpa{i}{j}.

.. math::

   \newcommand{\tmpb}[1]{
      \frac{1}{H_{#1} H_{#1}}
      \sum_k
      \frac{1}{H_{#1}}
      \pder{H_{#1}}{X^k}
      U^k
   }
   \sum_j
   \left(
      \vec{E}_j
      \otimes
      \vec{E}_j
   \right)
   \frac{1}{H_j H_j}
   \frac{1}{H_j}
   \sum_k
   \pder{H_j}{X^k}
   U^k.

.. math::

   &
   -
   \sum_{ijk}
   \left(
      \vec{E}^j
      \otimes
      \vec{E}_k
   \right)
   \delta_{ij}
   \frac{H_j}{2 H_k H_k}
   \pder{H_i}{X^k}
   U^i
   -
   \sum_{ijk}
   \left(
      \vec{E}^j
      \otimes
      \vec{E}_k
   \right)
   \delta_{ij}
   \frac{H_i}{2 H_k H_k}
   \pder{H_j}{X^k}
   U^i \\
   =
   &
   -
   \sum_{jk}
   \left(
      \vec{E}^j
      \otimes
      \vec{E}_k
   \right)
   \frac{H_j}{2 H_k H_k}
   \pder{H_j}{X^k}
   U^j
   -
   \sum_{ik}
   \left(
      \vec{E}^i
      \otimes
      \vec{E}_k
   \right)
   \frac{H_i}{2 H_k H_k}
   \pder{H_i}{X^k}
   U^i \\
   =
   &
   -
   \sum_{ij}
   \left(
      \vec{E}_j
      \otimes
      \vec{E}_i
   \right)
   \frac{1}{H_j H_j}
   \frac{H_j}{H_i H_i}
   \pder{H_j}{X^i}
   U^j \\
   =
   &
   -
   \sum_{ij}
   \left(
      \vec{E}_j
      \otimes
      \vec{E}_i
   \right)
   \tmpa{j}{i}.

.. math::

   \newcommand{\tmpd}[2]{
      \frac{1}{H_{#2} H_{#2}}
      \pder{U^{#1}}{X^{#2}}
   }
   \sum_{ij}
   \left(
      \vec{E}_j
      \otimes
      \vec{E}_i
   \right)
   \tmpd{i}{j}.

Note that the negative of the third element is the transpose of the first one.

In summary, the velocity gradient is given by the sum of

.. math::

   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      \tmpd{1}{1} & \tmpd{2}{1} & \tmpd{3}{1} \\
      \tmpd{1}{2} & \tmpd{2}{2} & \tmpd{3}{2} \\
      \tmpd{1}{3} & \tmpd{2}{3} & \tmpd{3}{3} \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{E}_1
      \\
      \vec{E}_2
      \\
      \vec{E}_3
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
      \tmpb{1} & 0        & 0        \\
      0        & \tmpb{2} & 0        \\
      0        & 0        & \tmpb{3} \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{E}_1
      \\
      \vec{E}_2
      \\
      \vec{E}_3
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
      0 & \tmpa{2}{1} - \tmpa{1}{2} & \tmpa{3}{1} - \tmpa{1}{3} \\
      \tmpa{1}{2} - \tmpa{2}{1} & 0 & \tmpa{3}{2} - \tmpa{2}{3} \\
      \tmpa{1}{3} - \tmpa{3}{1} & \tmpa{2}{3} - \tmpa{3}{2} & 0 \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{E}_1
      \\
      \vec{E}_2
      \\
      \vec{E}_3
   \end{pmatrix}.

==================
Strain-rate tensor
==================

The strain-rate tensor is defined as the symmetric part of it; namely the sum of

.. math::

   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      \tmpd{1}{1} & \frac{1}{2} \tmpd{2}{1} + \frac{1}{2} \tmpd{1}{2} & \frac{1}{2} \tmpd{3}{1} + \frac{1}{2} \tmpd{1}{3} \\
      sym.        & \tmpd{2}{2} & \frac{1}{2} \tmpd{3}{2} + \frac{1}{2} \tmpd{2}{3} \\
      sym.        & sym.        & \tmpd{3}{3} \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{E}_1
      \\
      \vec{E}_2
      \\
      \vec{E}_3
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
      \tmpb{1} &        0 &        0 \\
      sym.     & \tmpb{2} &        0 \\
      sym.     & sym.     & \tmpb{3} \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{E}_1
      \\
      \vec{E}_2
      \\
      \vec{E}_3
   \end{pmatrix}.

*******
Example
*******

.. mydetails:: Cylindrical coordinates

   .. include:: cylindrical.rst

.. mydetails:: Rectilinear coordinates

   .. include:: rectilinear.rst

.. mydetails:: Application

   .. include:: application.rst

