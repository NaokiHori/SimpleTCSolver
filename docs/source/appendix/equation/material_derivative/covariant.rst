###############
Covariant basis
###############

***********
Description
***********

Assigning

.. math::

   \vec{u}
   =
   \sum_i
   U^i
   \vec{E}_i

yields

.. math::

   \pder{}{t}
   \left(
      \sum_i
      U^i
      \vec{E}_i
   \right)
   +
   \sum_j
   V^j
   \pder{}{X^j}
   \left(
      \sum_i
      U^i
      \vec{E}_i
   \right)
   =
   \sum_i
   \pder{U^i}{t}
   \vec{E}_i
   +
   \sum_{ij}
   V^j
   \pder{U^i}{X^j}
   \vec{E}_i
   +
   \sum_{ij}
   V^j
   U^i
   \pder{\vec{E}_i}{X^j}.

Regarding the last term, by adopting :ref:`the relation <dedx>`:

.. math::

   \dedx{i}{j}{k},

the last term is reformulated as

.. math::

   &
   \sum_{ij}
   V^j
   U^i
   \frac{1}{H_i}
   \pder{H_i}{X^j}
   \vec{E}_i
   +
   \sum_{ij}
   V^i
   U^j
   \frac{1}{H_i}
   \pder{H_i}{X^j}
   \vec{E}_i
   -
   \sum_{ij}
   V^j
   U^j
   \frac{H_j}{H_i H_i}
   \pder{H_j}{X^i}
   \vec{E}_i \\
   =
   &
   \sum_{ij}
   \left(
      V^j
      U^i
      +
      V^i
      U^j
   \right)
   \frac{1}{H_i}
   \pder{H_i}{X^j}
   \vec{E}_i
   -
   \sum_{ij}
   V^j
   U^j
   \frac{H_j}{H_i H_i}
   \pder{H_j}{X^i}
   \vec{E}_i.

Explicitly, the material derivative is given by the sum of

.. math::

   \newcommand{\tmp}[1]{
      \pder{U^{#1}}{t}
   }
   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      \tmp{1}
      \\
      \tmp{2}
      \\
      \tmp{3}
   \end{pmatrix},

.. math::

   \newcommand{\tmp}[2]{
      V^{#2}
      \pder{U^{#1}}{X^{#2}}
   }
   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      \tmp{1}{1} + \tmp{1}{2} + \tmp{1}{3}
      \\
      \tmp{2}{1} + \tmp{2}{2} + \tmp{2}{3}
      \\
      \tmp{3}{1} + \tmp{3}{2} + \tmp{3}{3}
   \end{pmatrix},

.. math::

   \newcommand{\tmp}[2]{
      \left(
         V^{#2}
         U^{#1}
         +
         V^{#1}
         U^{#2}
      \right)
      \frac{1}{H_{#1}}
      \pder{H_{#1}}{X^{#2}}
   }
   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      \tmp{1}{1} + \tmp{1}{2} + \tmp{1}{3}
      \\
      \tmp{2}{1} + \tmp{2}{2} + \tmp{2}{3}
      \\
      \tmp{3}{1} + \tmp{3}{2} + \tmp{3}{3}
   \end{pmatrix},

.. math::

   \newcommand{\tmp}[2]{
      V^{#2}
      U^{#2}
      \frac{H_{#2}}{H_{#1} H_{#1}}
      \pder{H_{#2}}{X^{#1}}
   }
   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      - \tmp{1}{1} - \tmp{1}{2} - \tmp{1}{3}
      \\
      - \tmp{2}{1} - \tmp{2}{2} - \tmp{2}{3}
      \\
      - \tmp{3}{1} - \tmp{3}{2} - \tmp{3}{3}
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
      \pder{U^1}{t}
      +
      V^1 \pder{U^1}{X^1}
      +
      V^2 \pder{U^1}{X^2}
      +
      V^3 \pder{U^1}{X^3}
      -
      V^2 U^2 X^1
   \right),

.. math::

   \vec{E}_2
   \left(
      \pder{U^2}{t}
      +
      V^1 \pder{U^2}{X^1}
      +
      V^2 \pder{U^2}{X^2}
      +
      V^3 \pder{U^2}{X^3}
      +
      \frac{V^1 U^2}{X^1}
      +
      \frac{V^2 U^1}{X^1}
   \right),

.. math::

   \vec{E}_3
   \left(
      \pder{U^3}{t}
      +
      V^1 \pder{U^3}{X^1}
      +
      V^2 \pder{U^3}{X^2}
      +
      V^3 \pder{U^3}{X^3}
   \right).

=======================
Rectilinear coordinates
=======================

The sum of

.. math::

   \vec{E}_1
   \left(
      \pder{U^1}{t}
      +
      V^1 \pder{U^1}{X^1}
      +
      V^2 \pder{U^1}{X^2}
      +
      V^3 \pder{U^1}{X^3}
      +
      V^1 U^1 \frac{1}{H_1} \pder{H_1}{X^1}
   \right),

.. math::

   \vec{E}_2
   \left(
      \pder{U^2}{t}
      +
      V^1 \pder{U^2}{X^1}
      +
      V^2 \pder{U^2}{X^2}
      +
      V^3 \pder{U^2}{X^3}
      +
      V^2 U^2 \frac{1}{H_2} \pder{H_2}{X^2}
   \right),

.. math::

   \vec{E}_3
   \left(
      \pder{U^3}{t}
      +
      V^1 \pder{U^3}{X^1}
      +
      V^2 \pder{U^3}{X^2}
      +
      V^3 \pder{U^3}{X^3}
      +
      V^3 U^3 \frac{1}{H_3} \pder{H_3}{X^3}
   \right).

===========
Application
===========

The sum of

.. math::

   \vec{E}_1
   \left(
      \pder{U^1}{t}
      +
      V^1 \pder{U^1}{X^1}
      +
      V^2 \pder{U^1}{X^2}
      +
      V^3 \pder{U^1}{X^3}
      +
      V^1 U^1 \frac{1}{H_1} \pder{H_1}{X^1}
      -
      V^2 U^2 \frac{H_2}{H_1 H_1} \pder{H_2}{X^1}
   \right),

.. math::

   \vec{E}_2
   \left(
      \pder{U^2}{t}
      +
      V^1 \pder{U^2}{X^1}
      +
      V^2 \pder{U^2}{X^2}
      +
      V^3 \pder{U^2}{X^3}
      +
      V^1 U^2 \frac{1}{H_2} \pder{H_2}{X^1}
      +
      V^2 U^1 \frac{1}{H_2} \pder{H_2}{X^1}
   \right),

.. math::

   \vec{E}_3
   \left(
      \pder{U^3}{t}
      +
      V^1 \pder{U^3}{X^1}
      +
      V^2 \pder{U^3}{X^2}
      +
      V^3 \pder{U^3}{X^3}
   \right).

