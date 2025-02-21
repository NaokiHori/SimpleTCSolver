################
Normalised basis
################

***********
Description
***********

Assigning

.. math::

   \vec{u}
   =
   \sum_i
   \hat{U}^i
   \vec{\hat{E}}_i

yields

.. math::

   \pder{}{t}
   \left(
      \sum_i
      \hat{U}^i
      \vec{\hat{E}}_i
   \right)
   +
   \sum_j
   \frac{\hat{V}^j}{H_j}
   \pder{}{X^j}
   \left(
      \sum_i
      \hat{U}^i
      \vec{\hat{E}}_i
   \right)
   =
   \sum_i
   \pder{\hat{U}^i}{t}
   \vec{\hat{E}}_i
   +
   \sum_{ij}
   \frac{\hat{V}^j}{H_j}
   \pder{\hat{U}^i}{X^j}
   \vec{\hat{E}}_i
   +
   \sum_{ij}
   \frac{\hat{V}^j}{H_j}
   \hat{U}^i
   \pder{\vec{\hat{E}}_i}{X^j}.

By using :ref:`the relation <dehatdx>`:

.. math::

   \dehatdx{i}{j}{k},

the last term leads to the sum of

.. math::

   \sum_{ij}
   \frac{\hat{V}^j}{H_j}
   \hat{U}^i
   \frac{1}{H_i}
   \pder{H_j}{X^i}
   \vec{\hat{E}}_j
   =
   \sum_{ij}
   \frac{\hat{V}^i}{H_i}
   \hat{U}^j
   \frac{1}{H_j}
   \pder{H_i}{X^j}
   \vec{\hat{E}}_i

and

.. math::

   -
   \sum_{ijk}
   \frac{\hat{V}^j}{H_j}
   \hat{U}^i
   \delta_{ij}
   \frac{1}{H_k}
   \pder{H_i}{X^k}
   \vec{\hat{E}}_k
   =
   -
   \sum_{ij}
   \frac{\hat{V}^j}{H_j}
   \hat{U}^j
   \frac{1}{H_i}
   \pder{H_j}{X^i}
   \vec{\hat{E}}_i.

The explicit form is the sum of

.. math::

   \newcommand{\tmp}[1]{
      \pder{\hat{U}^{#1}}{t}
   }
   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
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
      \frac{\hat{V}^{#2}}{H_{#2}}
      \pder{\hat{U}^{#1}}{X^{#2}}
   }
   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
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
      \frac{\hat{V}^{#1}}{H_{#1}}
      \hat{U}^{#2}
      \frac{1}{H_{#2}}
      \pder{H_{#1}}{X^{#2}}
   }
   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
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
      \frac{\hat{V}^{#2}}{H_{#2}}
      \hat{U}^{#2}
      \frac{1}{H_{#1}}
      \pder{H_{#2}}{X^{#1}}
   }
   \begin{pmatrix}
      \vec{\hat{E}}_1
      &
      \vec{\hat{E}}_2
      &
      \vec{\hat{E}}_3
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

The sum of

.. math::

   \vec{\hat{E}}_1
   \left(
      \pder{\hat{U}^1}{t}
      +
      \hat{V}^1
      \pder{\hat{U}^1}{X^1}
      +
      \frac{\hat{V}^2}{X^1}
      \pder{\hat{U}^1}{X^2}
      +
      \hat{V}^3
      \pder{\hat{U}^1}{X^3}
      -
      \frac{\hat{V}^2}{X^1}
      \hat{U}^2
   \right),

.. math::

   \vec{\hat{E}}_2
   \left(
      \pder{\hat{U}^2}{t}
      +
      \hat{V}^1
      \pder{\hat{U}^2}{X^1}
      +
      \frac{\hat{V}^2}{X^1}
      \pder{\hat{U}^2}{X^2}
      +
      \hat{V}^3
      \pder{\hat{U}^2}{X^3}
      +
      \frac{\hat{V}^2}{X^1}
      \hat{U}^1
   \right),

.. math::

   \vec{\hat{E}}_3
   \left(
      \pder{\hat{U}^3}{t}
      +
      \hat{V}^1
      \pder{\hat{U}^3}{X^1}
      +
      \frac{\hat{V}^2}{X^1}
      \pder{\hat{U}^3}{X^2}
      +
      \hat{V}^3
      \pder{\hat{U}^3}{X^3}
   \right).

