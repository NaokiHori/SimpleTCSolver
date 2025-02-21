################################
Energy diffusion and dissipation
################################

The divergence of a second-order tensor is given by the sum of

.. math::

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
   \newcommand{\tmpc}[2]{
      \frac{1}{J}
      \pder{}{X^{#2}}
      \left(
         J
         S^{#1 #2}
      \right)
   }
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

I consider the inner product of this vector and the velocity vector:

.. math::

   \vec{u}
   =
   \sum_i
   U^i
   \vec{E}_i.

The first part yields

.. math::

   \frac{1}{J}
   \pder{}{X^1}
   \left(
      J
      U^1
      S^{1 1}
   \right)
   +
   \frac{1}{J}
   \pder{}{X^2}
   \left(
      J
      U^1
      S^{1 2}
   \right)
   +
   \frac{1}{J}
   \pder{}{X^3}
   \left(
      J
      U^1
      S^{1 3}
   \right)
   -
   \pder{U^1}{X^1}
   S^{1 1}
   -
   \pder{U^1}{X^2}
   S^{1 2}
   -
   \pder{U^1}{X^3}
   S^{1 3}

.. math::

   \frac{1}{J}
   \pder{}{X^1}
   \left(
      J
      U^2
      S^{2 1}
   \right)
   +
   \frac{1}{J}
   \pder{}{X^2}
   \left(
      J
      U^2
      S^{2 2}
   \right)
   +
   \frac{1}{J}
   \pder{}{X^3}
   \left(
      J
      U^2
      S^{2 3}
   \right)
   -
   \pder{U^2}{X^1}
   S^{2 1}
   -
   \pder{U^2}{X^2}
   S^{2 2}
   -
   \pder{U^2}{X^3}
   S^{2 3}

.. math::

   \frac{1}{J}
   \pder{}{X^1}
   \left(
      J
      U^3
      S^{3 1}
   \right)
   +
   \frac{1}{J}
   \pder{}{X^2}
   \left(
      J
      U^3
      S^{3 2}
   \right)
   +
   \frac{1}{J}
   \pder{}{X^3}
   \left(
      J
      U^3
      S^{3 3}
   \right)
   -
   \pder{U^3}{X^1}
   S^{3 1}
   -
   \pder{U^3}{X^2}
   S^{3 2}
   -
   \pder{U^3}{X^3}
   S^{3 3}

In each term, the first-three terms describe the diffusion of the energy, while the last-three terms are the dissipation of the energy.

The rest of the vector yields the sum of

.. math::

   S^{1 1}
   \left(
      \frac{2}{H_1}
      \pder{H_1}{X^1}
      U^1
      -
      \frac{H_1}{H_1 H_1}
      \pder{H_1}{X^1}
      U^1
      -
      \frac{H_1}{H_2 H_2}
      \pder{H_1}{X^2}
      U^2
      -
      \frac{H_1}{H_3 H_3}
      \pder{H_1}{X^3}
      U^3
   \right)

.. math::

   S^{1 2}
   \left(
      \frac{1}{H_1}
      \pder{H_1}{X^2}
      U^1
      +
      \frac{1}{H_2}
      \pder{H_2}{X^1}
      U^2
   \right)

.. math::

   S^{1 3}
   \left(
      \frac{1}{H_1}
      \pder{H_1}{X^3}
      U^1
      +
      \frac{1}{H_3}
      \pder{H_3}{X^1}
      U^3
   \right)

.. math::

   S^{2 1}
   \left(
      \frac{1}{H_1}
      \pder{H_1}{X^2}
      U^1
      +
      \frac{1}{H_2}
      \pder{H_2}{X^1}
      U^2
   \right)

.. math::

   S^{2 2}
   \left(
      \frac{2}{H_2}
      \pder{H_2}{X^2}
      U^2
      -
      \frac{H_2}{H_1 H_1}
      \pder{H_2}{X^1}
      U^1
      -
      \frac{H_2}{H_2 H_2}
      \pder{H_2}{X^2}
      U^2
      -
      \frac{H_2}{H_3 H_3}
      \pder{H_2}{X^3}
      U^3
   \right)

.. math::

   S^{2 3}
   \left(
      \frac{1}{H_3}
      \pder{H_3}{X^2}
      U^3
      +
      \frac{1}{H_2}
      \pder{H_2}{X^3}
      U^2
   \right)

.. math::

   S^{3 1}
   \left(
      \frac{1}{H_1}
      \pder{H_1}{X^3}
      U^1
      +
      \frac{1}{H_3}
      \pder{H_3}{X^1}
      U^3
   \right)

.. math::

   S^{3 2}
   \left(
      \frac{1}{H_2}
      \pder{H_2}{X^3}
      U^2
      +
      \frac{1}{H_3}
      \pder{H_3}{X^2}
      U^3
   \right)

.. math::

   S^{3 3}
   \left(
      \frac{2}{H_3}
      \pder{H_3}{X^3}
      U^3
      -
      \frac{H_3}{H_1 H_1}
      \pder{H_3}{X^1}
      U^1
      -
      \frac{H_3}{H_2 H_2}
      \pder{H_3}{X^2}
      U^2
      -
      \frac{H_3}{H_3 H_3}
      \pder{H_3}{X^3}
      U^3
   \right)

