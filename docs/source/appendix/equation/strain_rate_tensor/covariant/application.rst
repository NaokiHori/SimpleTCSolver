The velocity-gradient tensor is the sum of

.. math::

   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      \tmpd{1}{1}
      &
      \tmpd{2}{1}
      &
      \tmpd{3}{1}
      \\
      \tmpd{1}{2}
      &
      \tmpd{2}{2}
      &
      \tmpd{3}{2}
      \\
      \tmpd{1}{3}
      &
      \tmpd{2}{3}
      &
      \tmpd{3}{3}
      \\
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
      0
      &
      0
      &
      0
      \\
      0
      &
      \frac{1}{H_2 H_2} \frac{1}{H_2} \pder{H_2}{X^1} U^1
      &
      0
      \\
      0
      &
      0
      &
      0
      \\
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
      0
      &
      \tmpa{2}{1}
      &
      0
      \\
      -
      \tmpa{2}{1}
      &
      0
      &
      0
      \\
      0
      &
      0
      &
      0
      \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{E}_1
      \\
      \vec{E}_2
      \\
      \vec{E}_3
   \end{pmatrix}.

The strain-rate tensor is

.. math::

   \begin{pmatrix}
      \vec{E}_1
      &
      \vec{E}_2
      &
      \vec{E}_3
   \end{pmatrix}
   \begin{pmatrix}
      \frac{1}{H_1 H_1}
      \pder{U^1}{X^1}
      +
      \frac{1}{H_1 H_1}
      \frac{1}{H_1}
      \pder{H_1}{X^1}
      U^1
      &
      \frac{1}{2}
      \frac{1}{H_1 H_1}
      \pder{U^2}{X^1}
      +
      \frac{1}{2}
      \frac{1}{H_1 H_1}
      \pder{U^2}{X^1}
      &
      \frac{1}{2}
      \frac{1}{H_1 H_1}
      \pder{U^3}{X^1}
      +
      \frac{1}{2}
      \frac{1}{H_3 H_3}
      \pder{U^1}{X^3}
      \\
      sym.
      &
      \frac{1}{H_2 H_2}
      \pder{U^2}{X^2}
      +
      \frac{1}{H_2 H_2}
      \frac{1}{H_2}
      \pder{H_2}{X^1}
      U^1
      &
      \frac{1}{2}
      \frac{1}{H_2 H_2}
      \pder{U^3}{X^2}
      +
      \frac{1}{2}
      \frac{1}{H_3 H_3}
      \pder{U^2}{X^3}
      \\
      sym.
      &
      sym.
      &
      \frac{1}{H_3 H_3}
      \pder{U^3}{X^3}
   \end{pmatrix}
   \begin{pmatrix}
      \vec{E}_1
      \\
      \vec{E}_2
      \\
      \vec{E}_3
   \end{pmatrix}.

