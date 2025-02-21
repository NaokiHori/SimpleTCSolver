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
      \pder{U^1}{X^1}
      &
      \pder{U^2}{X^1}
      &
      \pder{U^3}{X^1}
      \\
      \frac{1}{X^1 X^1}
      \pder{U^1}{X^2}
      &
      \frac{1}{X^1 X^1}
      \pder{U^2}{X^2}
      &
      \frac{1}{X^1 X^1}
      \pder{U^3}{X^2}
      \\
      \pder{U^1}{X^3}
      &
      \pder{U^2}{X^3}
      &
      \pder{U^3}{X^3}
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
      \frac{1}{X^1 X^1}
      \frac{U^1}{X^1}
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
      \frac{U^2}{X^1}
      &
      0
      \\
      -
      \frac{U^2}{X^1}
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
      \pder{U^1}{X^1}
      &
      \frac{1}{2}
      \pder{U^2}{X^1}
      +
      \frac{1}{2}
      \frac{1}{X^1 X^1}
      \pder{U^1}{X^2}
      &
      \frac{1}{2}
      \pder{U^3}{X^1}
      +
      \frac{1}{2}
      \pder{U^1}{X^3}
      \\
      sym.
      &
      \frac{1}{X^1 X^1}
      \pder{U^2}{X^2}
      +
      \frac{1}{X^1 X^1}
      \frac{U^1}{X^1}
      &
      \frac{1}{2}
      \frac{1}{X^1 X^1}
      \pder{U^3}{X^2}
      +
      \frac{1}{2}
      \pder{U^2}{X^3}
      \\
      sym.
      &
      sym.
      &
      \pder{U^3}{X^3}
      \\
   \end{pmatrix}
   \begin{pmatrix}
      \vec{E}_1
      \\
      \vec{E}_2
      \\
      \vec{E}_3
   \end{pmatrix}.

