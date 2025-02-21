#######################
Rectilinear coordinates
#######################

.. note::

   Normalised representations will not be discussed since they coincide with the Cartesian one.

The rectilinear coordinates are defined as

.. math::

   & X^1 \equiv f_1 \left( x^1 \right),

   & X^2 \equiv f_2 \left( x^2 \right),

   & X^3 \equiv f_3 \left( x^3 \right).

:math:`f_i` are not explicitly discussed here, which are to be determined when the computational grid is designed.
Without loss of generality, I assume that the mapping functions :math:`f_i` monotonically increase.

The transformation matrix is

.. math::

   \begin{pmatrix}
      \pder{x^1}{X^1} & \pder{x^1}{X^2} & \pder{x^1}{X^3} \\
      \pder{x^2}{X^1} & \pder{x^2}{X^2} & \pder{x^2}{X^3} \\
      \pder{x^3}{X^1} & \pder{x^3}{X^2} & \pder{x^3}{X^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      \pder{x^1}{X^1} &               0 &               0 \\
                    0 & \pder{x^2}{X^2} &               0 \\
                    0 &               0 & \pder{x^3}{X^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      H_1 &   0 &   0 \\
        0 & H_2 &   0 \\
        0 &   0 & H_3 \\
   \end{pmatrix},

or

.. math::

   \begin{pmatrix}
      \pder{X^1}{x^1} & \pder{X^1}{x^2} & \pder{X^1}{x^3} \\
      \pder{X^2}{x^1} & \pder{X^2}{x^2} & \pder{X^2}{x^3} \\
      \pder{X^3}{x^1} & \pder{X^3}{x^2} & \pder{X^3}{x^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      \pder{X^1}{x^1} &               0 &               0 \\
                    0 & \pder{X^2}{x^2} &               0 \\
                    0 &               0 & \pder{X^3}{x^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      \frac{1}{H_1} &             0 &             0 \\
                  0 & \frac{1}{H_2} &             0 \\
                  0 &             0 & \frac{1}{H_3} \\
   \end{pmatrix}.

As a consequence, the basis vectors are related by

.. math::

   &
   \vec{E}_1
   =
   \pder{x^1}{X^1}
   \vec{e}_1
   =
   H_1
   \vec{e}_1,

   &
   \vec{E}_2
   =
   \pder{x^2}{X^2}
   \vec{e}_2
   =
   H_2
   \vec{e}_2,

   &
   \vec{E}_3
   =
   \pder{x^3}{X^3}
   \vec{e}_3
   =
   H_3
   \vec{e}_3,

or

.. math::

   &
   \vec{e}_1
   =
   \pder{X^1}{x^1}
   \vec{E}_1
   =
   \frac{1}{H_1}
   \vec{E}_1,

   &
   \vec{e}_2
   =
   \pder{X^2}{x^2}
   \vec{E}_2
   =
   \frac{1}{H_2}
   \vec{E}_2,

   &
   \vec{e}_3
   =
   \pder{X^3}{x^3}
   \vec{E}_3
   =
   \frac{1}{H_3}
   \vec{E}_3.

Although already introduced, the scale factors are

.. math::

   H_i
   =
   \sqrt{
      \vec{E}_i
      \cdot
      \vec{E}_i
   }
   =
   \pder{x^i}{X^i},

and their spatial derivatives are

.. math::

   \begin{pmatrix}
      \pder{H_1}{X^1} & \pder{H_2}{X^1} & \pder{H_3}{X^1} \\
      \pder{H_1}{X^2} & \pder{H_2}{X^2} & \pder{H_3}{X^2} \\
      \pder{H_1}{X^3} & \pder{H_2}{X^3} & \pder{H_3}{X^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      \pder{H_1}{X^1} & 0 & 0 \\
      0 & \pder{H_2}{X^2} & 0 \\
      0 & 0 & \pder{H_3}{X^3} \\
   \end{pmatrix}.

The Jacobian determinant is

.. math::

   J
   =
   H_1
   H_2
   H_3
   =
   \pder{x^1}{X^1}
   \pder{x^2}{X^2}
   \pder{x^3}{X^3}.

Although not limited, I consider the velocity vector

.. math::

   \vec{u}
   \equiv
   \tder{\vec{x}}{t}
   =
   \sum_i
   u^i
   \vec{e}_i

as an example to see how the components are related.
Then the contravariant components are

.. math::

   &
   U^1
   =
   \pder{X^1}{x^1}
   u^1
   =
   \frac{u^1}{H_1},

   &
   U^2
   =
   \pder{X^2}{x^2}
   u^2
   =
   \frac{u^2}{H_2},

   &
   U^3
   =
   \pder{X^3}{x^3}
   u^3
   =
   \frac{u^3}{H_3}.

