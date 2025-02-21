#######################
Cylindrical coordinates
#######################

The cylindrical coordinates are defined as

.. math::

   & X^1 \equiv \sqrt{x^1 x^1 + x^2 x^2},

   & X^2 \equiv \tan^{-1} \left( \frac{x^2}{x^1} \right),

   & X^3 \equiv x^3,

or equivalently

.. math::

   & x^1 = X^1 \cos X^2,

   & x^2 = X^1 \sin X^2,

   & x^3 = X^3.

The transformation matrix is

.. math::

   \begin{pmatrix}
      \pder{x^1}{X^1} & \pder{x^1}{X^2} & \pder{x^1}{X^3} \\
      \pder{x^2}{X^1} & \pder{x^2}{X^2} & \pder{x^2}{X^3} \\
      \pder{x^3}{X^1} & \pder{x^3}{X^2} & \pder{x^3}{X^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      \cos X^2 & - X^1 \sin X^2 & 0 \\
      \sin X^2 &   X^1 \cos X^2 & 0 \\
             0 &              0 & 1 \\
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
      \cos X^2
      &
      \sin X^2
      &
      0
      \\
      - \frac{\sin X^2}{X^1}
      &
      \frac{\cos X^2}{X^1}
      &
      0
      \\
      0
      &
      0
      &
      1
      \\
   \end{pmatrix}.

As a consequence, the basis vectors are related by

.. math::

   &
   \vec{E}_1
   =
   \vec{e}_1
   \cos X^2
   +
   \vec{e}_2
   \sin X^2,

   &
   \vec{E}_2
   =
   -
   \vec{e}_1
   X^1
   \sin X^2
   +
   \vec{e}_2
   X^1
   \cos X^2,

   &
   \vec{E}_3
   =
   \vec{e}_3,

or

.. math::

   &
   \vec{e}_1
   =
   \vec{E}_1
   \cos X^2
   -
   \vec{E}_2
   \frac{\sin X^2}{X^1},

   &
   \vec{e}_2
   =
   \vec{E}_1
   \sin X^2
   +
   \vec{E}_2
   \frac{\cos X^2}{X^1},

   &
   \vec{e}_3
   =
   \vec{E}_3.

The scale factors are

.. math::

   &
   H_1
   =
   \sqrt{
      \vec{E}_1
      \cdot
      \vec{E}_1
   }
   =
   1,

   &
   H_2
   =
   \sqrt{
      \vec{E}_2
      \cdot
      \vec{E}_2
   }
   =
   X^1,

   &
   H_3
   =
   \sqrt{
      \vec{E}_3
      \cdot
      \vec{E}_3
   }
   =
   1,

and their spatial derivatives are

.. math::

   \begin{pmatrix}
      \pder{H_1}{X^1} & \pder{H_2}{X^1} & \pder{H_3}{X^1} \\
      \pder{H_1}{X^2} & \pder{H_2}{X^2} & \pder{H_3}{X^2} \\
      \pder{H_1}{X^3} & \pder{H_2}{X^3} & \pder{H_3}{X^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      0 & 1 & 0 \\
      0 & 0 & 0 \\
      0 & 0 & 0 \\
   \end{pmatrix}.

The Jacobian determinant is

.. math::

   J
   =
   H_1
   H_2
   H_3
   =
   X^1.

Thus the normalised basis vectors are

.. math::

   &
   \vec{\hat{E}}_1
   =
   \vec{e}_1
   \cos X^2
   +
   \vec{e}_2
   \sin X^2,

   &
   \vec{\hat{E}}_2
   =
   -
   \vec{e}_1
   \sin X^2
   +
   \vec{e}_2
   \cos X^2,

   &
   \vec{\hat{E}}_3
   =
   \vec{e}_3,

or inversely

.. math::

   &
   \vec{e}_1
   =
   \vec{\hat{E}}_1
   \cos X^2
   -
   \vec{\hat{E}}_2
   \sin X^2,

   &
   \vec{e}_2
   =
   \vec{\hat{E}}_1
   \sin X^2
   +
   \vec{\hat{E}}_2
   \cos X^2,

   &
   \vec{e}_3
   =
   \vec{\hat{E}}_3.

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
   u^1
   \cos X^2
   +
   u^2
   \sin X^2,

   &
   U^2
   =
   -
   u^1
   \frac{\sin X^2}{X^1}
   +
   u^2
   \frac{\cos X^2}{X^1},

   &
   U^3
   =
   u^3,

while the normalised components are

.. math::

   &
   \hat{U}^1
   =
   u^1
   \cos X^2
   +
   u^2
   \sin X^2,

   &
   \hat{U}^2
   =
   -
   u^1
   \sin X^2
   +
   u^2
   \cos X^2,

   &
   \hat{U}^3
   =
   u^3.

