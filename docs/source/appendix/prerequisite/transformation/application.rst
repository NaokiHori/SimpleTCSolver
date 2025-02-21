###########
Application
###########

I consider the cylindrical coordinates

.. math::

   & \xi^1 \equiv \sqrt{x^1 x^1 + x^2 x^2},

   & \xi^2 \equiv \tan^{-1} \left( \frac{x^2}{x^1} \right),

   & \xi^3 \equiv x^3,

which are further stretched in a rectilinear way:

.. math::

   & X^1 \equiv f_1 \left( \xi^1 \right),

   & X^2 \equiv f_2 \left( \xi^2 \right),

   & X^3 \equiv f_3 \left( \xi^3 \right).

Note that :math:`f_1` is an arbitrary function and can be non-linear, while the other two functions :math:`f_2, f_3` are linear.

Adopting the chain rule:

.. math::

   \pder{x^i}{X^j}
   =
   \sum_k
   \pder{x^i}{\xi^k}
   \pder{\xi^k}{X^j}

yields the transformation matrix:

.. math::

   \begin{pmatrix}
      \pder{x^1}{X^1} & \pder{x^1}{X^2} & \pder{x^1}{X^3} \\
      \pder{x^2}{X^1} & \pder{x^2}{X^2} & \pder{x^2}{X^3} \\
      \pder{x^3}{X^1} & \pder{x^3}{X^2} & \pder{x^3}{X^3} \\
   \end{pmatrix}
   &
   =
   \begin{pmatrix}
      \pder{x^1}{\xi^1} & \pder{x^1}{\xi^2} & \pder{x^1}{\xi^3} \\
      \pder{x^2}{\xi^1} & \pder{x^2}{\xi^2} & \pder{x^2}{\xi^3} \\
      \pder{x^3}{\xi^1} & \pder{x^3}{\xi^2} & \pder{x^3}{\xi^3} \\
   \end{pmatrix}
   \begin{pmatrix}
      \pder{\xi^1}{X^1} & \pder{\xi^1}{X^2} & \pder{\xi^1}{X^3} \\
      \pder{\xi^2}{X^1} & \pder{\xi^2}{X^2} & \pder{\xi^2}{X^3} \\
      \pder{\xi^3}{X^1} & \pder{\xi^3}{X^2} & \pder{\xi^3}{X^3} \\
   \end{pmatrix}
   \\
   &
   =
   \begin{pmatrix}
      \cos \xi^2 & - \xi^1 \sin \xi^2 & 0 \\
      \sin \xi^2 &   \xi^1 \cos \xi^2 & 0 \\
               0 &                  0 & 1 \\
   \end{pmatrix}
   \begin{pmatrix}
      \pder{\xi^1}{X^1} &                 0 &                 0 \\
                      0 & \pder{\xi^2}{X^2} &                 0 \\
                      0 &                 0 & \pder{\xi^3}{X^3} \\
   \end{pmatrix}
   \\
   &
   =
   \begin{pmatrix}
      \pder{\xi^1}{X^1} \cos \xi^2 & - \pder{\xi^2}{X^2} \xi^1 \sin \xi^2 &                 0 \\
      \pder{\xi^1}{X^1} \sin \xi^2 &   \pder{\xi^2}{X^2} \xi^1 \cos \xi^2 &                 0 \\
                                 0 &                                    0 & \pder{\xi^3}{X^3} \\
   \end{pmatrix}.

Note that, because of the linearity of :math:`f_2` and :math:`f_3`, :math:`\pder{\xi^2}{X^2}` and :math:`\pder{\xi^3}{X^3}` (specifically :math:`\Delta \theta` and :math:`\Delta z`, respectively) are both constant in space.
:math:`\pder{\xi^1}{X^1}` (:math:`\Delta r`), on the other hand, can vary in :math:`X^1`.

The relation of the basis vectors are

.. math::

   &
   \vec{E}_1
   =
   \vec{e}_1
   \pder{\xi^1}{X^1} \cos \xi^2
   +
   \vec{e}_2
   \pder{\xi^1}{X^1} \sin \xi^2,

   &
   \vec{E}_2
   =
   -
   \vec{e}_1
   \pder{\xi^2}{X^2} \xi^1 \sin \xi^2
   +
   \vec{e}_2
   \pder{\xi^2}{X^2} \xi^1 \cos \xi^2,

   &
   \vec{E}_3
   =
   \vec{e}_3
   \pder{\xi^3}{X^3}.

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
   \pder{\xi^1}{X^1}
   \left(
      =
      \Delta r
   \right),

   &
   H_2
   =
   \sqrt{
      \vec{E}_2
      \cdot
      \vec{E}_2
   }
   =
   \pder{\xi^2}{X^2}
   \xi^1
   \left(
      =
      r \Delta \theta
   \right),

   &
   H_3
   =
   \sqrt{
      \vec{E}_3
      \cdot
      \vec{E}_3
   }
   =
   \pder{\xi^3}{X^3}
   \left(
      =
      \Delta z
   \right),

and thus the Jacobian determinant leads to

.. math::

   J
   =
   \xi^1
   \pder{\xi^1}{X^1}
   \pder{\xi^2}{X^2}
   \pder{\xi^3}{X^3}
   \left(
      =
      r
      \Delta r
      \Delta \theta
      \Delta z
   \right).

The spatial derivatives of :math:`H_i` are

.. math::

   \begin{pmatrix}
      \pder{H_1}{X^1} & \pder{H_2}{X^1} & \pder{H_3}{X^1} \\
      \pder{H_1}{X^2} & \pder{H_2}{X^2} & \pder{H_3}{X^2} \\
      \pder{H_1}{X^3} & \pder{H_2}{X^3} & \pder{H_3}{X^3} \\
   \end{pmatrix}
   =
   \begin{pmatrix}
      \pder{H_1}{X^1} & \pder{H_2}{X^1} & 0 \\
      0 & \pder{H_2}{X^2} & 0 \\
      0 & 0 & \pder{H_3}{X^3} \\
   \end{pmatrix},

which can be further simplified as

.. math::

   \begin{pmatrix}
      \pder{H_1}{X^1} & \pder{H_2}{X^1} & 0 \\
      0 & 0 & 0 \\
      0 & 0 & 0 \\
   \end{pmatrix},

by utilising the linearity of the mapping functions :math:`f_2, f_3`.

