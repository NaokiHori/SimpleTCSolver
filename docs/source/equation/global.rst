###############
Global relation
###############

***************
Advective terms
***************

I consider to add the incompressibility constraint multiplied by the squared velocity in the direction

.. math::

   k_i \der{u_j}{x_j}

to the advective terms, giving

.. math::

   &
   \vr \ur \frac{1}{\vr} \der{\kr}{\vr}
   +
   \ut \frac{1}{\vr} \der{\kr}{\vt}
   -
   \frac{\ur \ut^2}{\vr}
   +
   \uz \der{\kr}{\vz}
   +
   \kr \frac{1}{\vr} \der{}{\vr} \left( \vr \ur \right)
   +
   \kr \frac{1}{\vr} \der{\ut}{\vt}
   +
   \kr \der{\uz}{\vz} \\
   =
   &
   \frac{1}{\vr} \der{}{\vr} \left( \vr \ur \kr \right)
   +
   \frac{1}{\vr} \der{}{\vt} \left( \ut \kr \right)
   -
   \frac{\ur \ut^2}{\vr}
   +
   \der{}{\vz} \left( \uz \kr \right)

from the radial,

.. math::

   &
   \vr \ur \frac{1}{\vr} \der{\kt}{\vr}
   +
   \ut \frac{1}{\vr} \der{\kt}{\vt}
   +
   \frac{\ur \ut^2}{\vr}
   +
   \uz \der{\kt}{\vz}
   +
   \kt \frac{1}{\vr} \der{}{\vr} \left( \vr \ur \right)
   +
   \kt \frac{1}{\vr} \der{\ut}{\vt}
   +
   \kt \der{\uz}{\vz} \\
   =
   &
   \frac{1}{\vr} \der{}{\vr} \left( \vr \ur \kt \right)
   +
   \frac{1}{\vr} \der{}{\vt} \left( \ut \kt \right)
   +
   \frac{\ur \ut^2}{\vr}
   +
   \der{}{\vz} \left( \uz \kt \right)

from the azimuthal, and

.. math::

   &
   \vr \ur \frac{1}{\vr} \der{\kz}{\vr}
   +
   \ut \frac{1}{\vr} \der{\kz}{\vt}
   +
   \uz \der{\kz}{\vz}
   +
   \kz \frac{1}{\vr} \der{}{\vr} \left( \vr \ur \right)
   +
   \kz \frac{1}{\vr} \der{\ut}{\vt}
   +
   \kz \der{\uz}{\vz} \\
   =
   &
   \frac{1}{\vr} \der{}{\vr} \left( \vr \ur \kz \right)
   +
   \frac{1}{\vr} \der{}{\vt} \left( \ut \kz \right)
   +
   \der{}{\vz} \left( \uz \kz \right)

from the axial directions, respectively.

After integrated in the volume

.. math::

   \int \int \int q \vr d\vr d\vt d\vz,

all conservative terms vanish because of the impermeable and the periodic boundary conditions.
Also the centrifugal and the Coriolis contributions cancel out to each other, and as a result I confirm that the advective terms do not modify the total kinetic energy (as expected).

***************
Diffusive terms
***************

After integrated in the volume

.. math::

   \int \int \int q \vr d\vr d\vt d\vz,

all conservative terms vanish because of the impermeable and the periodic boundary conditions, except

.. math::

   \int \int \int
   \mu \frac{1}{\vr} \der{}{\vr} \left(
      \vr \der{\kt}{\vr}
   \right)
   \vr d\vr d\vt d\vz
   +
   \int \int \int
   \mu \frac{1}{\vr} \der{}{\vr} \left(
      \vr \der{\kz}{\vr}
   \right)
   \vr d\vr d\vt d\vz

which yields

.. math::

   \int \int
   \mu \left(
      \der{\kt}{\vr}
   \right)
   \vr d\vt d\vz
   +
   \int \int
   \mu \left(
      \der{\kz}{\vr}
   \right)
   \vr d\vt d\vz

or

.. math::

   \ut_o
   \int_{S_o}
   \mu \der{\ut}{\vr}
   dS_o
   -
   \ut_i
   \int_{S_i}
   \mu \der{\ut}{\vr}
   dS_i
   +
   \uz_o
   \int_{S_o}
   \mu \der{\uz}{\vr}
   dS_o
   -
   \uz_i
   \int_{S_i}
   \mu \der{\uz}{\vr}
   dS_i

where the integrals in the last equation are the surface integrals on the inner and the outer cylinder walls.
Thus the velocities are the given boundary conditions, while the velocity gradients are the system response.
Usually I assume that the walls do not move in the axial direction in this project, and a result the :math:`\uz` contributions are zero.

The non-conservative terms yield

.. math::

   -
   \int \int \int
   \epsilon
   \vr d\vr d\vt d\vz,

where

.. math::

   \epsilon
   \equiv
   &
   \mu \left( \der{\ur}{\vr} \right)^2
   +
   \mu \left( \frac{1}{\vr} \der{\ur}{\vt} \right)^2
   +
   \mu \left( \der{\ur}{\vz} \right)^2
   +
   \mu \left( \frac{\ur}{\vr} \right)^2
   +
   \mu \left( 2 \frac{\ur}{\vr} \frac{1}{\vr} \der{\ut}{\vt} \right) \\
   +
   &
   \mu \left( \der{\ut}{\vr} \right)^2
   +
   \mu \left( \frac{1}{\vr} \der{\ut}{\vt} \right)^2
   +
   \mu \left( \der{\ut}{\vz} \right)^2
   +
   \mu \left( \frac{\ut}{\vr} \right)^2
   -
   \mu \left( 2 \frac{\ut}{\vr} \frac{1}{\vr} \der{\ur}{\vt} \right) \\
   +
   &
   \mu \left( \der{\uz}{\vr} \right)^2
   +
   \mu \left( \frac{1}{\vr} \der{\uz}{\vt} \right)^2
   +
   \mu \left( \der{\uz}{\vz} \right)^2,

which is

.. math::

   &
   \mu \left( \der{\ur}{\vr} \right)^2
   +
   \mu \left( \frac{1}{\vr} \der{\ur}{\vt} - \frac{\ut}{\vr} \right)^2
   +
   \mu \left( \der{\ur}{\vz} \right)^2 \\
   +
   &
   \mu \left( \der{\ut}{\vr} \right)^2
   +
   \mu \left( \frac{1}{\vr} \der{\ut}{\vt} + \frac{\ur}{\vr} \right)^2
   +
   \mu \left( \der{\ut}{\vz} \right)^2 \\
   +
   &
   \mu \left( \der{\uz}{\vr} \right)^2
   +
   \mu \left( \frac{1}{\vr} \der{\uz}{\vt} \right)^2
   +
   \mu \left( \der{\uz}{\vz} \right)^2,

or in short

.. math::

   \mu L_{ij} L_{ij}.

In summary, for normal Taylor-Couette flows adopted in this project where walls move only in the azimuthal direction, I have

.. math::

   \ut_o
   \int_{S_o}
   \mu \der{\ut}{\vr}
   dS_o
   -
   \ut_i
   \int_{S_i}
   \mu \der{\ut}{\vr}
   dS_i
   =
   \int_{V}
   \mu L_{ij} L_{ij}
   dV,

indicating that the kinetic energy injection through the diffusive effect on the walls is balanced by the energy dissipation in the flow in a statistical sense.

