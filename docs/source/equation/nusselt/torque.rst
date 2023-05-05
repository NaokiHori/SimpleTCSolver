######
Torque
######

I consider to integrate the momentum equation in the azimuthal direction after :math:`r` is multiplied.

Advective terms are dropped for simplicity, and the diffusive terms yield

.. math::

   \int \int
   \mu
   \left\{
      \der{}{\vr} \left(
         \vr \der{\ut}{\vr}
      \right)
      +
      \der{}{\vt} \left(
         \frac{1}{\vr} \der{\ut}{\vt}
      \right)
      -
      \frac{\ut}{\vr}
      +
      2 \left( \frac{1}{\vr} \der{\ur}{\vt} \right)
   \right\}
   \vr d\vr d\vt.

Since

.. math::

      \der{}{\vr} \left(
         \vr \der{\ut}{\vr}
      \right)
      -
      \frac{\ut}{\vr}
      \equiv
      \frac{1}{\vr}
      \der{}{\vr} \left\{
         \vr^3 \der{}{\vr} \left(
            \frac{\ut}{\vr}
         \right)
      \right\},

the above equation results in

.. math::

   \left[
      \int
      \vr
      \mu
      \vr \der{}{\vr} \left(
         \frac{\ut}{\vr}
      \right)
      \vr d\vt
   \right]_{r_i}^{r_o}
   =
   \left[
      \int
      \vr
      \times
      \mu
      \left(
         \der{\ut}{\vr}
         -
         \frac{\ut}{\vr}
      \right)
      \times
      \vr d\vt
   \right]_{r_i}^{r_o},

which is equal to zero in a statistical sense, indicating that the torque evaluated on the inner and the outer cylinder walls match.
Also note that I used an identity

.. math::

   \vr
   \der{}{\vr}
   \left(
      \frac{\ut}{\vr}
   \right)
   \equiv
   \der{\ut}{\vr}
   -
   \frac{\ut}{\vr}

for later convenience.

Finally I notice the torque measured on the inner and the outer cylinders are given by

.. math::

   T \left( r = r_i \right)
   \equiv
   -
   \left[
      \int
      \vr
      \times
      \mu
      \left(
         \der{\ut}{\vr}
         -
         \frac{\ut}{\vr}
      \right)
      \times
      \vr d\vt
   \right]_{r = r_i}
   =
   -
   \int_{S_i}
   \vr_i
   \mu
   \left(
      \der{\ut}{\vr}
      -
      \frac{\ut}{\vr}
   \right)
   dS_i

and

.. math::

   T \left( r = r_o \right)
   \equiv
   -
   \left[
      \int
      \vr
      \times
      \mu
      \left(
         \der{\ut}{\vr}
         -
         \frac{\ut}{\vr}
      \right)
      \times
      \vr d\vt
   \right]_{r = r_o}
   =
   -
   \int_{S_o}
   \vr_o
   \mu
   \left(
      \der{\ut}{\vr}
      -
      \frac{\ut}{\vr}
   \right)
   dS_o,

respectively.
Note that the negative sign comes from the Newton's third law.

To find the reference value, I consider the laminar state.
Since I have

.. math::

   \vr^3 \der{}{\vr} \left(
      \frac{\ut}{\vr}
   \right)
   =
   - 2 A,

I have

.. math::

   \ut
   =
   \frac{A}{\vr}
   +
   B \vr,

where :math:`A` and :math:`B` are constants of integration determined by the boundary conditions:

.. math::

   A
   =
   \frac{
      \vr_i \vr_o \left(
         \ut_i \vr_o
         -
         \ut_o \vr_i
      \right)
   }{
      \vr_o^2 - \vr_i^2
   },

.. math::

   B
   =
   \frac{
      \ut_o \vr_o
      -
      \ut_i \vr_i
   }{
      \vr_o^2 - \vr_i^2
   }.

As a result, the laminar value :math:`T_{lam}` is

.. math::

   T_{lam}
   =
   \int_{S_i}
   \mu
   \frac{2A}{\vr_i}
   dS_i
   =
   \int_{S_o}
   \mu
   \frac{2A}{\vr_o}
   dS_o
   =
   2 \mu A \int_{\vt} d\vt.

Thus the Nusselt number, which is defined as the torque normalised by the laminar value, is

.. math::

   Nu T_{lam}
   \equiv
   T \left( r_i \right)
   \equiv
   T \left( r_o \right).

