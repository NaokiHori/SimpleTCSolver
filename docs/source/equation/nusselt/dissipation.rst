###########
Dissipation
###########

I consider to integrate the equation of the kinetic energy in the whole volume, giving

.. math::

   \left[
      \int_{\vt}
      \mu
      \left(
         \ut \vr \der{\ut}{\vr}
      \right)
      d\vt
   \right]_{r_i}^{r_o}
   =
   \int_V
   \epsilon
   dV

in a statistical sense, where

.. math::

   \epsilon
   \equiv
   \mu
   \left(
        L_{\vr \vr}^2
      + L_{\vr \vt}^2
      + L_{\vt \vr}^2
      + L_{\vt \vt}^2
   \right).

The left-hand-side term leads to

.. math::

   \ut_o
   \int_{S_o}
   \mu
   \der{\ut}{\vr}
   dS_o
   -
   \ut_i
   \int_{S_i}
   \mu
   \der{\ut}{\vr}
   dS_i,

which is the kinetic energy flux through the walls.

Now I revisit the definition of the Nusselt number

.. math::

   Nu
   T_{lam}
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

   Nu
   T_{lam}
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

or equivalently

.. math::

   \frac{\ut_i}{\vr_i}
   Nu
   T_{lam}
   =
   -
   \ut_i
   \int_{S_i}
   \mu
   \der{\ut}{\vr}
   dS_i
   +
   \ut_i
   \frac{\ut_i}{\vr_i}
   \int_{S_i}
   \mu
   dS_i

and

.. math::

   \frac{\ut_o}{\vr_o}
   Nu
   T_{lam}
   =
   -
   \ut_o
   \int_{S_o}
   \mu
   \der{\ut}{\vr}
   dS_o
   +
   \ut_o
   \frac{\ut_o}{\vr_o}
   \int_{S_o}
   \mu
   dS_o.

Subtracting the second equation from the first one results in

.. math::

   Nu T_{lam} \left(
      \frac{\ut_i}{\vr_i}
      -
      \frac{\ut_o}{\vr_o}
   \right)
   =
   \int_V \epsilon dV
   +
   \mu
   \left(
      \ut_i^2
      -
      \ut_o^2
   \right)
   \int d\vt,

which is the description of the Nusselt number based on the energy dissipation.

