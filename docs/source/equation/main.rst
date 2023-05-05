########
Equation
########

Governing equations are as follows.

.. _eq_mass:

Mass conservation

.. math::

   \frac{1}{\vr} \der{\vr \ur}{\vr}
   +
   \frac{1}{\vr} \der{\ut}{\vt}
   +
   \der{\uz}{\vz}
   =
   0.

Momentum balance

.. math::

   \der{\ur}{t}
   &
   +
   \ur \lrr
   +
   \ut \left( \lrt \right)
   +
   \uz \lrz
   =
   -
   \der{p}{\vr} \\
   &
   +
   \mu \left[
      \frac{1}{\vr} \der{}{\vr} \left(
         \vr \der{\ur}{\vr}
      \right)
      +
      \frac{1}{\vr} \der{}{\vt} \left(
         \frac{1}{\vr} \der{\ur}{\vt}
      \right)
      +
      \der{}{\vz} \left(
         \der{\ur}{\vz}
      \right)
      -
      \frac{\ur}{\vr^2}
      -
      \frac{2}{\vr^2} \der{\ut}{\vt}
   \right], \\
   \der{\ut}{t}
   &
   +
   \ur \ltr
   +
   \ut \left( \ltt \right)
   +
   \uz \ltz
   =
   -
   \frac{1}{\vr} \der{p}{\vt} \\
   &
   +
   \mu \left[
      \frac{1}{\vr} \der{}{\vr} \left(
         \vr \der{\ut}{\vr}
      \right)
      +
      \frac{1}{\vr} \der{}{\vt} \left(
         \frac{1}{\vr} \der{\ut}{\vt}
      \right)
      +
      \der{}{\vz} \left(
         \der{\ut}{\vz}
      \right)
      -
      \frac{\ut}{\vr^2}
      +
      \frac{2}{\vr^2} \der{\ur}{\vt}
   \right], \\
   \der{\uz}{t}
   &
   +
   \ur \lzr
   +
   \ut \lzt
   +
   \uz \lzz
   =
   -
   \der{p}{\vz} \\
   &
   +
   \mu \left[
      \frac{1}{\vr} \der{}{\vr} \left(
         \vr \der{\uz}{\vr}
      \right)
      +
      \frac{1}{\vr} \der{}{\vt} \left(
         \frac{1}{\vr} \der{\uz}{\vt}
      \right)
      +
      \der{}{\vz} \left(
         \der{\uz}{\vz}
      \right)
   \right].

Equation of the squared radial velocity

.. math::

   \der{\kr}{t}
   &
   +
   \ur \der{\kr}{\vr}
   +
   \ut \frac{1}{\vr} \der{\kr}{\vt}
   -
   \frac{\ur \ut^2}{\vr}
   +
   \uz \der{\kr}{\vz}
   =
   -
   \ur \der{p}{\vr} \\
   &
   +
   \mu \left[
      \frac{1}{\vr} \der{}{\vr} \left(
         \vr \der{\kr}{\vr}
      \right)
      +
      \frac{1}{\vr} \der{}{\vt} \left(
         \frac{1}{\vr} \der{\kr}{\vt}
      \right)
      +
      \der{}{\vz} \left(
         \der{\kr}{\vz}
      \right)
   \right] \\
   &
   -
   \mu \left[
      \left( \der{\ur}{\vr} \right)^2
      +
      \left( \frac{1}{\vr} \der{\ur}{\vt} \right)^2
      +
      \left( \der{\ur}{\vz} \right)^2
      +
      \left( \frac{\ur}{\vr} \right)^2
      +
      2 \frac{\ur}{\vr} \frac{1}{\vr} \der{\ut}{\vt}
   \right], \\
   \der{\kt}{t}
   &
   +
   \ur \der{\kt}{\vr}
   +
   \ut \frac{1}{\vr} \der{\kt}{\vt}
   +
   \frac{\ur \ut^2}{\vr}
   +
   \uz \der{\kt}{\vz}
   =
   -
   \ut \frac{1}{\vr} \der{p}{\vt} \\
   &
   +
   \mu \left[
      \frac{1}{\vr} \der{}{\vr} \left(
         \vr \der{\kt}{\vr}
      \right)
      +
      \frac{1}{\vr} \der{}{\vt} \left(
         \frac{1}{\vr} \der{\kt}{\vt}
      \right)
      +
      \der{}{\vz} \left(
         \der{\kt}{\vz}
      \right)
   \right] \\
   &
   -
   \mu \left[
      \left( \der{\ut}{\vr} \right)^2
      +
      \left( \frac{1}{\vr} \der{\ut}{\vt} \right)^2
      +
      \left( \der{\ut}{\vz} \right)^2
      +
      \left( \frac{\ut}{\vr} \right)^2
      -
      \left( 2 \frac{\ut}{\vr} \frac{1}{\vr} \der{\ur}{\vt} \right)
   \right], \\
   \der{\kz}{t}
   &
   +
   \ur \der{\kz}{\vr}
   +
   \ut \frac{1}{\vr} \der{\kz}{\vt}
   +
   \uz \der{\kz}{\vz}
   =
   -
   \uz \der{p}{\vz} \\
   &
   +
   \mu \left[
      \frac{1}{\vr} \der{}{\vr} \left(
         \vr \der{\kz}{\vr}
      \right)
      +
      \frac{1}{\vr} \der{}{\vt} \left(
         \frac{1}{\vr} \der{\kz}{\vt}
      \right)
      +
      \der{}{\vz} \left(
         \der{\kz}{\vz}
      \right)
   \right] \\
   &
   -
   \mu \left[
      \left( \der{\uz}{\vr} \right)^2
      +
      \left( \frac{1}{\vr} \der{\uz}{\vt} \right)^2
      +
      \left( \der{\uz}{\vz} \right)^2
   \right].

They are derived in the following pages.
Also the global energy balance and the resulting Nusselt number relations are discussed.

.. toctree::
   :maxdepth: 1

   tensor
   adv
   dif
   global
   nusselt/main

