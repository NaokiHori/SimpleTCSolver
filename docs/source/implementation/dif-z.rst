########################
Axial momentum diffusion
########################

I consider variables at :math:`\left( \zic, \zjc, \zkc \right)` if not specified.

The axial momentum diffusion is given as

.. math::

   \der{}{\vr} \left(
      \lzr
      +
      \lrz
   \right)
   +
   \frac{1}{\vr}
   \der{}{\vt} \left(
      \lzt
      +
      \ltz
   \right)
   +
   \der{}{\vz} \left(
      \lzz
      +
      \lzz
   \right)
   +
   \frac{1}{\vr} \left(
      \lzr
      +
      \lrz
   \right).

**************
Simplification
**************

I consider the sum of

.. math::

   \dder{}{\vr} \left(
      \dder{\ur}{\vz}
   \right)
   =
   \dder{}{\vz} \left(
      \dder{\ur}{\vr}
   \right),

.. math::

   \frac{1}{\vr}
   \dder{}{\vt} \left(
      \dder{\ut}{\vz}
   \right)
   =
   \dder{}{\vz}
   \left(
      \frac{1}{\vr}
      \dder{\ut}{\vt}
   \right),

.. math::

   \dder{}{\vz} \left(
      \dder{\uz}{\vz}
   \right),

.. math::

   \frac{1}{\vr}
   \dder{\ur}{\vz}
   =
   \dder{}{\vz} \left(
      \frac{\ur}{\vr}
   \right),

giving

.. math::

   \dder{}{\vz} \left(
      \dder{\ur}{\vr}
      +
      \frac{1}{\vr}
      \dder{\ut}{\vt}
      +
      \dder{\uz}{\vz}
      +
      \frac{\ur}{\vr}
   \right),

which vanishes thanks to the incompressibility constraint since

.. math::

   \vat{
      \dder{\ur}{\vr}
   }{\zic}
   +
   \frac{
      \vat{\ur}{\zic}
   }{
      \vat{\vr}{\zic}
   }
   =
   \vat{
      \left(
         \frac{1}{\vr}
         \dder{\vr \ur}{\vr}
      \right)
   }{\zic}

holds if

.. math::

   \vat{\ur}{\zic}
   \approx
   \frac{1}{2} \left(
      \vat{\ur}{\zim}
      +
      \vat{\ur}{\zip}
   \right)

and

.. math::

   \vat{r}{\zic}
   =
   \frac{1}{2} \left(
      \vat{r}{\zim}
      +
      \vat{r}{\zip}
   \right)

are satisfied.

These relations are used to deduce

.. math::

   \dder{}{\vr}
   \left(
      \dder{\uz}{\vr}
   \right)
   +
   \frac{1}{\vr}
   \left(
      \dder{\uz}{\vr}
   \right)
   =
   \frac{1}{\vr}
   \dder{}{\vr} \left(
      \vr \dder{\uz}{\vr}
   \right)

at :math:`\zic` (:math:`x` cell faces) and as a result I obtain

.. math::

   \newcommand{\first}{
      \frac{1}{\vr}
      \dder{}{\vr} \left(
         \vr \dder{\uz}{\vr}
      \right)
   }
   \newcommand{\second}{
      \frac{1}{\vr}
      \dder{}{\vt} \left(
         \frac{1}{\vr} \dder{\uz}{\vt}
      \right)
   }
   \newcommand{\third}{
      \dder{}{\vz} \left(
         \dder{\uz}{\vz}
      \right)
   }
   \first
   +
   \second
   +
   \third.

**************
Implementation
**************

.. math::

   +
   \frac{1}{Re}
   \first

.. myliteralinclude:: /../../src/domain/coordinate.c
   :language: c
   :tag: uzdifx

.. myliteralinclude:: /../../src/fluid/compute_rhs/uz.c
   :language: c
   :tag: diffused in r

.. math::

   +
   \frac{1}{Re}
   \second

.. myliteralinclude:: /../../src/fluid/compute_rhs/uz.c
   :language: c
   :tag: diffused in t

.. math::

   +
   \frac{1}{Re}
   \third

.. myliteralinclude:: /../../src/fluid/compute_rhs/uz.c
   :language: c
   :tag: diffused in z

