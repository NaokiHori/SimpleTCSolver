############################
Azimuthal momentum diffusion
############################

.. include:: /references.txt

I consider variables at :math:`\left( \tic, \tjc, \tkc \right)` if not specified.

The azimuthal momentum diffusion is given as

.. math::

   \der{}{\vr} \left(
      \vr \der{}{\vr} \left( \frac{\ut}{\vr} \right)
      +
      \frac{1}{\vr} \der{\ur}{\vt}
   \right)
   & +
   \frac{2}{\vr}
   \der{}{\vt} \left(
      \ltt
   \right)
   +
   \der{}{\vz} \left(
      \ltz
      +
      \lzt
   \right) \\
   & +
   \frac{2}{\vr}
   \left\{
      \vr \der{}{\vr} \left(
         \frac{\ut}{\vr}
      \right)
      +
      \frac{1}{\vr}
      \der{\ur}{\vt}
   \right\}.

**************
Simplification
**************

I can merge the first and the last term to obtain

.. math::

   \frac{1}{\vr^2} \dder{}{\vr} \left(
      \vr^3 \dder{}{\vr} \left( \frac{\ut}{\vr} \right)
      +
      \vr \dder{\ur}{\vt}
   \right)
   & +
   \frac{1}{\vr}
   \dder{}{\vt} \left(
      \frac{1}{\vr} \dder{\ut}{\vt}
   \right)
   +
   \frac{1}{\vr}
   \dder{}{\vt} \left(
      \frac{1}{\vr} \dder{\ut}{\vt}
   \right) \\
   & +
   \frac{2}{\vr}
   \dder{}{\vt} \left(
      \frac{\ur}{\vr}
   \right)
   +
   \dder{}{\vz} \left(
      \dder{\ut}{\vz}
      +
      \frac{1}{\vr}
      \dder{\uz}{\vt}
   \right).

Since discrete azimuthal derivative is commutable with the other operators, I obtain

.. math::

   \frac{1}{\vr^2} \dder{}{\vr} \left(
      \vr^3 \dder{}{\vr} \left( \frac{\ut}{\vr} \right)
   \right)
   & +
   \frac{1}{\vr} \dder{}{\vt} \frac{1}{\vr} \dder{}{\vr} \left(
      \vr \ur
   \right)
   +
   \frac{1}{\vr}
   \dder{}{\vt} \left(
      \frac{1}{\vr} \dder{\ut}{\vt}
   \right)
   +
   \frac{1}{\vr}
   \dder{}{\vt} \left(
      \frac{1}{\vr} \dder{\ut}{\vt}
   \right) \\
   & +
   \frac{2}{\vr}
   \dder{}{\vt} \left(
      \frac{\ur}{\vr}
   \right)
   +
   \dder{}{\vz} \left(
      \dder{\ut}{\vz}
   \right)
   +
   \frac{1}{\vr}
   \dder{}{\vt} \left(
      \dder{\uz}{\vz}
   \right).

Sum of the 2nd, 3rd, 7th terms leads

.. math::

   \frac{1}{\vr} \dder{}{\vt} \left(
      \frac{1}{\vr} \dder{\vr \ur}{\vr}
      +
      \frac{1}{\vr} \dder{\ut}{\vt}
      +
      \dder{\uz}{\vz}
   \right),

which vanishes because of the continuity.
Thus I am left with

.. math::

   \frac{1}{\vr^2} \dder{}{\vr} \left(
      \vr^3 \dder{}{\vr} \left( \frac{\ut}{\vr} \right)
   \right)
   +
   \frac{1}{\vr}
   \dder{}{\vt} \left(
      \frac{1}{\vr} \dder{\ut}{\vt}
   \right)
   +
   \frac{2}{\vr}
   \dder{}{\vt} \left(
      \frac{\ur}{\vr}
   \right)
   +
   \dder{}{\vz} \left(
      \dder{\ut}{\vz}
   \right),

or equivalently

.. math::

   \newcommand{\signfirst}{+}
   \newcommand{\first}{
      \frac{1}{\vr^2}
      \dder{}{\vr}
      \left(
         \vr^3
         \dder{}{\vr}
         \left(
            \frac{\ut}{\vr}
         \right)
      \right)
   }
   \newcommand{\signsecond}{+}
   \newcommand{\second}{
      \frac{1}{\vr}
      \dder{}{\vt}
      \left(
         \frac{1}{\vr}
         \dder{\ut}{\vt}
      \right)
   }
   \newcommand{\signthird}{+}
   \newcommand{\third}{
      \dder{}{\vz}
      \left(
         \dder{\ut}{\vz}
      \right)
   }
   \newcommand{\signfourth}{+}
   \newcommand{\fourth}{
      \frac{1}{\vr}
      \dder{}{\vt}
      \left(
         \frac{1}{\vr}
         \dder{}{\vr}
         \left(
            \vr \ur
         \right)
      \right)
   }
   \newcommand{\signfifth}{-}
   \newcommand{\fifth}{
      \dder{}{\vt}
      \left(
         \dder{}{\vr}
         \left(
            \frac{\ur}{\vr}
         \right)
      \right)
   }
   \signfirst
   \first
   \signsecond
   \second
   \signthird
   \third \\
   \signfourth
   \fourth
   \signfifth
   \fifth,

where the last two terms are obtained by following |BARBOSA2005|, which are mathematically (not numerically) equivalent to

.. math::

   + \frac{2}{\vr^2} \der{\ur}{\vt}.

**************
Implementation
**************

.. math::

   \signfirst
   \frac{1}{Re}
   \first

.. myliteralinclude:: /../../src/domain/init.c
   :language: c
   :tag: Laplacian w.r.t. uy in x

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: uy is diffused in x

.. math::

   \signsecond
   \frac{1}{Re}
   \second

.. myliteralinclude:: /../../src/domain/init.c
   :language: c
   :tag: Laplacian w.r.t. uy in y

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: uy is diffused in y

.. math::

   \signthird
   \frac{1}{Re}
   \third

.. myliteralinclude:: /../../src/domain/init.c
   :language: c
   :tag: Laplacian in z

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: uy is diffused in z

.. math::

   \signfourth
   \frac{1}{Re}
   \fourth

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: diffusive effect from Christoffel symbol, 1

.. math::

   \signfifth
   \frac{1}{Re}
   \fifth

.. myliteralinclude:: /../../src/fluid/compute_rhs/uy.c
   :language: c
   :tag: diffusive effect from Christoffel symbol, 2

