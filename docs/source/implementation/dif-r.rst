#########################
Radial momentum diffusion
#########################

.. include:: /references.txt

I consider variables at :math:`\left( \ric, \rjc, \rkc \right)` if not specified.

The radial momentum diffusion is given as

.. math::

   2 \der{}{\vr} \lrr
   & +
   \frac{1}{\vr}
   \der{}{\vt} \left\{
      % \lrt + \ltr
      \frac{1}{\vr} \der{\ur}{\vt}
      +
      \vr \der{}{\vr} \left( \frac{\ut}{\vr} \right)
   \right\}
   +
   \der{}{\vz} \left(
      \lrz
      +
      \lzr
   \right) \\
   & +
   \frac{2}{\vr}
   \left\{
      \vr \der{}{\vr} \left(
         \frac{\ur}{\vr}
      \right)
      -
      \frac{1}{\vr}
      \der{\ut}{\vt}
   \right\}.

**************
Simplification
**************

I consider the sum of

.. math::

   \dder{}{\vr} \left(
      \dder{\ur}{\vr}
   \right),

.. math::

   \frac{1}{\vr}
   \dder{}{\vt}
   \left\{
      \vr \dder{}{\vr} \left( \frac{\ut}{\vr} \right)
   \right\}
   =
   \dder{}{\vr}
   \left(
      \frac{1}{\vr}
      \dder{\ut}{\vt}
   \right),

.. math::

   \dder{}{\vz}
   \left(
      \dder{\uz}{\vr}
   \right)
   =
   \dder{}{\vr}
   \left(
      \dder{\uz}{\vz}
   \right),

.. math::

   \frac{1}{\vr}
   \left\{
      \vr \dder{}{\vr} \left(
         \frac{\ur}{\vr}
      \right)
   \right\}
   =
   \dder{}{\vr} \left(
      \frac{\ur}{\vr}
   \right),

giving

.. math::

   \dder{}{\vr} \left(
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
   }{\rim}
   +
   \frac{
      \vat{\ur}{\rim}
   }{
      \vat{\vr}{\rim}
   }
   =
   \vat{
      \left(
         \frac{1}{\vr}
         \dder{\vr \ur}{\vr}
      \right)
   }{\rim},

   \vat{
      \dder{\ur}{\vr}
   }{\rip}
   +
   \frac{
      \vat{\ur}{\rip}
   }{
      \vat{\vr}{\rip}
   }
   =
   \vat{
      \left(
         \frac{1}{\vr}
         \dder{\vr \ur}{\vr}
      \right)
   }{\rip},

hold if

.. math::

   \vat{\ur}{\rim}
   \approx
   \frac{1}{2} \left(
      \vat{\ur}{\rimm}
      +
      \vat{\ur}{\ric}
   \right),

   \vat{\ur}{\rip}
   \approx
   \frac{1}{2} \left(
      \vat{\ur}{\ric}
      +
      \vat{\ur}{\ripp}
   \right),

and

.. math::

   \vat{r}{\rim}
   =
   \frac{1}{2} \left(
      \vat{r}{\rimm}
      +
      \vat{r}{\ric}
   \right),

   \vat{r}{\rip}
   =
   \frac{1}{2} \left(
      \vat{r}{\ric}
      +
      \vat{r}{\ripp}
   \right),

are satisfied.

Thus I obtain

.. math::

   \newcommand{\signfirst}{+}
   \newcommand{\first}{
      \dder{}{\vr}
      \left(
         \frac{1}{\vr}
         \dder{}{\vr}
         \left(
            \vr \ur
         \right)
      \right)
   }
   \newcommand{\signsecond}{+}
   \newcommand{\second}{
      \frac{1}{\vr}
      \dder{}{\vt}
      \left(
         \frac{1}{\vr}
         \dder{\ur}{\vt}
      \right)
   }
   \newcommand{\signthird}{+}
   \newcommand{\third}{
      \dder{}{\vz}
      \left(
         \dder{\ur}{\vz}
      \right)
   }
   \newcommand{\signfourth}{-}
   \newcommand{\fourth}{
      \frac{1}{\vr}
      \dder{}{\vt}
      \left(
         \frac{1}{\vr}
         \dder{}{\vr}
         \left(
            \vr \ut
         \right)
      \right)
   }
   \newcommand{\signfifth}{+}
   \newcommand{\fifth}{
      \dder{}{\vt}
      \left(
         \dder{}{\vr}
         \left(
            \frac{\ut}{\vr}
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

   - \frac{2}{\vr^2} \der{\ut}{\vt}.

**************
Implementation
**************

.. math::

   \signfirst
   \frac{1}{Re}
   \first

.. myliteralinclude:: /../../src/domain/coordinate.c
   :language: c
   :tag: lapuxx

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: diffused in r

.. math::

   \signsecond
   \frac{1}{Re}
   \second

.. myliteralinclude:: /../../src/domain/coordinate.c
   :language: c
   :tag: lapuxy

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: diffused in t

.. math::

   \signthird
   \frac{1}{Re}
   \third

.. myliteralinclude:: /../../src/domain/coordinate.c
   :language: c
   :tag: lapuxz

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: diffused in z

.. math::

   \signfourth
   \frac{1}{Re}
   \fourth

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: additional diffusive term 0

.. math::

   \signfifth
   \frac{1}{Re}
   \fifth

.. myliteralinclude:: /../../src/fluid/compute_rhs/ux.c
   :language: c
   :tag: additional diffusive term 1

