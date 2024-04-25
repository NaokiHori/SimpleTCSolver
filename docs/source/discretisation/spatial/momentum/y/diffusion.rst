#########
Diffusion
#########

**********************
Wall-normal derivative
**********************

.. math::

    \dmomdif{1}{2}
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \vgt{1}{2}
    \right)
    +
    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \vgt{2}{1}
    \right)

=============
Diagonal part
=============

.. math::

    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \vgt{1}{2}
    \right)
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \dlxy
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: second-order derivative in x

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: diffusion in x, 0

=================
Non-diagonal part
=================

.. math::

    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \vgt{2}{1}
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: diffusion in x, 1

**********************
Stream-wise derivative
**********************

.. math::

    \dmomdif{2}{2}
    =
    2
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{2}{2}
    \right)

=============
Diagonal part
=============

.. math::

    2
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{2}{2}^{\prime}
    \right)
    =
    2
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \frac{1}{\sfact{2}}
        \dif{\vel{2}}{\gcs{2}}
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: second-order derivative in y

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: diffusion in y, 0

=================
Non-diagonal part
=================

.. math::

    2
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{2}{2}^{\prime\prime}
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: diffusion in y, 1

********************
Span-wise derivative
********************

.. math::

    \dmomdif{3}{2}
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \vgt{3}{2}
    \right)
    +
    \mu
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \vgt{2}{3}
    \right)

=============
Diagonal part
=============

.. math::

    \mu
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \dlzy
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: second-order derivative in z

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: diffusion in z, 0

=================
Non-diagonal part
=================

.. math::

    \mu
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \vgt{2}{3}
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: diffusion in z, 1

***************
Additional term
***************

.. math::

    \dmomdify
    =
    \mu
    \frac{1}{J}
    \ave{
        \dif{
            \left(
                \frac{J}{\sfact{1}}
            \right)
        }{\gcs{1}}
        \left(
            \vgt{2}{1}
            +
            \vgt{1}{2}
        \right)
    }{\gcs{1}}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: additional diffusion

