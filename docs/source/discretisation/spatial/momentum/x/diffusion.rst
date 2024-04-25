#########
Diffusion
#########

**********************
Wall-normal derivative
**********************

.. math::

    \dmomdif{1}{1}
    =
    2
    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \vgt{1}{1}
    \right)
    =
    2
    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \dlxx
    \right)

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: second-order derivative in x

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: diffusion in x

**********************
Stream-wise derivative
**********************

.. math::

    \dmomdif{2}{1}
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{2}{1}
    \right)
    +
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{1}{2}
    \right)

=============
Diagonal part
=============

.. math::

    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{2}{1}^{\prime}
    \right)
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \dlyxa
    \right)

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: second-order derivative in y

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: diffusion in y, 0

=================
Non-diagonal part
=================

.. math::

    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{2}{1}^{\prime\prime}
    \right)
    +
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{1}{2}
    \right)

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: diffusion in y, 1

********************
Span-wise derivative
********************

.. math::

    \dmomdif{3}{1}
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \vgt{3}{1}
    \right)
    +
    \mu
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \vgt{1}{3}
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
        \vgt{3}{1}
    \right)
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \dlzx
    \right)

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: second-order derivative in z

.. myliteralinclude:: /../../src/fluid/predict/ux.c
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
        \vgt{1}{3}
    \right)

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: diffusion in z, 1

***************
Additional term
***************

.. math::

    \dmomdifx
    =
    -
    2
    \mu
    \frac{1}{J}
    \ave{
        \dif{
            \left(
                \frac{J}{\sfact{1}}
            \right)
        }{\gcs{1}}
        \vgt{2}{2}
    }{\gcs{1}}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: additional diffusion

