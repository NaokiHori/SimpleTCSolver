#########
Diffusion
#########

**********************
Wall-normal derivative
**********************

.. math::

    \dmomdif{1}{3}
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \vgt{1}{3}
    \right)
    +
    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \vgt{3}{1}
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
        \vgt{1}{3}
    \right)
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \dlxz
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: second-order derivative in x

.. myliteralinclude:: /../../src/fluid/predict/uz.c
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
        \vgt{3}{1}
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: diffusion in x, 1

**********************
Stream-wise derivative
**********************

.. math::

    \dmomdif{2}{3}
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{2}{3}
    \right)
    +
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vgt{3}{2}
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
        \vgt{2}{3}
    \right)
    =
    \mu
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \dlyz
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: second-order derivative in y

.. myliteralinclude:: /../../src/fluid/predict/uz.c
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
        \vgt{3}{2}
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: diffusion in y, 1

********************
Span-wise derivative
********************

.. math::

    \dmomdif{3}{3}
    =
    2
    \mu
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \vgt{3}{3}
    \right)
    =
    2
    \mu
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \dlzz
    \right)

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: second-order derivative in z

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: diffusion in z

