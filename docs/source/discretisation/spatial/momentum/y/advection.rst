#########
Advection
#########

*********************
Wall-normal component
*********************

.. math::

    \dmomadv{2}{1}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: uy is advected in x

*********************
Stream-wise component
*********************

.. math::

    \dmomadv{2}{2}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: uy is advected in y

*******************
Span-wise component
*******************

.. math::

    \dmomadv{2}{3}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: uy is advected in z

********************
Additional component
********************

.. math::

    \dmomadvy

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: coriolis effect

