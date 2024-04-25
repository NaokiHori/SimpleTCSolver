#########
Advection
#########

*********************
Wall-normal component
*********************

.. math::

    \dmomadv{1}{1}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: ux is advected in x

*********************
Stream-wise component
*********************

.. math::

    \dmomadv{1}{2}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: ux is advected in y

*******************
Span-wise component
*******************

.. math::

    \dmomadv{1}{3}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: ux is advected in z

********************
Additional component
********************

.. math::

    \dmomadvx

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: centrifugal effect

