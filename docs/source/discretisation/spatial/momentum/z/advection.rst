#########
Advection
#########

*********************
Wall-normal component
*********************

.. math::

    \dmomadv{3}{1}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: uz is advected in x

*********************
Stream-wise component
*********************

.. math::

    \dmomadv{3}{2}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: uz is advected in y

*******************
Span-wise component
*******************

.. math::

    \dmomadv{3}{3}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: uz is advected in z

