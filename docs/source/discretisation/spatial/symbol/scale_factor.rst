#####################################
Scale factor and Jacobian determinant
#####################################

Discrete scale factors are given by

.. math::

    \sfact{1}
    &
    \equiv
    \Delta \vx,

    \sfact{2}
    &
    \equiv
    \Delta \vy,

    \sfact{3}
    &
    \equiv
    \Delta \vz,

for Cartesian coordinates (``is_curved == false``), while

.. math::

    \sfact{1}
    &
    \equiv
    \Delta \vr,

    \sfact{2}
    &
    \equiv
    \vr
    \Delta \vt,

    \sfact{3}
    &
    \equiv
    \Delta \vz,

for cylindrical coordinates (``is_curved == true``), respectively.

They are stored at the wall-normal cell faces and centers (except :math:`\sfact{3}`, which is constant across the whole domain):

.. myliteralinclude:: /../../src/domain.c
    :language: c
    :tag: x scale factors at x cell faces

.. myliteralinclude:: /../../src/domain.c
    :language: c
    :tag: x scale factors at x cell centers

.. myliteralinclude:: /../../src/domain.c
    :language: c
    :tag: y scale factors at x cell faces

.. myliteralinclude:: /../../src/domain.c
    :language: c
    :tag: y scale factors at x cell centers

Discrete Jacobian determinants, which are also defined at cell faces and centers, are simply the product of the local scale factors:

.. math::

    J
    \equiv
    \Pi_i
    \sfact{i}.

.. myliteralinclude:: /../../src/domain.c
    :language: c
    :tag: jacobian determinants at x cell faces

.. myliteralinclude:: /../../src/domain.c
    :language: c
    :tag: jacobian determinants at x cell centers

