###############
Correction step
###############

Updated but still non-solenoidal velocity field :math:`u_i^*` is corrected to be divergence-free by solving a Poisson equation with respect to the scalar potential :math:`\phi`:

.. math::

    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \frac{1}{\sfact{1}}
        \dif{\phi}{\gcs{1}}
    \right)
    +
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \frac{1}{\sfact{2}}
        \dif{\phi}{\gcs{2}}
    \right)
    +
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \frac{1}{\sfact{3}}
        \dif{\phi}{\gcs{3}}
    \right)
    =
    \frac{1}{\Delta t}
    \dif{\vel{i}^*}{\gcs{i}}.

The eigenvalues are computed as follows:

.. myliteralinclude:: /../../src/fluid/compute_potential.c
    :language: c
    :tag: azimuthal eigenvalues

.. myliteralinclude:: /../../src/fluid/compute_potential.c
    :language: c
    :tag: axial eigenvalues

In the radial direction, tri-diagonal matrices are to be solved for each :math:`y` and :math:`z`:

.. myliteralinclude:: /../../src/fluid/compute_potential.c
    :language: c
    :tag: solve tri-diagonal matrix

After computing the scalar potential :math:`\phi`, the velocity field is corrected using its gradient:

.. math::

    \pder{\phi}{x_i}
    \approx
    \vec{e_{\gcs{1}}}
    \frac{1}{\sfact{1}}
    \dif{\phi}{\gcs{1}}
    +
    \vec{e_{\gcs{2}}}
    \frac{1}{\sfact{2}}
    \dif{\phi}{\gcs{2}}
    +
    \vec{e_{\gcs{3}}}
    \frac{1}{\sfact{3}}
    \dif{\phi}{\gcs{3}}.

.. myliteralinclude:: /../../src/fluid/correct/ux.c
    :language: c
    :tag: correct x velocity

.. myliteralinclude:: /../../src/fluid/correct/uy.c
    :language: c
    :tag: correct y velocity

.. myliteralinclude:: /../../src/fluid/correct/uz.c
    :language: c
    :tag: correct z velocity

