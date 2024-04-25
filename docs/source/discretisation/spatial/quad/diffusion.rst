#########
Diffusion
#########

********
Momentum
********

By utilising :ref:`the relations listed in the prerequisite <energy_prerequisite>`, the discrete diffusive terms in the momentum balances yield the following relations.

From the wall-normal momentum balance, we obtain

.. math::

    &
    -
    \sumzc
    \sumyc
    \sumxc
    J
    \vgt{1}{1}
    \sst{1}{1},

    &
    -
    \sumzc
    \sumyf
    \sumxf
    J
    \left(
        \frac{1}{\sfact{2}}
        \dif{\vel{1}}{\gcs{2}}
    \right)
    \sst{2}{1},

    &
    -
    \sumzf
    \sumyc
    \sumxf
    J
    \vgt{3}{1}
    \sst{3}{1},

    &
    -
    \sumzc
    \sumyc
    \sumxc
    J
    \left\{
        \frac{1}{J}
        \dif{
            \left(
                \frac{J}{\sfact{1}}
            \right)
        }{\gcs{1}}
        \ave{
            \vel{1}
        }{\gcs{1}}
    \right\}
    \sst{2}{2}.

From the stream-wise momentum balance, we obtain

.. math::

    &
    -
    \sumzc
    \sumyf
    \vat{
       \left(
          \frac{J}{\sfact{1}}
          \vel{2}
          \sst{1}{2}
       \right)
    }{i = \frac{1}{2}}
    -
    \sumzc
    \sumyf
    \sumxf
    J
    \left(
       \frac{1}{\sfact{1}}
       \dif{\vel{2}}{\gcs{1}}
    \right)
    \sst{1}{2}
    +
    \sumzc
    \sumyf
    \vat{
       \left(
          \frac{J}{\sfact{1}}
          \vel{2}
          \sst{1}{2}
       \right)
    }{i = \ngp{1} + \frac{1}{2}},

    &
    -
    \sumzc
    \sumyc
    \sumxc
    J
    \left(
       \frac{1}{\sfact{2}}
       \dif{\vel{2}}{\gcs{2}}
    \right)
    \sst{2}{2},

    &
    -
    \sumzf
    \sumyf
    \sumxc
    J
    \vgt{3}{2}
    \sst{3}{2},

    &
    \sumzc
    \sumyf
    \vat{J}{\frac{1}{2}}
    \left\{
       \frac{1}{\vat{J}{\frac{1}{2}}}
       \vat{
          \dif{
             \left(
                \frac{J}{\sfact{1}}
             \right)
          }{
             \gcs{1}
          }
       }{\frac{1}{2}}
       \frac{\vat{\vel{2}}{1}}{2}
    \right\}
    \vat{
       \sst{2}{1}
    }{\frac{1}{2}}
    +
    \sumzc
    \sumyf
    \sum_{i = \frac{3}{2}}^{\ngp{1} - \frac{1}{2}}
    J
    \left\{
       \frac{1}{J}
       \dif{
          \left(
             \frac{J}{\sfact{1}}
          \right)
       }{
          \gcs{1}
       }
       \ave{\vel{2}}{\gcs{1}}
    \right\}
    \sst{2}{1}
    +
    \sumzc
    \sumyf
    \vat{J}{\ngp{1} + \frac{1}{2}}
    \left\{
       \frac{1}{\vat{J}{\ngp{1} + \frac{1}{2}}}
       \vat{
          \dif{
             \left(
                \frac{J}{\sfact{1}}
             \right)
          }{
             \gcs{1}
          }
       }{\ngp{1} + \frac{1}{2}}
       \frac{\vat{\vel{2}}{\ngp{1}}}{2}
    \right\}
    \vat{
       \sst{2}{1}
    }{\ngp{1} + \frac{1}{2}}.

From the span-wise momentum balance, we obtain

.. math::

    &
    -
    \sumzf
    \sumyc
    \sumxf
    J
    \vgt{1}{3}
    \sst{1}{3},

    &
    -
    \sumzf
    \sumyf
    \sumxc
    J
    \vgt{2}{3}
    \sst{2}{3},

    &
    -
    \sumzc
    \sumyc
    \sumxc
    J
    \vgt{3}{3}
    \sst{3}{3}.

In total, we have

.. math::

    &
    -
    \sumzc
    \sumyc
    \sumxc
    J
    \vgt{1}{1}
    \sst{1}{1}
    -
    \sumzc
    \sumyf
    \sumxf
    J
    \vgt{2}{1}
    \sst{2}{1}
    -
    \sumzf
    \sumyc
    \sumxf
    J
    \vgt{3}{1}
    \sst{3}{1},

    &
    -
    \sumzc
    \sumyf
    \sumxf
    J
    \vgt{1}{2}
    \sst{1}{2}
    -
    \sumzc
    \sumyc
    \sumxc
    J
    \vgt{2}{2}
    \sst{2}{2}
    -
    \sumzf
    \sumyf
    \sumxc
    J
    \vgt{3}{2}
    \sst{3}{2},

    &
    -
    \sumzf
    \sumyc
    \sumxf
    J
    \vgt{1}{3}
    \sst{1}{3}
    -
    \sumzf
    \sumyf
    \sumxc
    J
    \vgt{2}{3}
    \sst{2}{3}
    -
    \sumzc
    \sumyc
    \sumxc
    J
    \vgt{3}{3}
    \sst{3}{3}

as the dissipative terms, while

.. math::

    -
    \sumzc
    \sumyf
    \vat{
       \left(
          \frac{J}{\sfact{1}}
          \vel{2}
          \sst{1}{2}
       \right)
    }{i = \frac{1}{2}}
    +
    \sumzc
    \sumyf
    \vat{
       \left(
          \frac{J}{\sfact{1}}
          \vel{2}
          \sst{1}{2}
       \right)
    }{i = \ngp{1} + \frac{1}{2}}

as the energy transfer on the walls, respectively.

******
Scalar
******

The global diffusive effects of the scalar transport on the quadratic quantity lead to

.. math::

    \sumzc
    \sumyc
    \sumxc
    J
    T
    \left\{
        \dscalardif{1}
        \dscalardif{2}
        \dscalardif{3}
    \right\},

giving

.. math::

    -
    \sumzc
    \sumyc
    \sumxf
    J
    \kappa
    \frac{1}{\sfact{1}}
    \dif{T}{\gcs{1}}
    \frac{1}{\sfact{1}}
    \dif{T}{\gcs{1}}
    -
    \sumzc
    \sumyf
    \sumxc
    J
    \kappa
    \frac{1}{\sfact{2}}
    \dif{T}{\gcs{2}}
    \frac{1}{\sfact{2}}
    \dif{T}{\gcs{2}}
    -
    \sumzf
    \sumyc
    \sumxc
    J
    \kappa
    \frac{1}{\sfact{3}}
    \dif{T}{\gcs{3}}
    \frac{1}{\sfact{3}}
    \dif{T}{\gcs{3}}

as the dissipative terms, while

.. math::

    -
    \sumzc
    \sumyc
    \vat{
        \left(
            \frac{J}{\sfact{1}}
            \kappa
            T
            \frac{1}{\sfact{1}}
            \dif{T}{\gcs{1}}
        \right)
    }{\frac{1}{2}}
    +
    \sumzc
    \sumyc
    \vat{
        \left(
            \frac{J}{\sfact{1}}
            \kappa
            T
            \frac{1}{\sfact{1}}
            \dif{T}{\gcs{1}}
        \right)
    }{\ngp{1} + \frac{1}{2}}

as the transport on the walls.

