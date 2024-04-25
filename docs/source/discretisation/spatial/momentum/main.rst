
.. _discrete_momentum:

################
Momentum balance
################

The wall-normal, stream-wise, and span-wise momentum equations are defined at :math:`\gcs{1}`, :math:`\gcs{2}`, and :math:`\gcs{3}` cell faces, respectively.

Note that, to circumvent severe time-step constraints, diagonal terms of the diffusive terms can be treated implicitly in time, which makes the implementation a bit complex.

.. toctree::
    :maxdepth: 1

    x/main.rst
    y/main.rst
    z/main.rst

