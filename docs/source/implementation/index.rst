##############
Implementation
##############

.. include:: /references.txt

Here, how the governing equations are discretised is discussed.

.. note::

   Discretisations are extremely complicated compared to the Cartesian version because of the multiplications by :math:`r` or :math:`1 / r` and the additional terms corresponding to the Christoffel symbols.
   It is non-trivial to keep some basic mathematical properties of the equations (see e.g., |ECKHARDT2007|).
   In this library, I propose one solution which satisfies some properties while others are violated.

.. toctree::
   :caption: Implementation - advective terms
   :maxdepth: 1

   adv-r
   adv-t
   adv-z

.. toctree::
   :caption: Implementation - diffusive terms
   :maxdepth: 1

   dif-r
   dif-t
   dif-z

