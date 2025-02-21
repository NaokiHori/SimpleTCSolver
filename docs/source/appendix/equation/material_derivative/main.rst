###################
Material derivative
###################

By using

.. math::

   d\vec{u}
   =
   \sum_j
   \pder{\vec{u}}{Q^j}
   dQ^j,

the material derivative of a vector :math:`\vec{u}` with respect to :math:`t` (one of the variables in :math:`Q^j`) leads to

.. math::

   \mder{\vec{u}}
   &
   =
   \pder{\vec{u}}{t}
   +
   \sum_j
   \tder{X^j}{t}
   \pder{\vec{u}}{X^j} \\
   &
   =
   \pder{\vec{u}}{t}
   +
   \sum_j
   V^j
   \pder{\vec{u}}{X^j},

where

.. math::

   V^j
   \equiv
   \tder{X^j}{t}

is introduced to distinguish the advected (:math:`U^i`) and the advecting (:math:`V^i`) velocity components.

.. toctree::
   :maxdepth: 1

   covariant
   normal

