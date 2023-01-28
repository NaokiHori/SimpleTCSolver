
.. _implementation:

##############
Implementation
##############

.. include:: /references.txt

Here, how the governing equations are discretised is discussed.

.. note::

   Let me consider

   .. math::

      \der{}{\vr} \left( \frac{1}{\vr} \right),

   which is obviously equal to

   .. math::

      - \frac{1}{\vr^2}.

   This relation is, however, not true numerically, and instead I should have

   .. math::

      \vat{
         \left\{ \dder{}{\vr} \left( \frac{1}{\vr} \right) \right\}
      }{\ric}
      =
      - \frac{1}{\vat{\vr}{\rim} \vat{\vr}{\rip}}
      \ne
      - \frac{1}{\left( \vat{\vr}{\ric} \right)^2}

   when I adopt the classical second-order-accurate central-difference scheme.

   Hereafter :math:`\frac{\delta}{\delta \vr}` is used to imply a discrete derivative (i.e. differentiation), and the same for the other directions :math:`\vt`, :math:`\vz`.

   As discussed in `the Cartesian counterpart <https://naokihori.github.io/SimpleNavierStokesSolver/>`_, consistency (e.g. energy conservation property) in discrete space is fragile, and as a result I definitely violate the property which the orginal governing equations have if I accept this kind of commutation error.

   As seen in the above example, discretising the equations *perfectly* in cylindrical coordinates is extremely challenging because of these commutation errors and the additional terms corresponding to the Christoffel symbols, and as a consequence it is non-trivial to keep mathematical properties of the equations (see e.g., |ECKHARDT2007|).
   In this library, I propose one solution which satisfies some properties which are important *in my opinion*, while others are violated.

.. seealso::

   There is an idea called `mimetic finite difference`, which focuses on discretising the equations while keeping the properties they have in the continuous forms.
   Although this idea is indeed adopted to discretise the Navier-Stokes equations in cylindrical coordinates (e.g. |BARBOSA2005|, |OUD2016|), it generally depends on the identity

   .. math::

      \partial_j \partial_j u_i
      \equiv
      \partial_i \partial_j u_j
      -
      \epsilon_{ijk} \partial_j \epsilon_{klm} \partial_l u_m,

   which is only applicable when the dynamic viscosity :math:`\mu` is constant.

.. toctree::
   :caption: Implementation - advective terms
   :maxdepth: 1

   adv-add
   adv-r
   adv-t
   adv-z

.. toctree::
   :caption: Implementation - diffusive terms
   :maxdepth: 1

   dif-r
   dif-t
   dif-z

