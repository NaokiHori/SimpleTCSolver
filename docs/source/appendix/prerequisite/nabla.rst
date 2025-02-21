
.. _nabla:

##############
Nabla operator
##############

***********
Description
***********

I consider the following operator in the Cartesian coordinate:

.. math::

   \vec{\nabla}
   \equiv
   \sum_i
   \vec{e}_i
   \pder{}{x^i},

which frequently appears and can operate upon arbitrary order of tensors.

Now I aim at describing this relation on the general orthogonal coordinate systems.
Using :ref:`the basis vector transform <from_g_to_c>`

.. math::

   \fromgtoc

and the chain rule, I notice

.. math::

   \sum_i
   \vec{e}_i
   \pder{}{x^i}
   =
   \sum_{ijk}
   \vec{E}_j
   \pder{X^j}{x^i}
   \pder{X^k}{x^i}
   \pder{}{X^k}.

By using :ref:`the relation of the transformation matrices <jacobi_conv>`:

.. math::

   \jacobiconv,

I have

.. math::

   \sum_i
   \pder{X^j}{x^i}
   \pder{X^k}{x^i}
   =
   \frac{1}{H_j H_j}
   \frac{1}{H_k H_k}
   \sum_i
   \pder{x^i}{X^j}
   \pder{x^i}{X^k},

and by adopting :ref:`the relation of the metric tensor <metric_tensor>`:

.. math::

   \metrictensor,

this yields

.. math::

   \frac{1}{H_j H_j}
   \frac{1}{H_k H_k}
   H_j
   H_k
   \delta_{jk}
   =
   \frac{1}{H_j}
   \frac{1}{H_k}
   \delta_{jk}.

Thus

.. math::

   \sum_{jk}
   \vec{E}_j
   \frac{1}{H_j}
   \frac{1}{H_k}
   \delta_{jk}
   \pder{}{X^k}
   &
   =
   \sum_j
   \vec{E}_j
   \frac{1}{H_j H_j}
   \pder{}{X^j} \\
   &
   =
   \sum_j
   \vec{E}^j
   \pder{}{X^j}.

In summary,

.. math::

   \vec{\nabla}
   &
   \equiv
   \sum_i
   \vec{e}_i
   \pder{}{x^i} \\
   &
   =
   \sum_i
   \vec{E}^i
   \pder{}{X^i} \\
   &
   =
   \sum_i
   \vec{\hat{E}}_i
   \frac{1}{H_i}
   \pder{}{X^i}.

