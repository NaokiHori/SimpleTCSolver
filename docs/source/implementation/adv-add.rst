#########################################################
Additional terms - centrifugal and Coriolis contributions
#########################################################

.. include:: /references.txt

When written explicitly, there are 11 terms in total which describe the advection of momentum.

I can split them into two categorises, 1. conservative and 2. non-conservative terms.

#. Conservative terms

   These nine terms are written as

   .. math::

      \der{}{\vr},

   .. math::

      \frac{1}{\vr} \der{}{\vt},

   or

   .. math::

      \der{}{\vz},

   which are in conservative forms and the discretisations are (relatively) straightforward.

#. Non-conservative terms

   There are two additional terms which are not in conservative form in the momentum equations, which are the centrifugal force in the radial direction and the Coriolis force in the azimuthal direction.
   They appear because of the *distortion* of the coordinate system and behave like body forces.

   It is worthwhile to pay attention to their effects on the global energy balance.
   By computing the inner product of the momentum equation and the velocity vector, one can easily see that these terms do not contribute to the changes in the kinetic energy because

   .. math::

      \ur \times \frac{\ut \ut}{\vr}

   and

   .. math::

      - \ut \times \frac{\ut \ur}{\vr}

   cancel out to each other *locally* and as a result these two terms dissapear from the equation of the energy balance.

   Unfortunately it is non-trivial to satisfy this local conservation property after being discretised, since I adopt a staggered grid arrangement where each velocity components are defined at different positions.
   However, it is still crucial to mimic this property (or I would easily introduce artificial energy source or sink in my code), and many literatures stand on a similar perspective (e.g. |FUKAGATA2002|, |MORINISHI2004|, |OUD2016|) such that the local centrifugal and Coriolis contributions cancel out to each other at a cell vertex (see |FUKAGATA2002|).

   One minor problem is they treat these terms differently, and it alters the other properties of the equations.
   I adopt the scheme proposed by |MORINISHI2004| for the time being, which does not conserve the angular velocity flux |ECKHARDT2007| since the scheme violates one identity which is necessary to derive the conservation law.
   Actually it is possible to describe these fictitious terms such that both the energy and the angular velocity flux properties are discretely satisfied, so I might change the scheme in the future.

