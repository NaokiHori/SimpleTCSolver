########################
Velocity-gradient tensor
########################

For later convenience, 9 tensor components are pre-computed and stored to arrays.

**
11
**

.. math::

    \vat{
        \vgt{1}{1}
    }{
        \ccindex{i},
        \ccindex{j},
        \ccindex{k}
    }
    =
    \dlxx

.. myliteralinclude:: /../../src/fluid/predict/lxx.c
   :language: c
   :tag: compute lxx

**
21
**

.. math::

    \vat{
        \vgt{2}{1}
    }{
        \cpindex{i},
        \cpindex{j},
        \ccindex{k}
    }
    =
    \vat{
        \vgt{2}{1}^{\prime}
    }{
        \cpindex{i},
        \cpindex{j},
        \ccindex{k}
    }
    +
    \vat{
        \vgt{2}{1}^{\prime\prime}
    }{
        \cpindex{i},
        \cpindex{j},
        \ccindex{k}
    }
    =
    \dlyx

Two terms are computed and stored separately.

.. math::

    \vat{
        \vgt{2}{1}^{\prime}
    }{
        \cpindex{i},
        \cpindex{j},
        \ccindex{k}
    }
    =
    \dlyxa

.. myliteralinclude:: /../../src/fluid/predict/lyx.c
    :language: c
    :tag: compute dominant term of lyx

.. math::

    \vat{
        \vgt{2}{1}^{\prime\prime}
    }{
        \cpindex{i},
        \cpindex{j},
        \ccindex{k}
    }
    =
    \dlyxb

.. myliteralinclude:: /../../src/fluid/predict/lyx.c
    :language: c
    :tag: compute non-dominant term of lyx

**
31
**

.. math::

    \vat{
        \vgt{3}{1}
    }{
        \cpindex{i},
        \ccindex{j},
        \cpindex{k}
    }
    =
    \dlzx

.. myliteralinclude:: /../../src/fluid/predict/lzx.c
    :language: c
    :tag: compute lzx

**
12
**

.. math::

    \vat{
        \vgt{1}{2}
    }{
        \cpindex{i},
        \cpindex{j},
        \ccindex{k}
    }
    =
    \dlxy

.. myliteralinclude:: /../../src/fluid/predict/lxy.c
    :language: c
    :tag: compute lxy

**
22
**

.. math::

    \vat{
        \vgt{2}{2}
    }{
        \ccindex{i},
        \ccindex{j},
        \ccindex{k}
    }
    =
    \vat{
        \vgt{2}{2}^{\prime}
    }{
        \ccindex{i},
        \ccindex{j},
        \ccindex{k}
    }
    +
    \vat{
        \vgt{2}{2}^{\prime\prime}
    }{
        \ccindex{i},
        \ccindex{j},
        \ccindex{k}
    }
    =
    \dlyy

Two terms are computed and stored separately.

.. math::

    \vat{
        \vgt{2}{2}^{\prime}
    }{
        \ccindex{i},
        \ccindex{j},
        \ccindex{k}
    }
    =
    \dlyya

.. myliteralinclude:: /../../src/fluid/predict/lyy.c
    :language: c
    :tag: compute dominant term of lyy

.. math::

    \vat{
        \vgt{2}{2}^{\prime\prime}
    }{
        \ccindex{i},
        \ccindex{j},
        \ccindex{k}
    }
    =
    \dlyyb

.. myliteralinclude:: /../../src/fluid/predict/lyy.c
    :language: c
    :tag: compute non-dominant term of lyy

**
32
**

.. math::

    \vat{
        \vgt{3}{2}
    }{
        \ccindex{i},
        \cpindex{j},
        \cpindex{k}
    }
    =
    \dlzy

.. myliteralinclude:: /../../src/fluid/predict/lzy.c
   :language: c
   :tag: compute lzy

**
13
**

.. math::

    \vat{
        \vgt{1}{3}
    }{
        \cpindex{i},
        \ccindex{j},
        \cpindex{k}
    }
    =
    \dlxz

.. myliteralinclude:: /../../src/fluid/predict/lxz.c
    :language: c
    :tag: compute lxz

**
23
**

.. math::

    \vat{
        \vgt{2}{3}
    }{
        \ccindex{i},
        \cpindex{j},
        \cpindex{k}
    }
    =
    \dlyz

.. myliteralinclude:: /../../src/fluid/predict/lyz.c
    :language: c
    :tag: compute lyz

**
33
**

.. math::

    \vat{
        \vgt{3}{3}
    }{
        \ccindex{i},
        \ccindex{j},
        \ccindex{k}
    }
    =
    \dlzz

.. myliteralinclude:: /../../src/fluid/predict/lzz.c
    :language: c
    :tag: compute lzz

