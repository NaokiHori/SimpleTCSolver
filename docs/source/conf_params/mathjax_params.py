from mathjax_macro.general import add as add_general
from mathjax_macro.momentum import add as add_momentum
from mathjax_macro.scalar import add as add_scalar
from mathjax_macro.discrete import add as add_discrete
from mathjax_macro.derivation import add as add_derivation

macros = dict()

add_general(macros)
add_momentum(macros)
add_scalar(macros)
add_discrete(macros)
add_derivation(macros)

mathjax3_config = {"tex": {"macros": macros}}
