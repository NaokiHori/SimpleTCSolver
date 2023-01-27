# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

project = 'SimpleTCSolver'
copyright = '2023, Naoki Hori'
author = 'Naoki Hori'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.append(os.path.abspath("./ext"))
extensions = [
    "myliteralinclude",
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_theme_options = {
    "fixed_sidebar": "true",
    "github_user": "NaokiHori",
    "github_repo": "SimpleTCSolver",
    "github_type": "true",
}

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
mathjax3_config = {
    "TeX": {
        "Macros": {
            "vr": ["{r}"],
            "vt": ["{\\theta}"],
            "vz": ["{z}"],
            "ur": ["{u_{\\vr}}"],
            "ut": ["{u_{\\vt}}"],
            "uz": ["{u_{\\vz}}"],
            "dr": ["{\\Delta \\vr}"],
            "dt": ["{\\Delta \\vt}"],
            "dz": ["{\\Delta \\vz}"],
            "er": ["{\\underline{e}_{\\vr}}"],
            "et": ["{\\underline{e}_{\\vt}}"],
            "ez": ["{\\underline{e}_{\\vz}}"],
            "der": ["{\\frac{\\partial #1}{\\partial #2}}", 2],      # derivative
            "dder": ["{\\frac{\\delta #1}{\\delta #2}}", 2],         # discrete derivative
            "ave": ["{\\overline{#1}^{#2}}", 2], # interpolation
            "dif": ["{\\delta_{#2} {#1}}", 2],   # differentiation
            "vat": ["{\\left. {#1} \\right|_{#2}}", 2], # value at
            ## indices, pressure, x-face, y-face in two directions
            # p-i
            "pimm": ["i-1           "],
            "pim":  ["i-\\frac{1}{2}"],
            "pic":  ["i             "],
            "pip":  ["i+\\frac{1}{2}"],
            "pipp": ["i+1           "],
            # p-j
            "pjmm": ["j-1           "],
            "pjm":  ["j-\\frac{1}{2}"],
            "pjc":  ["j             "],
            "pjp":  ["j+\\frac{1}{2}"],
            "pjpp": ["j+1           "],
            # p-k
            "pkmm": ["k-1           "],
            "pkm":  ["k-\\frac{1}{2}"],
            "pkc":  ["k             "],
            "pkp":  ["k+\\frac{1}{2}"],
            "pkpp": ["k+1           "],
            # r-i
            "rimm": ["i-\\frac{1}{2}"],
            "rim":  ["i             "],
            "ric":  ["i+\\frac{1}{2}"],
            "rip":  ["i+1           "],
            "ripp": ["i+\\frac{3}{2}"],
            # r-j
            "rjmm": ["j-1           "],
            "rjm":  ["j-\\frac{1}{2}"],
            "rjc":  ["j             "],
            "rjp":  ["j+\\frac{1}{2}"],
            "rjpp": ["j+1           "],
            # r-k
            "rkmm": ["k-1           "],
            "rkm":  ["k-\\frac{1}{2}"],
            "rkc":  ["k             "],
            "rkp":  ["k+\\frac{1}{2}"],
            "rkpp": ["k+1           "],
            # t-i
            "timm": ["i-1           "],
            "tim":  ["i-\\frac{1}{2}"],
            "tic":  ["i             "],
            "tip":  ["i+\\frac{1}{2}"],
            "tipp": ["i+1           "],
            # t-j
            "tjmm": ["j-\\frac{1}{2}"],
            "tjm":  ["j             "],
            "tjc":  ["j+\\frac{1}{2}"],
            "tjp":  ["j+1           "],
            "tjpp": ["j+\\frac{3}{2}"],
            # t-k
            "tkmm": ["k-1           "],
            "tkm":  ["k-\\frac{1}{2}"],
            "tkc":  ["k             "],
            "tkp":  ["k+\\frac{1}{2}"],
            "tkpp": ["k+1           "],
            # z-i
            "zimm": ["i-1           "],
            "zim":  ["i-\\frac{1}{2}"],
            "zic":  ["i             "],
            "zip":  ["i+\\frac{1}{2}"],
            "zipp": ["i+1           "],
            # z-j
            "zjmm": ["j-1           "],
            "zjm":  ["j-\\frac{1}{2}"],
            "zjc":  ["j             "],
            "zjp":  ["j+\\frac{1}{2}"],
            "zjpp": ["j+1           "],
            # z-k
            "zkmm": ["k-\\frac{1}{2}"],
            "zkm":  ["k             "],
            "zkc":  ["k+\\frac{1}{2}"],
            "zkp":  ["k+1           "],
            "zkpp": ["k+\\frac{3}{2}"],
            # l_ij
            "lrr": ["\\der{\\ur}{\\vr}"],
            "lrt": ["\\frac{1}{\\vr} \\der{\\ur}{\\vt} - \\frac{\\ut}{\\vr}"],
            "lrz": ["\\der{\\ur}{\\vz}"],
            "ltr": ["\\der{\\ut}{\\vr}"],
            "ltt": ["\\frac{1}{\\vr} \\der{\\ut}{\\vt} + \\frac{\\ur}{\\vr}"],
            "ltz": ["\\der{\\ut}{\\vz}"],
            "lzr": ["\\der{\\uz}{\\vr}"],
            "lzt": ["\\frac{1}{\\vr} \\der{\\uz}{\\vt}"],
            "lzz": ["\\der{\\uz}{\\vz}"],
            # t_ij
            "trr": ["2 \\mu \\left( \\der{\\ur}{\\vr} \\right)"],
            "trt": ["\\mu \\left( \\frac{1}{\\vr} \\der{\\ur}{\\vt} + \\vr \\der{}{\\vr} \\left( \\frac{\\ut}{\\vr} \\right) \\right)"],
            "trz": ["\\mu \\left( \\der{\\ur}{\\vz} + \\der{\\uz}{\\vr} \\right)"],
            "ttr": ["\\trt"],
            "ttt": ["2 \\mu \\left( \\frac{1}{\\vr} \\der{\\ut}{\\vt} + \\frac{\\ur}{\\vr} \\right)"],
            "ttz": ["\\mu \\left( \\der{\\ut}{\\vz} + \\frac{1}{\\vr} \\der{\\uz}{\\vt} \\right)"],
            "tzr": ["\\trz"],
            "tzt": ["\\ttz"],
            "tzz": ["2 \\mu \\left( \\der{\\uz}{\\vz} \\right)"],
        }
    }
}

