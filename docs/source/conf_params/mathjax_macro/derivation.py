def add(macros: dict):
    macros["vectorcompbasis"] = (
        "\\vec{r}"
        "="
        "\\sum_i"
        "x^i"
        "\\vec{e}_i"
        "="
        "\\sum_i"
        "X^i"
        "\\vec{E}_i"
    )
    macros["fromctog"] = (
        "\\vec{E}_i"
        "="
        "\\sum_j"
        "\\pder{x^j}{X^i}"
        "\\vec{e}_j"
    )
    macros["fromgtoc"] = (
        "\\vec{e}_i"
        "="
        "\\sum_j"
        "\\pder{X^j}{x^i}"
        "\\vec{E}_j"
    )
    macros["jacobiconv"] = (
        "\\pder{x^i}{X^j}"
        "="
        "\\pder{X^j}{x^i}"
        "H_j H_j"
    )
    macros["orthogonal"] = (
        "\\vec{E}_i"
        "\\cdot"
        "\\vec{E}_j"
        "="
        "H_i"
        "H_j"
        "\\delta_{ij}"
    )
    macros["metrictensor"] = (
        "H_i"
        "H_j"
        "\\delta_{ij}"
        "="
        "\\sum_k"
        "\\pder{x^k}{X^i}"
        "\\pder{x^k}{X^j}"
    )
    macros["dedx"] = [
        (
            "\\pder{\\vec{E}_#1}{X^#2}"
            "="
            "-"
            "\\sum_#3"
            "\\delta_{#1 #2}"
            "\\frac{H_#2}{2 H_#3 H_#3}"
            "\\pder{H_#1}{X^#3}"
            "\\vec{E}_#3"
            "-"
            "\\sum_#3"
            "\\delta_{#1 #2}"
            "\\frac{H_#1}{2 H_#3 H_#3}"
            "\\pder{H_#2}{X^#3}"
            "\\vec{E}_#3"
            "+"
            "\\frac{1}{H_#1}"
            "\\pder{H_#1}{X^#2}"
            "\\vec{E}_#1"
            "+"
            "\\frac{1}{H_#2}"
            "\\pder{H_#2}{X^#1}"
            "\\vec{E}_#2"
        ),
        3
    ]
    macros["dehatdx"] = [
        (
            "\\pder{\\hat{E}_#1}{X^#2}"
            "="
            "-"
            "\\sum_#3"
            "\\delta_{#1 #2}"
            "\\frac{1}{H_#3}"
            "\\pder{H_#1}{X^#3}"
            "\\vec{\\hat{E}}_#3"
            "+"
            "\\frac{1}{H_#1}"
            "\\pder{H_#2}{X^#1}"
            "\\vec{\\hat{E}}_#2"
        ),
        3
    ]
    macros["christoffel"] = (
        "\\Gamma_{ij}^k"
        "\\equiv"
        "\\sum_l"
        "\\pder{}{X^i}"
        "\\left("
        "\\pder{x^l}{X^j}"
        "\\right)"
        "\\pder{X^k}{x^l}"
    )
    macros["fromchristoffeltoscalefactor"] = (
        "2 H_k H_k \\Gamma_{ij}^k"
        "="
        "-"
        "\\pder{}{X^k}"
        "\\left( \\delta_{ij} H_i H_j \\right)"
        "+"
        "\\pder{}{X^i}"
        "\\left( \\delta_{jk} H_j H_k \\right)"
        "+"
        "\\pder{}{X^j}"
        "\\left( \\delta_{ki} H_k H_i \\right)"
    )
    macros["sumofchristoffel"] = (
        "\\sum_j"
        "\\Gamma_{ij}^j"
        "\\mathcal{Q}"
        "="
        "\\sum_j"
        "\\frac{1}{H_j}"
        "\\pder{H_j}{X^i}"
        "\\mathcal{Q}"
        "="
        "\\frac{1}{J}"
        "\\pder{J}{X^i}"
        "\\mathcal{Q}"
    )

