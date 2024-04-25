def add(macros):
    # temporal derivatives
    macros["scalartemp"] = (
            "\\pder{\\scalar}{t}"
    )

    # conservative advective terms
    macros["scalaradv"] = [
            (
                "-"
                "\\frac{\\vel{#1}}{\\sfact{#1}}"
                "\\pder{\\scalar}{\\gcs{#1}}"
            ),
            1
    ]

    # conservative diffusive terms
    macros["scalardif"] = [
            (
                "+"
                "\\kappa"
                "\\frac{1}{J}"
                "\\pder{}{\gcs{#1}}"
                "\\left("
                "   \\frac{"
                "      J"
                "   }{"
                "      \\sfact{#1}"
                "   }"
                "   \\frac{1}{\\sfact{#1}}"
                "   \\pder{\\scalar}{\\gcs{#1}}"
                "\\right)"
            ),
            1
    ]

