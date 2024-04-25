def add(macros):
    # temporal derivatives
    macros["momtemp"] = [
            (
                "\\pder{\\vel{#1}}{t}"
            ),
            1
    ]

    # conservative advective terms
    # 1st argument: direction to which transported
    # 2nd argument: transported velocity
    macros["momadv"] = [
            (
                "-"
                "\\frac{\\vel{#1}}{\\sfact{#1}}"
                "\\pder{\\vel{#2}}{\\gcs{#1}}"
            ),
            2
    ]

    # additional advective term in x
    macros["momadvx"] = (
            "+"
            "\\frac{1}{J}"
            "\\pder{}{\\gcs{1}}"
            "\\left("
            "  \\frac{J}{\\sfact{1}}"
            "\\right)"
            "\\vel{2}"
            "\\vel{2}"
    )

    # additional advective term in y
    macros["momadvy"] = (
            "-"
            "\\frac{1}{J}"
            "\\pder{}{\\gcs{1}}"
            "\\left("
            "  \\frac{J}{\\sfact{1}}"
            "\\right)"
            "\\vel{2}"
            "\\vel{1}"
    )

    # pressure-gradient terms
    macros["mompre"] = [
            (
                "-"
                "\\frac{1}{h_{\\gcs{#1}}}"
                "\\pder{p}{\\gcs{#1}}"
            ),
            1
    ]

    # conservative diffusive terms
    # 1st argument: direction to which transported
    # 2nd argument: transported velocity
    macros["momdif"] = [
            (
                "+"
                "\\frac{1}{J}"
                "\\pder{}{\gcs{#1}}"
                "\\left("
                "   \\frac{"
                "      J"
                "   }{"
                "      \\sfact{#1}"
                "   }"
                "   \\sst{#1}{#2}"
                "\\right)"
            ),
            2
    ]

    # additional diffusive term in x
    macros["momdifx"] = [
            (
                "-"
                "\\frac{1}{J}"
                "\\pder{}{\\gcs{1}}"
                "\\left("
                "  \\frac{J}{\\sfact{1}}"
                "\\right)"
                "\\sst{2}{2}"
            )
    ]

    # additional diffusive term in y
    macros["momdify"] = [
            (
                "+"
                "\\frac{1}{J}"
                "\\pder{}{\\gcs{1}}"
                "\\left("
                "  \\frac{J}{\\sfact{1}}"
                "\\right)"
                "\\sst{2}{1}"
            )
    ]

