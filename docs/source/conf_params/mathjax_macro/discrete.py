def add(macros):
    # number of grid points
    macros["ngp"] = ["N_{#1}", 1]
    # summation symbols for different locations
    macros["sumxf"] = "\\sum_{i = \\frac{1}{2}}^{\\ngp{1} + \\frac{1}{2}}"
    macros["sumxc"] = "\\sum_{i = 1}^{\\ngp{1}}"
    macros["sumyf"] = "\\sum_{j = \\frac{1}{2}}^{\\ngp{2} - \\frac{1}{2}}"
    macros["sumyc"] = "\\sum_{j = 1}^{\\ngp{2}}"
    macros["sumzf"] = "\\sum_{k = \\frac{1}{2}}^{\\ngp{3} - \\frac{1}{2}}"
    macros["sumzc"] = "\\sum_{k = 1}^{\\ngp{3}}"
    # indices
    macros["cmindex"] = ["#1-\\frac{1}{2}", 1]
    macros["ccindex"] = ["#1             ", 1]
    macros["cpindex"] = ["#1+\\frac{1}{2}", 1]
    # velocity gradient
    macros["dlxx"] = (
            "\\frac{1}{\\sfact{1}}"
            "\\dif{\\vel{1}}{\\gcs{1}}"
    )
    macros["dlyx"] = (
            "\\dlyxa"
            "+"
            "\\dlyxb"
    )
    macros["dlyxa"] = (
            "\\frac{1}{\\sfact{2}}"
            "\\dif{\\vel{1}}{\\gcs{2}}"
    )
    macros["dlyxb"] = (
            "\\left\\{"
            "  \\begin{alignat}{2}"
            "    &"
            "    \\text{Negative wall:}"
            "    &"
            "    -"
            "    \\frac{1}{J}"
            "    \\dif{\\left( \\frac{J}{\\sfact{1}} \\right)}{\\gcs{1}}"
            "    \\frac{\\vat{\\vel{2}}{1, \\cpindex{j}, \\ccindex{k}}}{2}"
            "    \\\\"
            "    &"
            "    \\text{Positive wall:}"
            "    &"
            "    -"
            "    \\frac{1}{J}"
            "    \\dif{\\left( \\frac{J}{\\sfact{1}} \\right)}{\\gcs{1}}"
            "    \\frac{\\vat{\\vel{2}}{\\ngp{1}, \\cpindex{j}, \\ccindex{k}}}{2}"
            "    \\\\"
            "    &"
            "    \\text{Otherwise:}"
            "    &"
            "    -"
            "    \\frac{1}{J}"
            "    \\dif{\\left( \\frac{J}{\\sfact{1}} \\right)}{\\gcs{1}}"
            "    \\ave{\\vel{2}}{\\gcs{1}}"
            "    \\\\"
            "  \\end{alignat}"
            "\\right."
    )
    macros["dlzx"] = (
            "\\frac{1}{\\sfact{3}}"
            "\\dif{\\vel{1}}{\\gcs{3}}"
    )
    macros["dlxy"] = (
            "\\frac{1}{\\sfact{1}}"
            "\\dif{\\vel{2}}{\\gcs{1}}"
    )
    macros["dlyy"] = (
            "\\dlyya"
            "+"
            "\\dlyyb"
    )
    macros["dlyya"] = (
            "\\frac{1}{\\sfact{2}}"
            "\\dif{\\vel{2}}{\\gcs{2}}"
    )
    macros["dlyyb"] = (
            "\\frac{1}{J}"
            "\\dif{\\left( \\frac{J}{\\sfact{1}} \\right)}{\\gcs{1}}"
            "\\ave{\\vel{1}}{\\gcs{1}}"
    )
    macros["dlzy"] = (
            "\\frac{1}{\\sfact{3}}"
            "\\dif{\\vel{2}}{\\gcs{3}}"
    )
    macros["dlxz"] = (
            "\\frac{1}{\\sfact{1}}"
            "\\dif{\\vel{3}}{\\gcs{1}}"
    )
    macros["dlyz"] = (
            "\\frac{1}{\\sfact{2}}"
            "\\dif{\\vel{3}}{\\gcs{2}}"
    )
    macros["dlzz"] = (
            "\\frac{1}{\\sfact{3}}"
            "\\dif{\\vel{3}}{\\gcs{3}}"
    )
    # discrete incompressibility
    macros["ddiv"] = [
            "\\frac{1}{J}"
            "\\dif{}{\\gcs{1}}"
            "\\left("
            "   \\frac{J}{\\sfact{1}}"
            "   \\vel{1}"
            "\\right)"
            "+"
            "\\frac{1}{J}"
            "\\dif{}{\\gcs{2}}"
            "\\left("
            "   \\frac{J}{\\sfact{2}}"
            "   \\vel{2}"
            "\\right)"
            "+"
            "\\frac{1}{J}"
            "\\dif{}{\\gcs{3}}"
            "\\left("
            "   \\frac{J}{\\sfact{3}}"
            "   \\vel{3}"
            "\\right)"
            "="
            "0"
    ]
    # discrete advective terms in momentum balance
    macros["dmomadv"] = [
            "-"
            "\\frac{1}{J}"
            "\\ave{"
            "  \\ave{"
            "    \\frac{J}{\\sfact{#2}}"
            "    \\vel{#2}"
            "  }{\\gcs{#1}}"
            "  \\dif{\\vel{#1}}{\\gcs{#2}}"
            "}{\\gcs{#2}}"
            , 2
    ]
    # centrifugal term
    macros["dmomadvx"] = (
            "+"
            "\\frac{1}{J}"
            "\\ave{"
            "  \\dif{"
            "    \\left("
            "      \\frac{J}{\\sfact{1}}"
            "    \\right)"
            "  }{\\gcs{1}}"
            "  \\ave{\\vel{2}}{\\gcs{2}}"
            "  \\ave{\\vel{2}}{\\gcs{2}}"
            "}{\\gcs{1}}"
    )
    # coriolis term
    macros["dmomadvy"] = (
            "-"
            "\\frac{1}{J}"
            "\\ave{"
            "  \\dif{"
            "    \\left("
            "      \\frac{J}{\\sfact{1}}"
            "    \\right)"
            "  }{\\gcs{1}}"
            "  \\ave{\\vel{2}}{\\gcs{2}}"
            "  \\ave{\\vel{1}}{\\gcs{1}}"
            "}{\\gcs{2}}"
    )
    # discrete pressure-gradient terms in momentum balance
    macros["dmompre"] = [
            "-"
            "\\frac{1}{\\sfact{#1}}"
            "\\dif{p}{\\gcs{#1}}"
            , 1
    ]
    # discrete diffusive terms in momentum balance
    macros["dmomdif"] = [
            "+"
            "\\frac{1}{J}"
            "\\dif{}{\\gcs{#1}}"
            "\\left("
            "   \\frac{J}{\\sfact{#1}}"
            "   \\sst{#1}{#2}"
            "\\right)"
            , 2
    ]
    # additional diffusive term in x
    macros["dmomdifx"] = (
            "-"
            "\\frac{1}{J}"
            "\\ave{"
            "  \\dif{"
            "    \\left("
            "      \\frac{J}{\\sfact{1}}"
            "    \\right)"
            "  }{\\gcs{1}}"
            "  \\sst{2}{2}"
            "}{\\gcs{1}}"
    )
    # additional diffusive term in y
    macros["dmomdify"] = (
            "+"
            "\\frac{1}{J}"
            "\\ave{"
            "  \\dif{"
            "    \\left("
            "      \\frac{J}{\\sfact{1}}"
            "    \\right)"
            "  }{\\gcs{1}}"
            "  \\sst{2}{1}"
            "}{\\gcs{1}}"
    )

    # discrete advective terms in scalar transport
    macros["dscalaradv"] = [
            "-"
            "\\frac{1}{J}"
            "\\ave{"
            "  \\frac{J}{\\sfact{#1}}"
            "  \\vel{#1}"
            "  \\dif{\\scalar}{\\gcs{#1}}"
            "}{\\gcs{#1}}"
            , 1
    ]
    # discrete diffusive terms in scalar transport
    macros["dscalardif"] = [
            "+"
            "\\kappa"
            "\\frac{1}{J}"
            "\\dif{"
            "}{\\gcs{#1}}"
            "\\left("
            "  \\frac{J}{\\sfact{#1}}"
            "  \\frac{1}{\\sfact{#1}}"
            "  \\dif{\\scalar}{\\gcs{#1}}"
            "\\right)"
            , 1
    ]
