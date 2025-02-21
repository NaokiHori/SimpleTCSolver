def add(macros):
    # coordinate
    macros["vx"] = "{x}"
    macros["vy"] = "{y}"
    macros["vr"] = "{r}"
    macros["vt"] = "{\\theta}"
    macros["vz"] = "{z}"
    # a symbol used to denote general coordinate
    macros["gcs"] = ["{\\xi^{#1}}", 1]
    # velocity and scalar
    macros["vel"] = ["{u_{#1}}", 1]
    macros["scalar"] = "{T}"
    # quadratic quantities
    macros["quad"] = ["{k_{#1}}", 1]
    # scale factors
    macros["sfact"] = ["{h_{\\gcs{#1}}}", 1]
    # velocity-gradient tensor and shear-stress tensor
    macros["vgt"] = ["{l_{#1 #2}}", 2]
    macros["sst"] = ["{\\tau_{#1 #2}}", 2]
    # derivatives
    macros["pder"] = ["{\\frac{\\partial #1}{\\partial #2}}", 2]
    macros["dder"] = ["{\\frac{\\delta #1}{\\delta #2}}", 2]
    macros["tder"] = ["\\frac{d #1}{d #2}", 2]
    macros["mder"] = ["\\frac{D #1}{D t}", 1]
    # discrete operators
    macros["ave"] = ["{\\overline{#1}^{#2}}", 2]
    macros["dif"] = ["{\\delta_{#2} {#1}}", 2]
    macros["vat"] = ["{\\left. {#1} \\right|_{#2}}", 2]
