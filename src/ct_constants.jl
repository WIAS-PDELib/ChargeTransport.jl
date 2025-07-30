"""
    Constants

    Default physical constants (dimensionless) from PhysicalConstants.jl via LessUnitful.jl
"""
@kwdef struct Constants
    k_B = ph"BoltzmannConstant"          # (k_B) JK^{-1}
    Planck_constant = ph"PlanckConstant"  # (h)   Js
    m_e = ph"ElectronMass"               # (m_e) kg
    q = ph"ElementaryCharge"             # (e)   C
    ε_0 = ph"VacuumElectricPermittivity" # (ε_0) C/(V*m)
end

"""
    constants

    A globally available object containing the default constants
"""
const constants = Constants()


"""
    pdelib_constants

    constants with slightly modified values as used in pdelib (https://wias-berlin.de/software/index.jsp?id=pdelib)
"""
const pdelib_constants = Constants(
    k_B = 1.3806503e-23,
    q = 1.602176462e-19,
    ε_0 = 8.85418781762039e-12
)

"""
    teSCA_constants

    constants with slightly modified values used in WIAS-TeSCA (https://wias-berlin.de/software/index.jsp?id=TeSCA&lang=0&archive=true)
"""
const teSCA_constants = Constants(
    k_B = 1.380662e-23,
    q = 1.6021e-19,
    ε_0 = 8.85419e-12
)

"""
    unity_constants

    unit constants
"""
const unity_constants = Constants(
    k_B = 1.0,
    Planck_constant = 1.0,
    m_e = 1.0,
    q = 1.0,
    ε_0 = 1.0
)

# Numerical parameters
const tiny_penalty_value = 1.0e-10        # tiny penalty value
