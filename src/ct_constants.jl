"""
    ApplicationSpecificConstants

    set of special constants for comparison with other libraries

    To use constants from library xyz, perform the destruction operator
        `( ;  k_B, q, ε_0, eV ) = xyz_constants()`
    in your scripts.
"""
@kwdef struct ApplicationSpecificConstants
    k_B::Float64  # BoltzmannConstant
    q::Float64   # ElementaryCharge
    ε_0::Float64 # VacuumElectricPermittivity
    eV::Float64
end


# Set this as temporary values from pdelib
function pdelib_constants()
    @local_unitfactors V
    return ApplicationSpecificConstants(
        k_B = 1.3806503e-23,
        q = 1.602176462e-19,
        ε_0 = 8.85418781762039e-12,
        eV = q * V,
    )
end

function set_TeSCA_constants()
    @local_unitfactors V
    return ApplicationSpecificConstants(
        k_B = 1.380662e-23,
        q = 1.6021e-19,
        ε_0 = 8.85419e-12,
        eV = q * V
    )
end

# set unity for constants
function set_unity_constants()
    return ApplicationSpecificConstants(
        k_B = 1.0,
        q = 1.0,
        ε_0 = 1.0,
        eV = 1.0
    )
end

# Numerical parameters
const tiny_penalty_value = 1.0e-10        # tiny penalty value
