"""
    ApplicationSpecificConstants

    set of special constants for comparison with other libraries

    To use constants from library xyz, perform the destruction operator
        `( ;  kB, q, ε0, eV ) = xyz_constants()`
    in your scripts.
"""
@kwdef struct ApplicationSpecificConstants
    kB::Float64
    q::Float64
    ε0::Float64
    eV::Float64
end


# Set this as temporary values from pdelib
function pdelib_constants()
    return ApplicationSpecificConstants(
        kB = 1.3806503e-23,
        q = 1.602176462e-19,
        ε0 = 8.85418781762039e-12,
        eV = q * V,
    )
end

function set_TeSCA_constants()
    return ApplicationSpecificConstants(
        kB = 1.380662e-23,
        q = 1.6021e-19,
        ε0 = 8.85419e-12,
        eV = q * V
    )
end

# set unity for constants
function set_unity_constants()
    return ApplicationSpecificConstants(
        kB = 1.0,
        q = 1.0,
        ε0 = 1.0,
        eV = 1.0
    )
end
