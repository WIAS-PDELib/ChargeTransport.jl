using ForwardDiff
using SpecialFunctions
"""
$(TYPEDSIGNATURES)


The Boltzmann statistics function ``\\exp(x)``.
"""
function Boltzmann(x::Real)
    return exp(x)
end

"""
$(TYPEDSIGNATURES)


The Blakemore approximation ``1/(\\exp(-x) + γ)`` with variable real scalar ``γ``, see
J. S. Blakemore. “The Parameters of Partially Degenerate Semiconductors”. In: Proceedings of
the Physical Society. Section A 65 (1952), pp. 460–461.

"""
function Blakemore(x::Real, γ::Real)
    return 1 / (exp(-x) + γ)
end


# The Blakemore approximation ``1/(\\exp(-x) + γ)`` with ``γ = 0.27``
function Blakemore(x::Real)
    return Blakemore(x, 0.27)
end


"""
$(TYPEDSIGNATURES)


The Fermi-Dirac integral of order ``-1`` which reads ``1/(\\exp(-x) + 1)``, see J.S. Blakemore,
Approximations for Fermi-Dirac integrals, especially the function ``F_{1/2} (\\eta)`` used to
describe electron density in a semiconductor, Solid-State Electronics 25 (11) (1982) 1067 – 1076.
"""
function FermiDiracMinusOne(x::Real)
    return Blakemore(x, 1.0)
end


"""
$(TYPEDSIGNATURES)


The incomplete Fermi-Dirac integral of order 1/2, implemented according to [Bednarczyk1978,
"The Approximation of the Fermi-Dirac integral ``F_{1/2}(\\eta)``"].
"""
function FermiDiracOneHalfBednarczyk(x::Real)

    a = x^4 + 33.6 * x * (1.0 - 0.68 * exp(-0.17 * (x + 1)^2)) + 50
    return 1.0 / (3 / 4 * sqrt(pi) * a^(-3 / 8) + exp(-x))

end

"""
$(TYPEDSIGNATURES)


The incomplete Fermi-Dirac integral of order 1/2, implemented according to the software
package TeSCA, see https://wias-berlin.de/software/index.jsp?lang=1&id=TeSCA.

Modified to use log1p(x)=log(1+x).
"""
function FermiDiracOneHalfTeSCA(x::Real)
    if x < 1.6107
        z = log1p(exp(x))
        return (1 + 0.16 * z) * z
    elseif 1.6107 <= x <= 344.7
        z = log1p(exp(x^(3 / 4)))
        return 0.3258 - (0.0321 - 0.7523 * z) * z
    else
        z = x^(3 / 4)
        return 0.3258 - (0.0321 - 0.7523 * z) * z
    end
end


"""
$(TYPEDSIGNATURES)


Degenerate limit of incomplete Fermi-Dirac integral of order 1/2.
"""
function degenerateLimit(x)
    return x < 0 ? NaN : 4 / (3 * sqrt(pi)) * x^(3 / 2)
end

function gaussFermi_sHatEquals2(x)
    sHat=2
    H = 0.7466947053286296
    K = 0.5565965819658305
    xp(x) = abs(x)
    if(xp(x) > sHat^2)
        G = exp((sHat*sHat/2-xp(x))) / (1.0 + exp(K * (sHat*sHat-xp(x)) ))
    else
        G = 0.5 * erfc(xp(x)/(sHat*sqrt(2)) * H)
    end
    # G = xp > 4 ? exp(-abs(2-xp)) * 1.0 / (1.0 + exp(-0.5565965819658305 * abs(-xp + 4.0))) : 0.5 * erfc(xp/(2sqrt(2)) * 0.7466947053286296)
    if(x>0)
        G = 1-G
    end
    return G
end

function gaussFermi_sHatEquals4(x)
    sHat=4
    H = 0.896563421365437
    K = 0.28181676579245685
    xp(x) = abs(x)
    if(xp(x) > sHat^2)
        G = exp((sHat*sHat/2-xp(x))) / (1.0 + exp(K * (sHat*sHat-xp(x)) ))
    else
        G = 0.5 * erfc(xp(x)/(sHat*sqrt(2)) * H)
    end
    if(x>0)
        G = 1-G
    end
    return G
end

function gaussFermi_sHatEquals6(x)
    sHat=6
    H = 0.9431587310721202
    K = 0.1683730458855388
    xp(x) = abs(x)
    if(xp(x) > sHat^2)
        G = exp((sHat*sHat/2-xp(x))) / (1.0 + exp(K * (sHat*sHat-xp(x)) ))
    else
        G = 0.5 * erfc(xp(x)/(sHat*sqrt(2)) * H)
    end
    if(x>0)
        G = 1-G
    end
    return G
end


function gaussFermi_sHatEquals8(x)
    sHat=8
    H = 0.9636458688943249
    K = 0.11249596843700593
    xp(x) = abs(x)
    if(xp(x) > sHat^2)
        G = exp((sHat*sHat/2-xp(x))) / (1.0 + exp(K * (sHat*sHat-xp(x)) ))
    else
        G = 0.5 * erfc(xp(x)/(sHat*sqrt(2)) * H)
    end
    if(x>0)
        G = 1-G
    end
    return G
end


function gaussFermi_sHatEquals10(x)
    sHat=10
    H = 0.9745474864965935
    K = 0.08091495970175866
    xp(x) = abs(x)
    if(xp(x) > sHat^2)
        G = exp((sHat*sHat/2-xp(x))) / (1.0 + exp(K * (sHat*sHat-xp(x)) ))
    else
        G = 0.5 * erfc(xp(x)/(sHat*sqrt(2)) * H)
    end
    if(x>0)
        G = 1-G
    end
    return G
end


function gaussFermi_sHatEquals30(x)
    sHat=30
    H = 0.9959651235434415
    K = 0.013889673859928786
    xp(x) = abs(x)
    if(xp(x) > sHat^2)
        G = exp((sHat*sHat/2-xp(x))) / (1.0 + exp(K * (sHat*sHat-xp(x)) ))
    else
        G = 0.5 * erfc(xp(x)/(sHat*sqrt(2)) * H)
    end
    if(x>0)
        G = 1-G
    end
    return G
end

"""
$(TYPEDSIGNATURES)


Plot different distribution integrals.
"""
function plotDistributions(; Plotter = nothing)

    Plotter.close()

    x = -5:0.1:700

    Plotter.semilogy(x, FermiDiracOneHalfBednarczyk.(x), label = "\$F_{1/2}  \$ (Bednarczyk)")
    Plotter.semilogy(x, FermiDiracOneHalfTeSCA.(x), label = "\$F_{1/2} \$ (TeSCA)")
    Plotter.semilogy(x, Boltzmann.(x), label = "Boltzmann")
    Plotter.semilogy(x, ones(size(x)) / 0.27, "--", label = "\$1/\\gamma=3.\\overline{703}\$", color = (0.6, 0.6, 0.6, 1))
    Plotter.semilogy(x, Blakemore.(x), label = "Blakemore (\$\\gamma=0.27\$)")
    Plotter.semilogy(x, degenerateLimit.(x), label = "degenerate limit")

    Plotter.xlabel("\$\\eta\$")
    Plotter.ylabel("\$\\mathcal{F}(\\eta)\$")
    Plotter.title("Distributions")
    Plotter.legend()
    Plotter.grid()

    return Plotter.show()
end

"""
$(TYPEDSIGNATURES)


Plot diffusion enhancements.
"""
function plotDiffusionEnhancements(; Plotter = nothing)

    Plotter.close()

    x = -5:0.01:10

    f = ChargeTransport.FermiDiracOneHalfBednarczyk; df = x -> ForwardDiff.derivative(f, x)
    Plotter.semilogy(x, f.(x) ./ df.(x), label = "\$F_{1/2}\$")

    f = ChargeTransport.Boltzmann; df = x -> ForwardDiff.derivative(f, x)
    Plotter.semilogy(x, f.(x) ./ df.(x), label = "Boltzmann")

    f = ChargeTransport.Blakemore; df = x -> ForwardDiff.derivative(f, x)
    Plotter.semilogy(x, f.(x) ./ df.(x), label = "Blakemore (\$\\gamma=0.27\$)")

    f = ChargeTransport.degenerateLimit; df = x -> ForwardDiff.derivative(f, x)
    Plotter.semilogy(x, f.(x) ./ df.(x), label = "degenerate limit")

    Plotter.xlabel("\$\\eta\$")
    Plotter.ylabel("\$g(\\eta)\$")
    Plotter.title("Diffusion Enhancements")
    Plotter.legend()
    Plotter.grid()

    return Plotter.show()
end
