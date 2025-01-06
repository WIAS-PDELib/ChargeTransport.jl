using ForwardDiff

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
    return sqrt(pi) / (2 * (3 / 4 * sqrt(pi) * a^(-3 / 8) + exp(-x)))
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

    Plotter.show()
    return nothing
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

    Plotter.show()
    return nothing
end
