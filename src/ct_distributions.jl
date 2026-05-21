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

### Approximation of Gauss-Fermi integral using Paasch's method [J. Appl. Phys. 107, 104501 (2010)]
### Paasch-Scheinert functions
"""
$(TYPEDSIGNATURES)

Paasch-Scheinert approximation of Gauss-Fermi integral parameterized by s = sigma/k_BT (the gaussian width in units of thermal energy).
For details on function see Paasch and Scheinert, J. Appl. Phys 107, 104501 (2010)
"""
struct GaussFermiPaasch{T} <: Function
    s::T
end
function (GaussFermiPaasch::GaussFermiPaasch{T})(x::Real) where {T}
    function H(s::Real)
        return sqrt(2) / s * erfcinv(exp(-(s^2) / 2))
    end
    function K(s::Real)
        return 2 * (1 - H(s) / s * sqrt(2 / pi) * exp(1 / 2 * s^2 * (1 - H(s)^2)))
    end
    s = GaussFermiPaasch.s
    if (abs(x) > s^2)
        G = exp((s * s / 2 - abs(x))) / (1.0 + exp(K(s) * (s * s - abs(x))))
    else
        G = 0.5 * erfc(abs(x) / (s * sqrt(2)) * H(s))
    end
    if (x > 0)
        G = 1 - G
    end
    return G
end
"""
$(TYPEDEF)

Struct containing information for numerical integration of Gauss-Fermi integral. 
Can be constructed automatically from ctsys.data using constructGaussFermiSimpson13(data, itrap::Int64).
"""
# Approximate Gauss-Fermi integral using Simpsons rule
struct GaussFermiSimpson13 <: Function
    s::Float64
    Et::Float64
    kBT::Float64
    nPoints::Int64
end
"""
$(TYPEDSIGNATURES)

Simpson's 1/3 quadrature rule approximation of Gauss-Fermi integral parameterized by s = sigma/k_BT (the gaussian width in units of thermal energy).
"""
function (GaussFermiSimpson13::GaussFermiSimpson13)(x::Real)

    @inline function gaussfermi_G(E::Float64, sigma::Float64, Et::Float64, kBT::Float64, x::Real)
        return exp(-(E - Et) * (E - Et) / (2 * sigma * sigma)) * FermiDiracMinusOne(x - (E - Et) / kBT)
    end

    s = GaussFermiSimpson13.s
    Et = GaussFermiSimpson13.Et
    kBT = GaussFermiSimpson13.kBT
    nP = Int64(2 * ceil((GaussFermiSimpson13.nPoints) / 2.0) + 1)

    Elower = Et - 6 * s
    Eupper = Et + 6 * s

    dE = (Eupper - Elower) / (nP - 1.0)

    v = (2.0 * sum(2.0 * gaussfermi_G(E, s, Et, kBT, x) + gaussfermi_G(E + dE, s, Et, kBT, x) for E in (Elower + dE):(2 * dE):(Eupper - 2 * dE)) + gaussfermi_G(Elower, s, Et, kBT, x) + gaussfermi_G(Eupper, s, Et, kBT, x)) * dE / (3.0 * sqrt(2.0π) * s)

    return v
end
"""
$(TYPEDSIGNATURES)


Plot different distribution integrals.
"""
function plotDistributions(; Plotter = nothing)

    Plotter.close()

    x = -5:0.1:10

    Plotter.semilogy(x, FermiDiracOneHalfBednarczyk.(x), label = "\$F_{1/2}  \$ (Bednarczyk)")
    Plotter.semilogy(x, FermiDiracOneHalfTeSCA.(x), label = "\$F_{1/2} \$ (TeSCA)")
    Plotter.semilogy(x, Boltzmann.(x), label = "Boltzmann")
    Plotter.semilogy(x, ones(size(x)) / 0.27, "--", label = "\$1/\\gamma=3.\\overline{703}\$", color = (0.6, 0.6, 0.6, 1))
    Plotter.semilogy(x, Blakemore.(x), label = "Blakemore (\$\\gamma=0.27\$)")
    Plotter.semilogy(x, degenerateLimit.(x), label = "degenerate limit")
    Plotter.semilogy(x, GaussFermiPaasch(1).(x), label = "Gauss-Fermi (Paasch), ŝ = 1")
    Plotter.semilogy(x, GaussFermiSimpson13(1, 0, 1, 1000).(x), label = "Gauss-Fermi (Simpson 1/3), ŝ = 1", linestyle = "dotted")
    Plotter.semilogy(x, GaussFermiPaasch(10).(x), label = "Gauss-Fermi (Paasch), ŝ = 10")
    Plotter.semilogy(x, GaussFermiSimpson13(10, 0, 1, 1000).(x), label = "Gauss-Fermi (Simpson 1/3), ŝ = 10", linestyle = "dotted")

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

    f = ChargeTransport.GaussFermiPaasch(1); df = x -> ForwardDiff.derivative(f, x)
    Plotter.semilogy(x, f.(x) ./ df.(x), label = "Gauss-Fermi (Paasch), ŝ = 1")

    f = ChargeTransport.GaussFermiSimpson13(1, 0, 1, 1000); df = x -> ForwardDiff.derivative(f, x)
    Plotter.semilogy(x, f.(x) ./ df.(x), label = "Gauss-Fermi (Simpson 1/3), ŝ = 1", linestyle = "dotted")

    f = ChargeTransport.GaussFermiPaasch(10); df = x -> ForwardDiff.derivative(f, x)
    Plotter.semilogy(x, f.(x) ./ df.(x), label = "Gauss-Fermi (Paasch), ŝ = 10")

    f = ChargeTransport.GaussFermiSimpson13(10, 0, 1, 1000); df = x -> ForwardDiff.derivative(f, x)
    Plotter.semilogy(x, f.(x) ./ df.(x), label = "Gauss-Fermi (Simpson 1/3), ŝ = 10", linestyle = "dotted")

    Plotter.xlabel("\$\\eta\$")
    Plotter.ylabel("\$g(\\eta)\$")
    Plotter.title("Diffusion Enhancements")
    Plotter.legend()
    Plotter.grid()

    return Plotter.show()
end
