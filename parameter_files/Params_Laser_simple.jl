# Parameters for a simple laser simulation
@kwdef struct Params_Laser_simple

    #####################################################################
    ############################ parameters ############################

    ########## charge carriers ##########
    iphin = 1 # electron quasi Fermi potential
    iphip = 2 # hole quasi Fermi potential
    numberOfCarriers = 2

    #general parameters
    T = 300.0 * K             # temperature
    λ = 0.88e-6 * m             # wavelength   =880 nm
    U0 = 1.41 * V             # first maxBias for first solution
    U = [1.45, 1.49, 1.53, 1.57, 1.61, 1.65, 1.69, 1.73, 1.77, 1.81] * V   # further bias values
    P = [0.0, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0] * 1.0e-3 * W   # nu=1:1, P_1 = 300 mW
    speedoflight = 299_792_458 * m / s

    ########## device geometry ##########

    # region numbers
    regionDonor1 = 1               # n doped regions
    regionDonor2 = 2
    regionDonor3 = 3
    regionAcceptor1 = 4               # p doped regions
    regionAcceptor2 = 5
    regionsDonor = [regionDonor1, regionDonor2, regionDonor3]
    regionsAcceptor = [regionAcceptor1, regionAcceptor2]

    # boundary region numbers
    bregionDonor1 = 1    # bottom boundary
    bregionAcceptor2 = 2    # top boundary
    bregionNoFlux = 3

    # contact voltage
    topVoltageAcceptor2 = U0


    ########## laser model parameters ##########

    # refractive index parameters
    nTilde_d0 = 0.0e-23 * 1.0e-6     # =0 means no nonlinearity, converges quickly
    # nTilde_d0    =  7.0e-23*1e-6
    γ_n0 = 0.5

    #gain parameters
    g_00 = 3.0e+5

    #absorption parameters
    α_00 = 0.0e2
    # α_00   =  30e2
    f_cn0 = 4.0e-22
    f_cp0 = 1.2e-21

    #recombination parameters
    B_rad0 = 3.0e-16
    C_augn0 = 1.0e-42
    C_augp0 = 1.0e-42
    τ_n0 = 1.4e-9
    τ_p0 = 1.4e-9

    #overall mobility
    μn0 = 1000.0e-4
    μp0 = 100.0e-4


    ##                parameters for every material
    ##               1       2       3       4       5

    doping = [1.0, 0.1, 0.05, +0.1, +1.0] * 1.0e+24

    #relative dielectric permittivity
    εr = [11.48, 12.02, 12.85, 12.02, 11.48]

    ## band edge energies
    EV = [-1.72, -1.602, -1.415, -1.602, -1.72] * eV         # Ev = Ec - Eg
    EC = [0.302, 0.243, -0.009, 0.243, 0.302] * eV
    #EG         = [  2.022,  1.845,  1.406,  1.845,  2.022 ] *  eV

    ## effective densities of states
    NC = [13.83, 0.77, 0.47, 0.77, 13.83] * 1.0e+24 / (m^3)
    NV = [15.28, 14.23, 10.0, 14.23, 15.28] * 1.0e+24 / (m^3)

    ## recombination parameters
    r0 = [1.0, 1.0, 1.0, 1.0, 1.0] * B_rad0
    Auger_Cn = [1.0, 1.0, 1.0, 1.0, 1.0] * C_augn0
    Auger_Cp = [1.0, 1.0, 1.0, 1.0, 1.0] * C_augp0
    τn = [1.0, 1.0, 1.0, 1.0, 1.0] * τ_n0
    τp = [1.0, 1.0, 1.0, 1.0, 1.0] * τ_p0

    ## mobilities
    μn = [0.101, 1.191, 5.108, 1.191, 0.101] * μn0 * (m^2) / (V * s)
    μp = [0.492, 1.225, 3.037, 1.225, 0.492] * μp0 * (m^2) / (V * s)


    ## refractive index parameters
    nTilde = [3.258, 3.382, 3.628, 3.382, 3.258]                                  # n_0
    nTilde_d = [0.0, 0.0, 1.0, 0.0, 0.0] * nTilde_d0                      # n_d
    γn = [1.0, 1.0, 1.0, 1.0, 1.0] * γ_n0                           # γ_n

    ## gain parameters
    gain0 = [0.0, 0.0, 1.0, 0.0, 0.0] * g_00
    α0 = [0.0, 0.0, 1.0, 0.0, 0.0] * α_00
    fcnalf = [1.0, 1.0, 1.0, 1.0, 1.0] * f_cn0
    fcpalf = [1.0, 1.0, 1.0, 1.0, 1.0] * f_cp0


    #interior doping, in donor and acceptor regions
    Nd = vcat(doping[1:3], [0, 0])
    Na = vcat([0, 0, 0], doping[4:5])

    #layer thicknesses
    heights = [1.0, 0.5, 0.05, 0.5, 1.0] * 1.0e-6           #vertical x5
    #width of lateral regions
    widths = [10.0, 2.0, 10.0] * 1.0e-6                       #horizontal

    heightPoints = [1.0, 1.5, 1.55, 2.05, 3.05] * 1.0e-6
    widthPoints = [1.0, 1.2, 2.2] * 1.0e-5


    #material in every zone from bottom to top
    mat = [
        1        1        1 ;
        2        2        2 ;
        3        3        3 ;
        4        4        4 ;
        0        5        0
    ]

end
