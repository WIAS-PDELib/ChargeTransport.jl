
module PSC_InterfaceSpecies

using VoronoiFVM
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot
using DelimitedFiles

function main(;n = 13, Plotter = PyPlot, plotting = false, verbose = false, test = false, unknown_storage=:dense)

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## region numbers
    regionAcceptor          = 1                           # p doped region
    regionIntrinsic         = 2                           # intrinsic region
    regionDonor             = 3                           # n doped region
    regions                 = [regionAcceptor, regionIntrinsic, regionDonor]
    numberOfRegions         = length(regions)

    ## boundary region numbers
    bregionAcceptor         = 1
    bregionDonor            = 2
    bregionJunction1        = 3
    bregionJunction2        = 4
    bregions                = [bregionAcceptor, bregionDonor, bregionJunction1, bregionJunction2]
    numberOfBoundaryRegions = length(bregions)

    ## grid (the nearer to interface, the finer)
    h_pdoping               = 3.00e-6 * cm + 1.0e-7 * cm
    h_intrinsic             = 3.00e-5 * cm 
    h_ndoping               = 8.50e-6 * cm + 1.0e-7 * cm

    x0                      = 0.0 * cm 
    δ                       = 4*n        # the larger, the finer the mesh
    t                       = 0.5*(cm)/δ # tolerance for geomspace and glue (with factor 10)
    k                       = 1.5        # the closer to 1, the closer to the boundary geomspace works

    coord_p_u               = collect(range(x0, h_pdoping/2, step=h_pdoping/(0.9*δ)))
    coord_p_g               = geomspace(h_pdoping/2, 
                                        h_pdoping, 
                                        h_pdoping/(1.2*δ), 
                                        h_pdoping/(1.2*δ), 
                                        tol=t)
    coord_i_g1              = geomspace(h_pdoping, 
                                        h_pdoping+h_intrinsic/k, 
                                        h_intrinsic/(7.1*δ), 
                                        h_intrinsic/(7.1*δ), 
                                        tol=t)
    coord_i_g2              = geomspace(h_pdoping+h_intrinsic/k, 
                                        h_pdoping+h_intrinsic,               
                                        h_intrinsic/(7.1*δ),    
                                        h_intrinsic/(7.1*δ), 
                                        tol=t)
    coord_n_g               = geomspace(h_pdoping+h_intrinsic,               
                                        h_pdoping+h_intrinsic+h_ndoping/2, 
                                        h_ndoping/(2.0*δ),   
                                        h_ndoping/(2.0*δ),      
                                        tol=t)
    coord_n_u               = collect(range(h_pdoping+h_intrinsic+h_ndoping/2, h_pdoping+h_intrinsic+h_ndoping, step=h_pdoping/(1.0*δ)))

    coord                   = glue(coord_p_u,coord_p_g,  tol=10*t) 
    icoord_p                = length(coord)
    coord                   = glue(coord,    coord_i_g1, tol=10*t)
    coord                   = glue(coord,    coord_i_g2, tol=10*t) 
    icoord_pi               = length(coord)
    coord                   = glue(coord,    coord_n_g,  tol=10*t)
    coord                   = glue(coord,    coord_n_u,  tol=10*t)
    grid                    = ExtendableGrids.simplexgrid(coord)
    numberOfNodes           = length(coord)

    ## cellmask! for defining the subregions and assigning region number (doping profiles do not intersect)
    cellmask!(grid, [0.0 * μm],                 [h_pdoping],                           regionAcceptor, tol = 1.0e-18)   # p-doped region   = 1
    cellmask!(grid, [h_pdoping],                [h_pdoping + h_intrinsic],             regionIntrinsic, tol = 1.0e-18)  # intrinsic region = 2
    cellmask!(grid, [h_pdoping + h_intrinsic],  [h_pdoping + h_intrinsic + h_ndoping], regionDonor, tol = 1.0e-18)      # n-doped region   = 3

    ## bfacemask! for ``active'' boundary regions, i.e. internal interfaces. On the outer boundary regions, the 
    ## conditions will be formulated later
    bfacemask!(grid, [h_pdoping],               [h_pdoping],                           bregionJunction1)  # first  inner interface
    bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic],             bregionJunction2)  # second inner interface

    # if plotting
    #     GridVisualize.gridplot(grid, Plotter = Plotter)
    #     Plotter.title("Grid")
    #     Plotter.figure()
    # end
 
    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    iphin               = 1 # electron quasi Fermi potential
    iphip               = 2 # hole quasi Fermi potential
    iphin_b1            = 3
    iphip_b1            = 4
    ipsi                = 4

    ichargeCarriers     = [iphin, iphip]   # this is an Array of indices of the charge carriers 
    numberOfCarriers    = length(ichargeCarriers) # electrons, holes and anion vacancies

    ## temperature
    T                   =  300.0                *  K

    ## band edge energies    
    Ec_a                = -3.0                  *  eV 
    Ev_a                = -5.1                  *  eV 

    Ec_i                = -3.8                  *  eV 
    Ev_i                = -5.4                  *  eV 

    Ec_d                = -3.8                  *  eV 
    Ev_d                = -6.2                  *  eV 

    EC                  = [Ec_a, Ec_i, Ec_d] 
    EV                  = [Ev_a, Ev_i, Ev_d]
    
    ## effective densities of state
    Nc_a                = 1.0e20                / (cm^3)
    Nv_a                = 1.0e20                / (cm^3)

    Nc_i                = 1.0e19                / (cm^3)
    Nv_i                = 1.0e19                / (cm^3)

    Nc_d                = 1.0e19                / (cm^3)
    Nv_d                = 1.0e19                / (cm^3)

    NC                  = [Nc_a, Nc_i, Nc_d]
    NV                  = [Nv_a, Nv_i, Nv_d]

    ## mobilities 
    μn_a                = 0.1                   * (cm^2) / (V * s)  
    μp_a                = 0.1                   * (cm^2) / (V * s)  

    μn_i                = 2.00e1                * (cm^2) / (V * s)  
    μp_i                = 2.00e1                * (cm^2) / (V * s)

    μn_d                = 1.0e-3                * (cm^2) / (V * s) 
    μp_d                = 1.0e-3                * (cm^2) / (V * s) 

    μn                  = [μn_a, μn_i, μn_d] 
    μp                  = [μp_a, μp_i, μp_d] 

    ## relative dielectric permittivity  
    ε_a                 = 4.0                   *  1.0  
    ε_i                 = 23.0                  *  1.0 
    ε_d                 = 3.0                   *  1.0 

    ε                   = [ε_a, ε_i, ε_d] 

    ## radiative recombination
    r0_a                = 6.3e-11               * cm^3 / s 
    r0_i                = 3.6e-12               * cm^3 / s  
    r0_d                = 6.8e-11               * cm^3 / s
        
    r0                  = [r0_a, r0_i, r0_d]
        
    ## life times and trap densities 
    τn_a                = 1.0e-6                * s 
    τp_a                = 1.0e-6                * s
        
    τn_i                = 1.0e-7                * s
    τp_i                = 1.0e-7                * s
    τn_d                = τn_a
    τp_d                = τp_a
        
    τn                  = [τn_a, τn_i, τn_d]
    τp                  = [τp_a, τp_i, τp_d]
        
    ## SRH trap energies (needed for calculation of trap_density! (SRH))
    Ei_a                = -4.05                 * eV   
    Ei_i                = -4.60                 * eV   
    Ei_d                = -5.00                 * eV   

    EI                  = [Ei_a, Ei_i, Ei_d]
        
    ## doping
    Nd                  = 2.089649130192123e17  / (cm^3) 
    Na                  = 4.529587947185444e18  / (cm^3) 
    C0                  = 1.0e18                / (cm^3) 

    ## contact voltages
    voltageAcceptor     =  1.0                  * V 

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## initialize Data instance and fill in data
    data                                = Data(grid, numberOfCarriers)

    ## possible choices: model_stationary, model_transient
    data.model_type                     = model_stationary

    ## possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA FermiDiracMinusOne, Blakemore
    data.F                              = [Boltzmann, Boltzmann]

    data.bulk_recombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip, 
                                                                  bulk_recomb_Auger = false,
                                                                  bulk_recomb_radiative = true,
                                                                  bulk_recomb_SRH = true)

    ## possible choices: ohmic_contact, schottky_contact (outer boundary) and interface_model_none,
    ## interface_model_surface_recombination (inner boundary).
    data.boundary_type[bregionAcceptor]  = ohmic_contact  
    data.boundary_type[bregionJunction1] = interface_model_discont_qF
    data.boundary_type[bregionJunction2] = interface_model_discont_qF                 
    data.boundary_type[bregionDonor]     = ohmic_contact   

    data.isContinuous[iphin]             = false
    data.isContinuous[iphip]             = false
    
    ## possible choices: scharfetter_gummel, scharfetter_gummel_graded, excess_chemical_potential,
    ## excess_chemical_potential_graded, diffusion_enhanced, generalized_sg
    data.flux_approximation              = excess_chemical_potential
    
    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                                       = Params(grid, numberOfCarriers)

    params.temperature                                           = T
    params.UT                                                    = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                                  = -1
    params.chargeNumbers[iphip]                                  =  1

    for ireg in 1:numberOfRegions ## interior region data

        params.dielectricConstant[ireg]                          = ε[ireg]

        ## effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]                      = NC[ireg]
        params.densityOfStates[iphip, ireg]                      = NV[ireg]

        params.bandEdgeEnergy[iphin, ireg]                       = EC[ireg] 
        params.bandEdgeEnergy[iphip, ireg]                       = EV[ireg] 

        params.mobility[iphin, ireg]                             = μn[ireg]
        params.mobility[iphip, ireg]                             = μp[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg]                      = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]             = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]             = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg]          = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg]          = trap_density!(iphip, ireg, data, EI[ireg])
    end

    ## outer boundary region data
    params.bDensityOfStates[iphin, bregionAcceptor]              = Nc_a
    params.bDensityOfStates[iphip, bregionAcceptor]              = Nv_a

    params.bDensityOfStates[iphin, bregionDonor]                 = Nc_d
    params.bDensityOfStates[iphip, bregionDonor]                 = Nv_d

    params.bBandEdgeEnergy[iphin, bregionAcceptor]               = Ec_a
    params.bBandEdgeEnergy[iphip, bregionAcceptor]               = Ev_a

    params.bBandEdgeEnergy[iphin, bregionDonor]                  = Ec_d
    params.bBandEdgeEnergy[iphip, bregionDonor]                  = Ev_d

    # ##############################################################
    d      = 6.28 * 10e-8 * cm # lattice size perovskite
    delta1 = 0.0 * eV
    delta2 = 0.0 * eV

    ## inner boundary region data
    params.bDensityOfStates[iphin_b1, bregionJunction1]          = d * params.densityOfStates[iphin, regionIntrinsic] 
    params.bDensityOfStates[iphip_b1, bregionJunction1]          = d * params.densityOfStates[iphip, regionIntrinsic] 

    params.bBandEdgeEnergy[iphin_b1, bregionJunction1]           = params.bandEdgeEnergy[iphin, regionIntrinsic] + delta1
    params.bBandEdgeEnergy[iphip_b1, bregionJunction1]           = params.bandEdgeEnergy[iphip, regionIntrinsic] + delta2

    params.bMobility[iphin_b1, bregionJunction1]                 = params.mobility[iphin, regionIntrinsic]
    params.bMobility[iphip_b1, bregionJunction1]                 = params.mobility[iphip, regionIntrinsic]

    ##############################################################
    
    ## interior doping
    params.doping[iphin,  regionDonor]                           = Nd 
    params.doping[iphip,  regionAcceptor]                        = Na                     
    ## boundary doping
    params.bDoping[iphin, bregionDonor]                          = Nd 
    params.bDoping[iphip, bregionAcceptor]                       = Na    

    data.params                                                  = params
    ctsys                                                        = System(grid, data, unknown_storage=unknown_storage)

    ## print all params stored in ctsys.data.params
    if test == false
        show_params(ctsys)
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define outer boundary conditions")
    end
    ################################################################################

    ## set zero voltage ohmic contacts for each charge carrier at all outer boundaries.
    set_ohmic_contact!(ctsys, bregionAcceptor, 0.0)
    set_ohmic_contact!(ctsys, bregionDonor,    0.0)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control                   = VoronoiFVM.NewtonControl()
    control.verbose           = verbose
    control.max_iterations    = 300
    control.tol_absolute      = 1.0e-10
    control.tol_relative      = 1.0e-10
    control.handle_exceptions = true
    control.tol_round         = 1.0e-10
    control.max_round         = 5

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    control.damp_initial      = 0.05
    control.damp_growth       = 1.21 # >= 1
    control.max_round         = 5

    ## initialize solution and starting vectors
    initialGuess              = unknowns(ctsys)
    solution                  = unknowns(ctsys)

    solution                  = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

    initialGuess             .= solution 

    if plotting == true

        vis = GridVisualizer(Plotter = PyPlot, layout=(3,1))

        subgrids = VoronoiFVM.subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
        phin_sol = VoronoiFVM.views(solution, data.chargeCarrierList[iphin], subgrids, ctsys.fvmsys)
        phip_sol = VoronoiFVM.views(solution, data.chargeCarrierList[iphip], subgrids, ctsys.fvmsys)
        psi_sol  = VoronoiFVM.views(solution, data.index_psi, subgrids, ctsys.fvmsys)
    
        for i = 1:length(phin_sol)
            scalarplot!(vis[1, 1], subgrids[i], phin_sol[i], clear = false, color=:green)
            scalarplot!(vis[1, 1], subgrids[i], phip_sol[i], clear = false, color=:red)
            scalarplot!(vis[1, 1], subgrids[i], psi_sol[i],  clear = false, color=:blue)
            if i == 3
                scalarplot!(vis[1, 1], subgrids[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                scalarplot!(vis[1, 1], subgrids[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$",  color=:red)
                scalarplot!(vis[1, 1], subgrids[i], psi_sol[i],  clear = false, label = "\$ \\psi \$",color=:blue)
            end
        end

        sol_ref_EQ = readdlm("data/reference-sol-EQ.dat")
        PyPlot.plot(sol_ref_EQ[:, 1], sol_ref_EQ[:, 2], linestyle="--", color = "black")
        PyPlot.plot(sol_ref_EQ[:, 1], sol_ref_EQ[:, 3], linestyle="--", color = "black")
        PyPlot.plot(sol_ref_EQ[:, 1], sol_ref_EQ[:, 4], linestyle="--", color = "black")
        Plotter.legend(fancybox = true, loc = "best", fontsize=11)
        Plotter.title("Solution in EQ")

    end
    #writedlm("reference-sol-EQ.dat", [coord solution'])

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    # Set calculation type to outOfEquilibrium for starting with respective simulation.
    ctsys.data.calculation_type      = outOfEquilibrium
    control.damp_initial          = 0.5
    control.damp_growth           = 1.21
    control.max_round             = 7
    control.max_iterations    = 300

    maxBias    = voltageAcceptor # bias goes until the given contactVoltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 101)
    IV         = zeros(0)

    for Δu in biasValues

        println("Δu  = ", Δu )

        ## set non equilibrium boundary conditions
        set_ohmic_contact!(ctsys, bregionAcceptor, Δu)

        solve!(solution, initialGuess, ctsys, control = control, tstep = Inf)

        subgrids = VoronoiFVM.subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
        phin_sol = VoronoiFVM.views(solution, data.chargeCarrierList[iphin], subgrids, ctsys.fvmsys)
        phip_sol = VoronoiFVM.views(solution, data.chargeCarrierList[iphip], subgrids, ctsys.fvmsys)
        psi_sol  = VoronoiFVM.views(solution, data.index_psi, subgrids, ctsys.fvmsys)
    
        for i = 1:length(phin_sol)
            scalarplot!(vis[2, 1], subgrids[i], phin_sol[i], clear = false, color=:green)
            scalarplot!(vis[2, 1], subgrids[i], phip_sol[i], clear = false, color=:red)
            scalarplot!(vis[2, 1], subgrids[i], psi_sol[i],  clear = false, color=:blue)
            if i == 3
                scalarplot!(vis[2, 1], subgrids[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                scalarplot!(vis[2, 1], subgrids[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$",  color=:red)
                scalarplot!(vis[2, 1], subgrids[i], psi_sol[i],  clear = false, label = "\$ \\psi \$",color=:blue)
            end
        end

        initialGuess .= solution

        factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)
        tf      = testfunction(factory, [1], [2])
        I       = integrate(ctsys.fvmsys, tf, solution)

        val = 0.0
        for ii = 1:length(I)-1
            val = val + I[ii]
        end

        push!(IV, abs(val) )

    end # bias loop

    # writedlm("reference-sol.dat", [coord solution'])
    # res = [biasValues IV]
    # writedlm("reference-IV.dat", res)

    if test == false
        println("*** done\n")
    end
    
    if plotting == true


        subgrids = VoronoiFVM.subgrids(data.chargeCarrierList[iphin], ctsys.fvmsys)
        phin_sol = VoronoiFVM.views(solution, data.chargeCarrierList[iphin], subgrids, ctsys.fvmsys)
        phip_sol = VoronoiFVM.views(solution, data.chargeCarrierList[iphip], subgrids, ctsys.fvmsys)
        psi_sol  = VoronoiFVM.views(solution, data.index_psi, subgrids, ctsys.fvmsys)
    
        for i = 1:length(phin_sol)
            scalarplot!(vis[2, 1], subgrids[i], phin_sol[i], clear = false, color=:green)
            scalarplot!(vis[2, 1], subgrids[i], phip_sol[i], clear = false, color=:red)
            scalarplot!(vis[2, 1], subgrids[i], psi_sol[i],  clear = false, color=:blue)
            if i == 3
                scalarplot!(vis[2, 1], subgrids[i], phin_sol[i], clear = false, label = "\$ \\varphi_n \$", color=:green)
                scalarplot!(vis[2, 1], subgrids[i], phip_sol[i], clear = false, label = "\$ \\varphi_p \$",  color=:red)
                scalarplot!(vis[2, 1], subgrids[i], psi_sol[i],  clear = false, label = "\$ \\psi \$",color=:blue)
            end
        end

        sol_ref = readdlm("data/reference-sol.dat")
        PyPlot.plot(sol_ref[:, 1], sol_ref[:, 2], linestyle="--", color = "black")
        PyPlot.plot(sol_ref[:, 1], sol_ref[:, 3], linestyle="--", color = "black")
        PyPlot.plot(sol_ref[:, 1], sol_ref[:, 4], linestyle="--", color = "black")
        Plotter.legend(fancybox = true, loc = "best", fontsize=11)
        Plotter.title("Solution with Bias")
        ##########################################################
        scalarplot!(vis[3, 1], biasValues, IV, clear = false, color=:green)
        IV_ref         = readdlm("data/reference-IV.dat")
        PyPlot.plot(IV_ref[:, 1], IV_ref[:, 2], linestyle="--", color = "black")

    end

    testval = VoronoiFVM.norm(ctsys.fvmsys, solution, 2)
    return testval


end #  main

function test()
    #testval = 49.92777921026983
    main(test = true, unknown_storage=:dense) ≈ testval #&& main(test = true, unknown_storage=:sparse) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
