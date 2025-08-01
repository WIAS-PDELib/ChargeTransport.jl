#=
# PSC device on 2D domain (unstructured grid).
([source code](SOURCE_URL))

Simulating a three layer PSC device PCBM | MAPI | Pedot with mobile ions. The simulations are
performed in 2D on an unstructured grid, out of equilibrium and with abrupt interfaces.

=#

ENV["LC_NUMERIC"] = "C" # put this in to work with Triangulate.jl, which is originally written in c++

module PSC_2D_unstructuredGrid

    using ChargeTransport
    using ExtendableGrids
    using GridVisualize

    ## For using this example one additionally needs to add Triangulate.
    ## SimplexGridFactory is a wrapper for using this meshgenerator.
    using SimplexGridFactory
    using Triangulate

    ## problem with linux, when including PyPlot not until the end: "ERROR: LoadError: InitError: could not load library "/home/abdel/.julia/artifacts/8cc532f6a1ace8d1b756fc413f4ab340195ec3c3/lib/libgio-2.0.so"/home/abdel/.julia/artifacts/8cc532f6a1ace8d1b756fc413f4ab340195ec3c3/lib/libgobject-2.0.so.0: undefined symbol: g_uri_ref"
    ## It seems that this problem is common: https://discourse.julialang.org/t/could-not-load-library-librsvg-very-strange-error/21276
    using PyPlot

    function main(
            Plotter = PyPlot, ; plotting = false, verbose = false, test = false,
            parameter_set = Params_PSC_PCBM_MAPI_Pedot, # choose the parameter set
        )

        PyPlot.close("all")

        ################################################################################
        if test == false
            println("Define physical parameters and model")
        end
        ################################################################################

        @local_unitfactors μm cm s ns V K ps Hz W

        # parameter with variation
        p = parameter_set(numberOfBoundaryRegions = 5)
        bregionNoFlux = 5

        height = 5.0e-6 * cm

        ## contact voltage
        voltageAcceptor = 1.2 * V

        ## primary data for I-V scan protocol
        scanrate = 0.4 * V / s
        number_tsteps = 31
        endVoltage = voltageAcceptor # bias goes until the given voltage at acceptor boundary

        ## with fixed timestep sizes we can calculate the times a priori
        tend = endVoltage / scanrate
        tvalues = range(0, stop = tend, length = number_tsteps)

        if test == false
            println("*** done\n")
        end
        ################################################################################
        if test == false
            println("Set up grid and regions")
        end
        ################################################################################

        b = SimplexGridBuilder(Generator = Triangulate)

        ## specify boundary nodes
        length_0 = point!(b, 0.0, 0.0)
        length_n = point!(b, p.h_ndoping, 0.0)
        length_ni = point!(b, p.h_ndoping + p.h_intrinsic, 0.0)
        length_nip = point!(b, p.h_total, 0.0)

        height_0 = point!(b, 0.0, height)
        height_n = point!(b, p.h_ndoping, height)
        height_ni = point!(b, p.h_ndoping + p.h_intrinsic, height)
        height_nip = point!(b, p.h_total, height)

        ## specify boundary regions
        ## metal interface
        facetregion!(b, p.bregionDonor)
        facet!(b, length_0, height_0)
        facetregion!(b, p.bregionAcceptor)
        facet!(b, length_nip, height_nip)

        ## no flux
        facetregion!(b, bregionNoFlux)
        facet!(b, length_0, length_nip)
        facetregion!(b, bregionNoFlux)
        facet!(b, height_0, height_nip)

        ## inner interface
        facetregion!(b, p.bregionJ1)
        facet!(b, length_n, height_n)
        facetregion!(b, p.bregionJ2)
        facet!(b, length_ni, height_ni)

        ## cell regions
        cellregion!(b, p.regionDonor)
        regionpoint!(b, p.h_ndoping / 2, height / 2)
        cellregion!(b, p.regionIntrinsic)
        regionpoint!(b, p.h_ndoping + p.h_intrinsic / 2, height / 2)
        cellregion!(b, p.regionAcceptor)
        regionpoint!(b, p.h_ndoping + p.h_intrinsic + p.h_pdoping / 2, height / 2)

        options!(b, maxvolume = 1.0e-16)

        grid = simplexgrid(b)

        if plotting
            GridVisualize.gridplot(grid, Plotter = Plotter, resolution = (600, 400), linewidth = 0.5, legend = :lt)
            Plotter.title("Grid")
        end

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Define System and fill in information about model")
        end
        ################################################################################

        ## Initialize Data instance and fill in data
        data = Data(grid, p.numberOfCarriers)

        ## Possible choices: Stationary, Transient
        data.modelType = Transient

        ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA, FermiDiracMinusOne, Blakemore
        data.F = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]

        data.bulkRecombination = set_bulk_recombination(;
            iphin = p.iphin, iphip = p.iphip,
            bulk_recomb_Auger = false,
            bulk_recomb_radiative = true,
            bulk_recomb_SRH = true
        )

        ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
        ## InterfaceRecombination (inner boundary).
        data.boundaryType[p.bregionAcceptor] = OhmicContact
        data.boundaryType[p.bregionDonor] = OhmicContact

        ## Present ionic vacancies in perovskite layer
        enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionIntrinsic])

        ## Choose flux discretization scheme: ScharfetterGummel, ScharfetterGummelGraded,
        ## ExcessChemicalPotential, ExcessChemicalPotentialGraded, DiffusionEnhanced, GeneralizedSG
        data.fluxApproximation .= ExcessChemicalPotential

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Define Params and fill in physical parameters")
        end
        ################################################################################

        data.params = Params(p)
        ctsys = System(grid, data, unknown_storage = :sparse)

        ## print data
        if test == false
            show_params(ctsys)
            println("*** done\n")
        end
        ################################################################################
        if test == false
            println("Define control parameters for Solver")
        end
        ################################################################################

        control = SolverControl()
        control.verbose = verbose
        control.maxiters = 300
        control.max_round = 5

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Compute solution in thermodynamic equilibrium for Boltzmann")
        end
        ################################################################################

        solution = equilibrium_solve!(ctsys, control = control)
        inival = solution

        if plotting # currently, plotting the solution was only tested with PyPlot.
            ipsi = data.index_psi
            X = grid[Coordinates][1, :]
            Y = grid[Coordinates][2, :]

            Plotter.figure()
            Plotter.surf(X[:], Y[:], solution[ipsi, :])
            Plotter.title("Electrostatic potential \$ \\psi \$ in Equilibrium")
            Plotter.xlabel("length [m]")
            Plotter.ylabel("height [m]")
            Plotter.zlabel("potential [V]")
            Plotter.tight_layout()
            ################
            Plotter.figure()
            Plotter.surf(X[:], Y[:], solution[iphin, :])
            Plotter.title("quasi Fermi potential \$ \\varphi_n \$ in Equilibrium")
            Plotter.xlabel("length [m]")
            Plotter.ylabel("height [m]")
            Plotter.zlabel("potential [V]")
            Plotter.tight_layout()
        end

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("I-V Measurement Loop")
        end
        ################################################################################

        ## for saving I-V data
        IV = zeros(0) # for IV values
        biasValues = zeros(0) # for bias values

        for istep in 2:number_tsteps

            t = tvalues[istep]       # Actual time
            Δu = t * scanrate         # Applied voltage
            Δt = t - tvalues[istep - 1] # Time step size

            ## Apply new voltage; set non equilibrium boundary conditions
            set_contact!(ctsys, p.bregionAcceptor, Δu = Δu)

            if test == false
                println("time value: t = $(t) s")
            end

            solution = solve(ctsys, inival = inival, control = control, tstep = Δt)

            ## get I-V data
            current = get_current_val(ctsys, solution, inival, Δt)

            push!(IV, current)
            push!(biasValues, Δu)

            inival = solution
        end # time loop

        if test == false
            println("*** done\n")
        end

        if plotting
            Plotter.figure()
            Plotter.surf(X[:], Y[:], solution[p.ipsi, :])
            Plotter.title("Electrostatic potential \$ \\psi \$ at end time")
            Plotter.xlabel("length [m]")
            Plotter.ylabel("height [m]")
            Plotter.zlabel("potential [V]")
            ## ################
            Plotter.figure()
            Plotter.surf(X[:], Y[:], solution[p.iphin, :])
            Plotter.title("quasi Fermi potential \$ \\varphi_n \$ at end time")
            Plotter.xlabel("length [m]")
            Plotter.ylabel("height [m]")
            Plotter.zlabel("potential [V]")
            ## ################
            Plotter.figure()
            Plotter.plot(biasValues, IV .* (cm)^2 / height, label = "", linewidth = 3, marker = "o")
            PyPlot.grid()
            Plotter.ylabel("total current [A]") #
            Plotter.xlabel("Applied Voltage [V]")
        end

        testval = sum(filter(!isnan, solution)) / length(solution) # when using sparse storage, we get NaN values in solution
        return testval

    end #  main

    function test()
        testval = -0.5694033507574118
        return main(test = true) ≈ testval
    end

    if test == false
        println("This message should show when this module is successfully recompiled.")
    end

end # module
