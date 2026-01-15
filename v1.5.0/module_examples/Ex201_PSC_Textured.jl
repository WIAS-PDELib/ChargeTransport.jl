#=
# 2D textured perovskite solar cell.
([source code](@__SOURCE_URL__))

Simulating a three layer textured PSC device with mobile ions.
The simulations are performed in 2D on an unstructured grid.
This is a minimal example related to the manuscript
https://doi.org/10.48550/arXiv.2506.10691.

=#

ENV["LC_NUMERIC"] = "C" # put this in to work with Triangulate.jl, which is originally written in c++

module Ex201_PSC_Textured

    using ChargeTransport
    using ExtendableGrids
    using VoronoiFVM

    ## For using this example one additionally needs to add Triangulate.
    ## SimplexGridFactory is a wrapper for using this meshgenerator.
    using SimplexGridFactory
    using Triangulate

    function generate_grid(; parameter_set, amplitude = 4.0e-7)

        @local_unitfactors nm

        # use the destructuring operator to extract all the necessary parameters
        (;
            h_ETL, h_HTL, regionETL, regionPero, regionHTL, bregionLeft,
            bregionRight, bregionJ1, bregionJ2, heightDev,
        ) = parameter_set()

        # length of perovskite region depends on amplitude of cos
        h_active = 400 * nm - amplitude / 2 # perovskite
        heightLayers = [h_ETL, h_ETL + h_active, h_ETL + h_active + h_HTL]
        h_total = heightLayers[end]

        maxvolume = 2.0e-15
        eps = 2.0e-10

        ## ETL1  C60
        heightETL1Bottom = 0.0 + 3.0e1 * eps

        heightETL1Top = heightLayers[1] - 20.0 * eps
        heightETL1Top2 = heightLayers[1] - 20.0 * eps
        heightETL1Top3 = heightLayers[1] - 20.0 * eps
        heightETL1Top4 = heightLayers[1] - 20.0 * eps
        heightETL1Top5 = heightLayers[1] - 20.0 * eps
        heightETL1Top6 = heightLayers[1] - 20.0 * eps

        ## pero
        heightPeroBottom = heightLayers[1] + 20.0 * eps
        heightPeroBottom2 = heightLayers[1] + 20.0 * eps
        heightPeroBottom3 = heightLayers[1] + 20.0 * eps
        heightPeroBottom4 = heightLayers[1] + 20.0 * eps
        heightPeroBottom5 = heightLayers[1] + 20.0 * eps
        heightPeroBottom6 = heightLayers[1] + 20.0 * eps

        heightPeroTop = heightLayers[2] - 20.0 * eps
        heightPeroTop2 = heightLayers[2] - 20.0 * eps
        heightPeroTop3 = heightLayers[2] - 20.0 * eps
        heightPeroTop4 = heightLayers[2] - 20.0 * eps
        heightPeroTop5 = heightLayers[2] - 20.0 * eps
        heightPeroTop6 = heightLayers[2] - 20.0 * eps

        heightPero1MaxBottom = 60.0 * nm - 5.0 * nm
        heightPero1MaxTop = 60.0 * nm + 5.0 * nm
        heightPero2MaxBottom = 200.0 * nm - 5.0 * nm
        heightPero2MaxTop = 200.0 * nm + 5.0 * nm
        heightPero3MaxBottom = 350.0 * nm - 5.0 * nm
        heightPero3MaxTop = 350.0 * nm + 5.0 * nm

        # HTL
        heightHTLBottom = heightLayers[2] + 20 * eps
        heightHTLBottom2 = heightLayers[2] + 20 * eps
        heightHTLBottom3 = heightLayers[2] + 20 * eps
        heightHTLBottom4 = heightLayers[2] + 20 * eps
        heightHTLBottom5 = heightLayers[2] + 20 * eps
        heightHTLBottom6 = heightLayers[2] + 20 * eps

        heightHTLTop = heightLayers[3] - 3.0e1 * eps

        #############################################
        ampl = 0.5 * amplitude; phase = pi; period = 7.5e-7

        fCos(x, heightLayer) = ampl .* cos.(phase .+ 2 .* pi .* x ./ period) .+ ampl .+ heightLayer
        fPlanar(x, heightLayer) = heightLayer

        XFine = collect(range(0.0, heightDev, length = 140))

        b = SimplexGridBuilder(Generator = Triangulate)
        ###############################################################
        ## first region (SnO2)
        ###############################################################

        ## nodes
        height_0 = point!(b, heightDev, 0.0)
        height_n = point!(b, heightDev, heightLayers[1])

        length_0 = point!(b, 0.0, 0.0)
        length_n = point!(b, 0.0, heightLayers[1])

        ## facets
        facetregion!(b, bregionLeft); facet!(b, height_0, length_0)
        facetregion!(b, bregionJ1); facet!(b, height_n, length_n)
        facetregion!(b, 0); facet!(b, length_0, length_n)
        facet!(b, height_0, height_n)

        ## refinement
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fPlanar(XFine[ix - 1], heightETL1Bottom))
            B = point!(b, XFine[ix], fPlanar(XFine[ix], heightETL1Bottom))
            facetregion!(b, 0); facet!(b, A, B)
        end

        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fPlanar(XFine[ix - 1], heightETL1Top2))
            B = point!(b, XFine[ix], fPlanar(XFine[ix], heightETL1Top2))
            facetregion!(b, 0); facet!(b, A, B)
        end

        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fPlanar(XFine[ix - 1], heightETL1Top3))
            B = point!(b, XFine[ix], fPlanar(XFine[ix], heightETL1Top3))
            facetregion!(b, 0); facet!(b, A, B)
        end

        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fPlanar(XFine[ix - 1], heightETL1Top4))
            B = point!(b, XFine[ix], fPlanar(XFine[ix], heightETL1Top4))
            facetregion!(b, 0); facet!(b, A, B)
        end

        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fPlanar(XFine[ix - 1], heightETL1Top5))
            B = point!(b, XFine[ix], fPlanar(XFine[ix], heightETL1Top5))
            facetregion!(b, 0); facet!(b, A, B)
        end

        ## regions
        cellregion!(b, regionETL)
        regionpoint!(b, 0.0, 0.6 * heightLayers[1])
        regionpoint!(b, 0.0, 0.99999 * heightETL1Bottom)

        regionpoint!(b, 0.0, 1.00001 * heightETL1Top)
        regionpoint!(b, 0.0, 1.00001 * heightETL1Top2)
        regionpoint!(b, 0.0, 1.00001 * heightETL1Top3)
        regionpoint!(b, 0.0, 1.00001 * heightETL1Top4)
        regionpoint!(b, 0.0, 1.00001 * heightETL1Top5)
        regionpoint!(b, 0.0, 1.00001 * heightETL1Top6)


        ###############################################################
        ## second region (perovskite)
        ###############################################################

        ## nodes
        length_ni = point!(b, 0.0, heightLayers[2])
        height_ni = point!(b, heightDev, heightLayers[2])

        ## facets
        facet!(b, length_n, length_ni)
        facet!(b, height_n, height_ni)

        ## sinusoidal boundary
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightLayers[2]))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightLayers[2]))
            facetregion!(b, bregionJ2); facet!(b, A, B)
        end

        ## refinement

        ## bottom
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fPlanar(XFine[ix - 1], heightPeroBottom))
            B = point!(b, XFine[ix], fPlanar(XFine[ix], heightPeroBottom))
            facetregion!(b, 0); facet!(b, A, B)
        end

        ## bottom
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fPlanar(XFine[ix - 1], heightPeroBottom3))
            B = point!(b, XFine[ix], fPlanar(XFine[ix], heightPeroBottom3))
            facet!(b, A, B)
        end

        ## bottom
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fPlanar(XFine[ix - 1], heightPeroBottom4))
            B = point!(b, XFine[ix], fPlanar(XFine[ix], heightPeroBottom4))
            facet!(b, A, B)
        end

        ## bottom
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fPlanar(XFine[ix - 1], heightPeroBottom5))
            B = point!(b, XFine[ix], fPlanar(XFine[ix], heightPeroBottom5))
            facet!(b, A, B)
        end

        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightPeroTop))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightPeroTop))
            facet!(b, A, B)
        end

        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightPeroTop2))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightPeroTop2))
            facet!(b, A, B)
        end

        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightPeroTop3))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightPeroTop3))
            facet!(b, A, B)
        end

        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightPeroTop4))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightPeroTop4))
            facet!(b, A, B)
        end


        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightPeroTop5))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightPeroTop5))
            facet!(b, A, B)
        end


        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightPeroTop6))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightPeroTop6))
            facet!(b, A, B)
        end

        ## regions
        cellregion!(b, regionPero)
        regionpoint!(b, 0.0, 0.9999 * heightPeroBottom)
        regionpoint!(b, 0.0, 0.9999 * heightPeroBottom2)
        regionpoint!(b, 0.0, 0.9999 * heightPeroBottom3)
        regionpoint!(b, 0.0, 0.9999 * heightPeroBottom4)
        regionpoint!(b, 0.0, 0.9999 * heightPeroBottom5)
        regionpoint!(b, 0.0, 0.9999 * heightPeroBottom6)
        regionpoint!(b, 0.0, 1.5 * heightLayers[1])
        regionpoint!(b, 0.0, 1.0001 * heightPeroTop)
        regionpoint!(b, 0.0, 1.0001 * heightPeroTop)
        regionpoint!(b, 0.0, 1.0001 * heightPeroTop2)
        regionpoint!(b, 0.0, 1.0001 * heightPeroTop3)
        regionpoint!(b, 0.0, 1.0001 * heightPeroTop4)
        regionpoint!(b, 0.0, 1.0001 * heightPeroTop5)
        regionpoint!(b, 0.0, 1.0001 * heightPeroTop6)

        regionpoint!(b, 0.0, 0.9999 * heightPero1MaxBottom)
        regionpoint!(b, 0.0, 0.9999 * heightPero1MaxTop)
        regionpoint!(b, 0.0, 0.9999 * heightPero2MaxBottom)
        regionpoint!(b, 0.0, 0.9999 * heightPero2MaxTop)
        regionpoint!(b, 0.0, 0.9999 * heightPero3MaxBottom)
        regionpoint!(b, 0.0, 1.0001 * heightPero3MaxBottom)
        regionpoint!(b, 0.0, 1.0001 * heightPero3MaxTop)

        ###############################################################
        ## third region (PTAA)
        ###############################################################

        ## nodes
        length_nip = point!(b, 0.0, heightLayers[3])
        height_nip = point!(b, heightDev, heightLayers[3])

        ## facets
        facet!(b, length_ni, length_nip)
        facet!(b, height_ni, height_nip)

        ## sinusoidal boundary
        for ix in 2:length(XFine)

            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightLayers[3]))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightLayers[3]))

            facetregion!(b, bregionRight); facet!(b, A, B)

        end

        # refinement

        ## bottom
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightHTLBottom))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightHTLBottom))
            facetregion!(b, 0); facet!(b, A, B)
        end


        ## bottom
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightHTLBottom3))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightHTLBottom3))
            facet!(b, A, B)
        end

        ## bottom
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightHTLBottom4))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightHTLBottom4))
            facet!(b, A, B)
        end

        ## bottom
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightHTLBottom5))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightHTLBottom5))
            facet!(b, A, B)
        end

        ## top
        for ix in 2:length(XFine)
            A = point!(b, XFine[ix - 1], fCos(XFine[ix - 1], heightHTLTop))
            B = point!(b, XFine[ix], fCos(XFine[ix], heightHTLTop))
            facet!(b, A, B)
        end


        ## regions
        cellregion!(b, regionHTL)
        regionpoint!(b, 0.0, 0.9999 * heightHTLBottom)
        regionpoint!(b, 0.0, 0.9999 * heightHTLBottom2)
        regionpoint!(b, 0.0, 0.9999 * heightHTLBottom3)
        regionpoint!(b, 0.0, 0.9999 * heightHTLBottom4)
        regionpoint!(b, 0.0, 0.9999 * heightHTLBottom5)
        regionpoint!(b, 0.0, 0.9999 * heightHTLBottom6)

        regionpoint!(b, 0.0, 1.00001 * heightHTLTop)
        regionpoint!(b, 0.0, heightLayers[2] + 0.5 * h_HTL)

        ###############################################################
        ## final
        ###############################################################

        options!(b, maxvolume = maxvolume)

        grid = simplexgrid(b)

        bfacemask!(grid, [0.0, 0.0], [0.0, heightLayers[3]], 0, tol = 1.0e-20)
        bfacemask!(grid, [heightDev, 0.0], [heightDev, heightLayers[3]], 0, tol = 1.0e-20)

        return grid

    end

    # you can also use other Plotters, if you add them to the example file
    # you can set verbose also to true to display some solver information
    function main(;
            Plotter = nothing, # only Plotter = PythonPlot or Plotter = PyPlot are supported in this example
            verbose = false, test = false,
            amplitude = 0.5e-7,
            parameter_set = Params_PSC_C60_TripleCation_PTAA, # choose the parameter set
            vacancyEnergyCalculation = false,                 # assume the vacancy energy level is either given or not
            vETL = 2000 * ufac"cm" / ufac"s",                 # surface reco velocity at ETL
            vHTL = 500 * ufac"cm" / ufac"s",                  # surface reco velocity at HTL
        )

        if Plotter !== nothing && (nameof(Plotter) ∉ [:PyPlot, :PythonPlot])
            error("Plotting in Ex201_PSC_Textured is only possible for Plotter = PythonPlot")
        end

        if Plotter !== nothing
            Plotter.rc("font", family = "sans-serif", size = 14)
            Plotter.rc("mathtext", fontset = "dejavusans")
            Plotter.close("all")
        end

        ################################################################################
        if test == false
            println("Define physical parameters and model")
        end
        ################################################################################

        @local_unitfactors V cm m s W nm

        (; q) = ChargeTransport.constants

        ######### primary data for I-V scan protocol ##########
        endVoltage = 1.2 * V
        tPrecond = 0.001
        scanrate = 1.0e3 * V / s
        tend = endVoltage / scanrate

        ## Define scan protocol function
        function scanProtocolPrecond(t)
            biasVal = 0.0
            if data.calculationType == InEquilibrium
                biasVal = 0.0
            else
                if 0.0 <= t && t <= tPrecond
                    biasVal = endVoltage
                elseif tPrecond <= t  && t <= tPrecond + tend
                    biasVal = endVoltage .- scanrate * (t - tPrecond)
                elseif tend + tPrecond < t  && t <= tPrecond + 2 * tend
                    biasVal = scanrate * (t - endVoltage / scanrate - tPrecond) # tend = endVoltage / scanrate
                else
                    biasVal = 0.0
                end
            end

            # we need this such that during the computation of the equilibrium solution, we are safe that no bias is applied.
            if !data.generationComplete
                biasVal = 0.0
            end
            return biasVal
        end

        # Apply zero voltage on left boundary and scan protocol on right boundary
        contactVoltageFunction = [zeroVoltage, scanProtocolPrecond]

        # parameter
        p = parameter_set()

        ################################################################################
        if test == false
            println("Set up grid and regions")
        end
        ################################################################################

        grid = generate_grid(parameter_set = parameter_set, amplitude = amplitude)

        if Plotter !== nothing
            gridplot(grid; Plotter, resolution = (600, 400), linewidth = 0.5, legend = :rc)
        end

        if !test
            println("*** done\n")
        end

        function BeerLamb(x, y)

            ampl = 0.5 * amplitude
            phase = pi
            period = 7.5e-7

            # length of perovskite region depends on amplitude of cos
            h_active = 400 * nm - amplitude / 2 # perovskite
            heightLayers = [p.h_ETL, p.h_ETL + h_active, p.h_ETL + h_active + p.h_HTL]

            # optical Parameters
            Fph = p.incidentPhotonFlux[p.regionPero]
            ag = p.absorption[p.regionPero]
            inv = p.invertedIllumination
            genPeak = ampl .* cos.(phase .+ 2 .* pi .* x ./ period) .+ ampl .+ heightLayers[2]

            if heightLayers[1] <= y <= ampl .* cos.(phase .+ 2 .* pi .* x ./ period) .+ ampl .+ heightLayers[2]
                G = Fph .* ag .* exp.(- inv .* ag .* (y .- genPeak))
            else
                G = 0
            end

            return G
        end

        generationData = zeros(length(grid[Coordinates][1, :]))

        i = 0
        for  coord in eachcol(grid[Coordinates])
            i = i + 1
            x0 = coord[1]
            y0 = coord[2]

            generationData[i] = BeerLamb(x0, y0)

        end

        ## Plot of photogeneration
        # XX = grid[Coordinates][1, :]
        # YY = grid[Coordinates][2, :]

        # figure()
        # tricontourf(XX ./ nm, YY ./ nm, generationData, levels = 40)

        # colorbar()

        ################################################################################
        if test == false
            println("Define System and fill in information about model")
        end
        ################################################################################

        ## Initialize Data instance and fill in data
        data = Data(grid, p.numberOfCarriers, contactVoltageFunction = contactVoltageFunction, generationData = generationData)

        ## Possible choices: Stationary, Transient
        data.modelType = Transient

        data.bulkRecombination = set_bulk_recombination(;
            iphin = p.iphin, iphip = p.iphip,
            bulk_recomb_Auger = false,
            bulk_recomb_radiative = true,
            bulk_recomb_SRH = true
        )

        ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
        data.boundaryType[p.bregionLeft] = MixedOhmicSchottkyContact
        data.boundaryType[p.bregionRight] = MixedOhmicSchottkyContact
        data.boundaryType[p.bregionJ1] = InterfaceRecombination
        data.boundaryType[p.bregionJ2] = InterfaceRecombination

        ## currently Beer-Lambert not working for textured devices
        data.generationModel = GenerationUserDefined

        ## Present ionic vacancies in perovskite layer
        ## by default the statistics function is set to FermiDiracMinusOne to limit ion depletion
        enable_ionic_carrier!(data, ionicCarrier = p.iphia, regions = [p.regionPero])

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Define Params and fill in physical parameters")
        end
        ################################################################################

        data.params = Params(p)

        if !vacancyEnergyCalculation
            data.params.bandEdgeEnergy[p.iphia, :] .= p.Ea2D
        end

        data.params.recombinationSRHvelocity[p.iphin, p.bregionJ1] = vETL
        data.params.recombinationSRHvelocity[p.iphip, p.bregionJ1] = vETL

        data.params.recombinationSRHvelocity[p.iphin, p.bregionJ2] = vHTL
        data.params.recombinationSRHvelocity[p.iphip, p.bregionJ2] = vHTL

        ctsys = System(grid, data, unknown_storage = :sparse)

        if test == false
            println("*** done\n")
        end
        ################################################################################
        if test == false
            println("Define control parameters for Solver")
        end
        ################################################################################

        control = ChargeTransport.SolverControl()
        control.verbose = verbose
        control.damp_initial = 0.5
        control.damp_growth = 1.61 # >= 1
        control.maxiters = 300
        control.max_round = 5
        control.abstol = 1.0e-7
        control.reltol = 1.0e-7
        control.tol_round = 1.0e-7
        control.Δu_opt = Inf

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Compute solution in thermodynamic equilibrium for Boltzmann")
        end
        ################################################################################

        #data.params.bandEdgeEnergy[3, 2] = -5.267 * eV
        ## Solve for the equilibrium solution in several steps:
        #   1) Poisson problem is solved,
        #   2) Full system is solved in equilibrium with photogeneration turned on,
        #   Since vacancyEnergyCalculation = true, we repeat 1) and 2) until Ea is found (which corresponds to given initial vacancy density p.Ca)
        solEQ = equilibrium_solve!(ctsys, control = control, vacancyEnergyCalculation = vacancyEnergyCalculation)
        inival = solEQ

        if test == false
            println("*** done\n")
        end
        ################################################################################
        if test == false
            println("Loop to increase bias (so that following scan protocol starts near VOC)")
        end
        ################################################################################

        biasLoop = collect(reverse(-range(0.0, endVoltage, length = 21)))

        sol0p0 = ChargeTransport.unknowns(ctsys)
        for istep in 1:length(biasLoop)

            ## turn slowly voltage on
            set_contact!(ctsys, p.bregionRight, Δu = biasLoop[istep])

            ## keep some big time stepping here, so we preserve the vacancy density
            sol0p0 = ChargeTransport.solve(ctsys, inival = inival, control = control, tstep = 1.0e3)
            inival = sol0p0

        end # bias loop

        if !test
            println("*** done\n")
        end

        # Now we have the correct initial condition to start the scan protocol
        ################################################################################
        if test == false
            println("Preconditioning")
        end
        ################################################################################

        control.Δt = 1.0e-4
        control.Δt_min = 1.0e-5
        control.Δt_max = 5.0e-4
        control.Δt_grow = 1.1

        solPrecond = ChargeTransport.solve(ctsys, inival = inival, times = (0.0, tPrecond), control = control)

        if !test
            println("*** done\n")
        end

        ################################################################################
        if !test
            println("Reverse IV Measurement loop")
        end
        ################################################################################

        control.Δt = 4.0e-3 / scanrate
        control.Δt_min = 4.0e-3 / scanrate
        control.Δt_max = 6.0e-2 / scanrate
        control.Δt_grow = 1.15

        solRev = ChargeTransport.solve(ctsys, inival = solPrecond.u[end], times = (tPrecond, tPrecond + tend), control = control)

        if !test
            println("*** done\n")
        end

        ################################################################################
        if !test
            println("Forward IV Measurement loop")
        end
        ################################################################################

        solForw = ChargeTransport.solve(ctsys, inival = solRev.u[end], times = (tPrecond + tend, tPrecond + 2 * tend), control = control)

        if Plotter !== nothing
            subg = subgrid(grid, [p.regionPero])

            # https://github.com/j-fu/GridVisualize.jl/blob/1f2b299a436b7750702ccca282fa14152d80ebf9/src/pyplot.jl#L86
            function tridata(grid::ExtendableGrid)
                coord = grid[Coordinates] ./ nm
                cellnodes = Matrix(grid[CellNodes])
                return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
            end

            nn = get_density(solForw.u[1], p.regionPero, ctsys, p.iphin)
            np = get_density(solForw.u[1], p.regionPero, ctsys, p.iphip)
            na = get_density(solForw.u[1], p.regionPero, ctsys, p.iphia)

            @show minimum(nn), maximum(nn)
            @show minimum(np), maximum(np)
            @show minimum(na), maximum(na)

            figsize = (7.2, 5.6)
            vmin = 5.0e17; vmax = 6.0e20
            vminIon = 1.0e22; vmaxIon = 5.0e23


            Plotter.figure(figsize = figsize)
            Plotter.tripcolor(tridata(subg)..., vcat(nn...), norm = Plotter.matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
            Plotter.xlabel(" \$x\$ [nm]", fontsize = 17)
            Plotter.ylabel(" \$y\$ [nm]", fontsize = 17)
            Plotter.axis([-20, 770, 20, 800])
            Plotter.title("Electron density (beginning forward scan)")
            Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
            Plotter.tight_layout()

            current_figure = Plotter.gcf()
            display(current_figure)

            #################
            Plotter.figure(figsize = figsize)
            Plotter.tripcolor(tridata(subg)..., vcat(np...), norm = Plotter.matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
            Plotter.xlabel(" \$x\$ [nm]", fontsize = 17)
            Plotter.ylabel(" \$y\$ [nm]", fontsize = 17)
            Plotter.axis([-20, 770, 20, 800])
            Plotter.title("Hole density (beginning forward scan)")
            Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
            Plotter.tight_layout()

            current_figure = Plotter.gcf()
            display(current_figure)

            #################
            Plotter.figure(figsize = figsize)
            Plotter.tripcolor(tridata(subg)..., vcat(na...), norm = Plotter.matplotlib.colors.LogNorm(vmin = vminIon, vmax = vmaxIon), shading = "gouraud", rasterized = true)
            Plotter.xlabel(" \$x\$ [nm]", fontsize = 17)
            Plotter.ylabel(" \$y\$ [nm]", fontsize = 17)
            Plotter.axis([-20, 770, 20, 800])
            Plotter.title("Ion density (beginning forward scan)")
            Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
            Plotter.tight_layout()

            current_figure = Plotter.gcf()
            display(current_figure)

            nn = get_density(solForw.u[end], p.regionPero, ctsys, p.iphin)
            np = get_density(solForw.u[end], p.regionPero, ctsys, p.iphip)
            na = get_density(solForw.u[end], p.regionPero, ctsys, p.iphia)

            println(" ")
            @show minimum(nn), maximum(nn)
            @show minimum(np), maximum(np)
            @show minimum(na), maximum(na)

            Plotter.figure(figsize = figsize)
            Plotter.tripcolor(tridata(subg)..., vcat(nn...), norm = Plotter.matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
            Plotter.xlabel(" \$x\$ [nm]", fontsize = 17)
            Plotter.ylabel(" \$y\$ [nm]", fontsize = 17)
            Plotter.axis([-20, 770, 20, 800])
            Plotter.title("Electron density (end forward scan)")
            Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
            Plotter.tight_layout()

            current_figure = Plotter.gcf()
            display(current_figure)

            #################
            Plotter.figure(figsize = figsize)
            Plotter.tripcolor(tridata(subg)..., vcat(np...), norm = Plotter.matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax), shading = "gouraud", rasterized = true)
            Plotter.xlabel(" \$x\$ [nm]", fontsize = 17)
            Plotter.ylabel(" \$y\$ [nm]", fontsize = 17)
            Plotter.axis([-20, 770, 20, 800])
            Plotter.title("Hole density (end forward scan)")
            Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
            Plotter.tight_layout()

            current_figure = Plotter.gcf()
            display(current_figure)

            #################
            Plotter.figure(figsize = figsize)
            Plotter.tripcolor(tridata(subg)..., vcat(na...), norm = Plotter.matplotlib.colors.LogNorm(vmin = vminIon, vmax = vmaxIon), shading = "gouraud", rasterized = true)
            Plotter.xlabel(" \$x\$ [nm]", fontsize = 17)
            Plotter.ylabel(" \$y\$ [nm]", fontsize = 17)
            Plotter.axis([-20, 770, 20, 800])
            Plotter.title("Ion density (end forward scan)")
            Plotter.colorbar(orientation = "vertical", label = " Density [\$\\mathrm{m}^{-3}\$]", extend = "both")
            Plotter.tight_layout()

            current_figure = Plotter.gcf()
            display(current_figure)
        end

        if !test
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Forward IV curve calculation")
        end
        ################################################################################

        IV = zeros(0) # for saving I-V data
        IVn = zeros(0); IVp = zeros(0)
        IVa = zeros(0); IVψ = zeros(0)

        ######################
        ISRHn = zeros(0); IRadn = zeros(0)
        ISRHp = zeros(0); IRadp = zeros(0)
        ISRnL = zeros(0)
        ISRnR = zeros(0)
        IGen = zeros(0)

        tvalues = solForw.t
        number_tsteps = length(tvalues)
        biasValues = data.contactVoltageFunction[2].(tvalues[2:end])

        factory = ChargeTransport.TestFunctionFactory(ctsys)
        tf = ChargeTransport.testfunction(factory, [p.bregionLeft], [p.bregionRight])

        for istep in 2:number_tsteps

            Δt = tvalues[istep] - tvalues[istep - 1] # Time step size
            inival = solForw.u[istep - 1]
            solution = solForw.u[istep]

            IEdge = VoronoiFVM.integrate_∇TxFlux(ctsys.fvmsys, tf, solution)
            IEdgeOld = VoronoiFVM.integrate_∇TxFlux(ctsys.fvmsys, tf, inival)

            push!(IVn, IEdge[p.iphin]); push!(IVp, IEdge[p.iphip])
            push!(IVψ, (IEdge[p.ipsi] - IEdgeOld[p.ipsi]) / Δt)

            push!(IVa, IEdge[p.iphia])
            push!(IV, IVψ[istep - 1] + IVn[istep - 1] + IVp[istep - 1] + IVa[istep - 1])

            IntSRH = ChargeTransport.integrate(ctsys, SRHRecombination!, solution)
            IntRad = ChargeTransport.integrate(ctsys, RadiativeRecombination!, solution)
            IntSR = ChargeTransport.integrate(ctsys, SRRecombination!, solution, boundary = true)
            IntGen = ChargeTransport.integrate(ctsys, Photogeneration!, solution)

            IntSRHnSum = 0.0; IntRadnSum = 0.0
            IntSRHpSum = 0.0; IntRadpSum = 0.0

            for ii in 1:p.numberOfRegions
                IntSRHnSum = IntSRHnSum - IntSRH[p.iphin, ii]
                IntRadnSum = IntRadnSum - IntRad[p.iphin, ii]

                IntSRHpSum = IntSRHpSum + IntSRH[p.iphip, ii]
                IntRadpSum = IntRadpSum + IntRad[p.iphip, ii]
            end

            IntGenp = IntGen[p.iphip, p.regionPero]
            IntSRnL = - IntSR[p.iphin, p.bregionJ1]
            IntSRnR = - IntSR[p.iphin, p.bregionJ2]

            push!(ISRHn, IntSRHnSum); push!(ISRHp, IntSRHpSum)
            push!(IRadn, IntRadnSum); push!(IRadp, IntRadpSum)
            push!(ISRnL, IntSRnL)
            push!(ISRnR, IntSRnR)
            push!(IGen, IntGenp)

        end

        if Plotter !== nothing

            Plotter.figure()
            tEnd = tPrecond + 2 * tend

            tt = collect(range(0.0, tEnd, length = 201))
            T = data.contactVoltageFunction[2]
            Plotter.plot(tt, T.(tt), marker = "o")
            Plotter.xlabel("time [s]")
            Plotter.ylabel("voltage [V]")
            Plotter.tight_layout()

            current_figure = Plotter.gcf()
            display(current_figure)

            Plotter.figure()
            Plotter.plot(biasValues, -IV .* (cm^2) .* 1.0e3 ./ p.heightDev, linewidth = 5, color = "blue", label = "forward")

            Plotter.grid()
            Plotter.legend()
            Plotter.xlabel("bias [V]", fontsize = 17)
            Plotter.ylabel("current density [mAcm\$^{-2} \$]", fontsize = 17)
            Plotter.tick_params(which = "both", labelsize = 18)
            Plotter.tight_layout()

            current_figure = Plotter.gcf()
            display(current_figure)

            #####################################
            Plotter.figure()
            Plotter.semilogy(biasValues, ISRHn .* (cm^2) .* 1.0e3 ./ p.heightDev, linewidth = 5, color = "blue", label = "SRH")
            Plotter.semilogy(biasValues, IRadn .* (cm^2) .* 1.0e3 ./ p.heightDev, linewidth = 5, color = "red", label = "rad")
            Plotter.semilogy(biasValues, IGen .* (cm^2) .* 1.0e3 ./ p.heightDev, linewidth = 5, color = "gold", label = "Gen")
            Plotter.semilogy(biasValues, ISRnL .* (cm^2) .* 1.0e3 ./ p.heightDev, linewidth = 5, color = "black", label = "SR, left")
            Plotter.semilogy(biasValues, ISRnR .* (cm^2) .* 1.0e3 ./ p.heightDev, linewidth = 5, color = "darkgreen", label = "SR, right")

            Plotter.grid()
            Plotter.legend()
            Plotter.xlabel("bias [V]", fontsize = 17)
            Plotter.ylabel("current density [mAcm\$^{-2} \$]", fontsize = 17)
            Plotter.tick_params(which = "both", labelsize = 18)
            Plotter.tight_layout()

            current_figure = Plotter.gcf()
            display(current_figure)
        end

        if test == false
            IV = -IV
            bias = biasValues

            powerDensity = bias .* (IV)           # power density function
            MaxPD, indexPD = findmax(powerDensity)

            open_circuit = compute_open_circuit_voltage(bias, IV)

            IncLightPowerDens = 1000.0 * W / m^2

            fillfactor = 100 * (bias[indexPD] * IV[indexPD]) / (IV[1] * open_circuit)

            efficiency = 100 * bias[indexPD] * (IV[indexPD] / p.heightDev) / (IncLightPowerDens)
            JSC = IV[1] .* (cm)^(2) .* 1.0e3 ./ p.heightDev

            println(" ")
            println("The JSC                  is $(round(JSC, digits = 3)) mAcm^{-2}.")
            println("The fill factor          is $(round(fillfactor, digits = 2)) %.")
            println("The efficiency           is $(round(efficiency, digits = 2)) %.")
            println("The open circuit voltage is $(round(open_circuit, digits = 4)) V.")
            println(" ")
        end
        if test == false
            integral = integrated_density(ctsys, sol = solEQ, icc = p.iphia, ireg = p.regionPero)
            mOmega = data.regionVolumes[p.regionPero]

            println("Calculated average vacancy density is: ", integral / mOmega)
            println(" ")
            vacancyEnergy = data.params.bandEdgeEnergy[p.iphia, p.regionPero] / q
            println("Value for vacancy energy is: ", vacancyEnergy, " eV. Save this value for later use.")
            println(" ")
        end

        testval = sum(filter(!isnan, solForw.u[end])) / length(solForw.u[end]) # when using sparse storage, we get NaN values in solution
        return testval

    end #  main

    function test()
        testval = -0.5987199947920246
        return main(test = true) ≈ testval
    end

end # module
