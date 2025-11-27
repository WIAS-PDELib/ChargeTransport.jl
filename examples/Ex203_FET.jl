#=
# Field effect transistor.
([source code](@__SOURCE_URL__))

We consider an n-channel Metal-Oxide-Semiconductor (MOS) field effect transistor.
The material is silicon.
=#

module Ex203_FET

using ChargeTransport
using ExtendableGrids
#using LessUnitful: @unitfactors # instead of @local_unitfactors

#using PlutoVista
#using PlutoUI
using GLMakie
using GridVisualize # gridplot function from ChargeTransport not working in 2D
#using PyPlot




function main( ; plotting = true, Plotter = GLMakie, test = false)

    if plotting
        Plotter.close("all")
    end

     # unit factors
    @local_unitfactors μm cm s ns V K

    # constants
    constants = ChargeTransport.constants
    (; q, k_B, ε_0) = constants

    eV = q * V

    ########## charge carriers ##########
    iphin = 1                            # quasi Fermi potential for electrons
    iphip = 2                            # quasi Fermi potential for holes
    numberOfCarriers = 2    

	########## device geometry ##########
	# region numbers
	region_gate = 1
	region_drain = 2
	region_source = 3
	region_bulk = 4
	
	# boundary region numbers
	bregion_gate = 1
	bregion_drain = 2
	bregion_source = 3
	bregion_bulk = 4
	bregion_neutral = 5

	# depth of the device  [cm]
	zaus=15e-4 * cm	 
		
	########## physical values of Si at room temperature ##########
	# https://de.wikipedia.org/wiki/Silicium
	# https://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html#Masses
	Ec = 1.107 * eV                      # conduction band-edge energy
    Ev = 0.0 * eV                        # valence band-edge energy, referenz
	
    Nc = 3.2e19 / (cm^3)               # conduction band density of states
    Nv = 1.8e19 / (cm^3)               # valence band densitiy of states
    
	# https://www.ioffe.ru/SVA/NSM/Semicond/Si/electric.html#Hall
	# Todo: probieren
    mun = 1000.0 * (cm^2) / (V * s)     # electron mobility
    mup = 350.0 * (cm^2) / (V * s)      # hole mobility
	# https://en.wikipedia.org/wiki/Relative_permittivity
    εr = 11.68 * 1.0			        # relative dielectric permittivity of Si
    T = 300.0 * K                       # room temperature

	# Recombination parameters
    Auger = 1.0e-29 * cm^6 / s
    SRH_TrapDensity = 1.0e10 / cm^3
    SRH_LifeTime = 1.0 * ns
    Radiative = 1.0e-10 * cm^3 / s

	# doping (Vereinfacht, nur eine Dotierung pro Region)
	Na_gate = 1e16 / cm^3
	Nd_drain  = 1e19 / cm^3
	Nd_source   = 1e19 / cm^3
	Na_bulk    = 1e15 / cm^3

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # Gridpoints x-direction
	x0 = -1.5 * μm              # beginning left contact (drain)
	x1 = -1.0 * μm              # end left contact (drain)
    x2 = -0.9 * μm              # beginning gate contact
    x3 = 0.9 * μm               # end gate contact
    x4 = 1.0 * μm               # beginning right contact (source)
    x5 = 1.5 * μm               # end right contact (source)     

    # Gridpoints y-direction
    y0 = -2.4 * μm              # bottom of the device (bulk)
    y1 = -0.2 * μm              # beginning n-channel (gate region)
    y2 = 0.0 * μm               # top of the device

    # Refinement x-direction
	X1 = geomspace(x0, x1, 2.0e-7, 5.0e-8)
    X2 = collect(x1:0.2*μm:x4)
    X3 = geomspace(x4, x5, 5.0e-8, 2.0e-7)
	X_temp = glue(X1, X2)
    X = glue(X_temp, X3)
	
	# Refinement y-direction
	Y1 = collect(-2.4:0.2:-1.2) .* μm
	Y2 = collect(-1.2:0.1:-0.1) .* μm
	Y3 = geomspace(-0.1* μm, 0.0* μm, 4.0e-8, 4.0e-8)
	Y12 = glue(Y1,Y2)
	Y = glue(Y12,Y3)
	
	grid = simplexgrid(X, Y) # hier sind jetzt die nx=15 drinn

	# cell regions
	cellmask!(grid, [x1, y1], [x4, y2], region_gate)
	cellmask!(grid, [x0, y1], [x1, y2], region_drain)
	cellmask!(grid, [x4, y1], [x5, y2], region_source)
	cellmask!(grid, [x0, y0], [x5, y1], region_bulk)
	
	# boundary regions 
    bfacemask!(grid, [x0, y0], [x5, y2], bregion_neutral)
	bfacemask!(grid, [x2, y2], [x3, y2], bregion_gate)
	bfacemask!(grid, [x0, y2], [x1, y2], bregion_drain)
	bfacemask!(grid, [x4, y2], [x5, y2], bregion_source)
	bfacemask!(grid, [x0, y0], [x5, y0], bregion_bulk)
    

    if plotting
        GridVisualize.gridplot(grid, Plotter = Plotter, resolution = (1000, 500), title = "Grid")
    end

    params = Params(grid[NumCellRegions], grid[NumBFaceRegions], numberOfCarriers)

	params.temperature = T
	params.chargeNumbers[iphin] = -1
	params.chargeNumbers[iphip] = 1

	for ireg in 1:grid[NumCellRegions] # region data

		params.dielectricConstant[ireg] = εr * ε_0

		# effective DOS, band-edge energy and mobilities
		params.densityOfStates[iphin, ireg] = Nc
		params.densityOfStates[iphip, ireg] = Nv
		params.bandEdgeEnergy[iphin, ireg] = Ec
		params.bandEdgeEnergy[iphip, ireg] = Ev
		params.mobility[iphin, ireg] = mun
		params.mobility[iphip, ireg] = mup

		# recombination parameters, wird das immer übergeben? auch wenn Recombination = false?
		params.recombinationRadiative[ireg] = Radiative
		params.recombinationSRHLifetime[iphin, ireg] = SRH_LifeTime
		params.recombinationSRHLifetime[iphip, ireg] = SRH_LifeTime
		params.recombinationSRHTrapDensity[iphin, ireg] = SRH_TrapDensity
		params.recombinationSRHTrapDensity[iphip, ireg] = SRH_TrapDensity
		params.recombinationAuger[iphin, ireg] = Auger
		params.recombinationAuger[iphip, ireg] = Auger
	end

	# doping                                                        #     n    p
	params.doping[iphin, region_drain] = Nd_drain 					#   [Nd  0.0;
	params.doping[iphin, region_source] = Nd_source                #    Nd  0.0;
	params.doping[iphip, region_gate] = Na_gate                    #    0.0  Na;
	params.doping[iphip, region_bulk] = Na_bulk                    #    0.0  na]

    # Initialize Data instance
	data = Data(grid, numberOfCarriers)

    data.modelType = Stationary

	# statistics
	# Manual Seite 1-3: Boltzmann or Fermi-Dirac, glaube der Default ist Bolzmann (siehe Seite 79)
	#data.F = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]
	data.F .= Boltzmann
	
    # recombination
	# avalanche generation switched of, therefore I switched of everything
	# sind die anderen per Default an oder aus??
	# im Output steht teilweise, dass wir Recombination haben?!
    data.bulkRecombination = set_bulk_recombination(;
        iphin = iphin, iphip = iphip,
        bulk_recomb_Auger = false,
        bulk_recomb_radiative = false,
        bulk_recomb_SRH = false
    )

    # boundary model
	# in Tesca 999 corresponds to Neumann conditions

	# in Tesca Gabez, also gate contacts (natural boundary conditions)
	data.boundaryType[bregion_gate] = OhmicContact # eigentlich Robin!!

	# in Tesca Dibez, also Dirichlet boundary conditions
	data.boundaryType[bregion_drain] = OhmicContact
	data.boundaryType[bregion_source] = OhmicContact
	data.boundaryType[bregion_bulk] = OhmicContact
	
    # flux discretization - depends on statistic (Boltzmann -> Scharfetter Gummel)
    data.fluxApproximation .= ScharfetterGummel

	data.params = params

    # Definition ChargeTransport System
    ctsys = System(grid, data, unknown_storage = :sparse)

    control = SolverControl()
	control.verbose = true # or false

	# Das sind die Werte von der Diode, anpassen?!
	control.maxiters = 50
	control.abstol = 1.0e-7
	control.reltol = 1.0e-7
	control.tol_round = 1.0e-7
	control.damp_initial = 0.5
	control.max_round = 3

    #=
    # Solution in equilibrium
    solution_eq = equilibrium_solve!(ctsys, control = control)

    # Todo: das noch überarbeiten!!
    # Todo: Dichteplots hinzufügen
    if plotting
        # psi aus Lösung rausholen
        psi_eq = solution_eq[3, :]

        # Get coordinates aus dem Grid
        Xpsi_eq = grid[Coordinates][1, :]
        Ypsi_eq = grid[Coordinates][2, :]

        Plotter.figure()
        Plotter.surf(Xpsi_eq[:], Ypsi_eq[:], psi_eq[:])

        Plotter.title("Equilibrium Solution")
        Plotter.xlabel("length [m]")
        Plotter.ylabel("width [m]")
        Plotter.zlabel("potential [V]")
        Plotter.tight_layout()
        #Plotter.gcf()

        # scalarplot(grid, psi_eq, Plotter = GLMakie)
    end
    =#

    #= hier ist noch ein convergence error - erstmal auf Gleichgewichtslösung konzentrieren
    biasValues_gate = range(0.0, stop = 5.0, length = 10)
	biasValues_drain = range(0.0, stop = 5.0, length = 20)

	# IV Werte von Drain
	IV = zeros(0)
	
	# Startwert
	inival = copy(solution_eq)
	solution = copy(solution_eq)

	# Bias Loop
	for Δu_gate in biasValues_gate
		
    	println("bias value at gate: Δu = ", Δu_gate, " V")

		# Müssen source und drain überhaupt angegeben werden??
		set_contact!(ctsys, bregion_source, Δu = 0.0)
		set_contact!(ctsys, bregion_bulk, Δu = 0.0)
		set_contact!(ctsys, bregion_gate, Δu = Δu_gate)

		# Startwerte nach hier verschieben??

		for Δu_drain in biasValues_drain

			println("bias value at drain: Δu = ", Δu_drain, " V")

			set_contact!(ctsys, bregion_drain, Δu = Δu_drain)
			
			solution = solve(ctsys; inival = inival, control = control)
			inival .= solution

			## get I-V data
	        current = get_current_val(ctsys, solution)
	        push!(IV, abs.(zaus * current)) # zaus=wide of device 

		end # bias drain
			
	end # bias gate

    =#

end # function main

end # module




