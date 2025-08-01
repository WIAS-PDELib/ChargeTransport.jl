#=
Simulating a simple laser structure with 5 layers.
The layers are defined by their material properties and thicknesses.
The simulation will solve the charge transport equations across the layers,
taking into account the stimulated recombination in the active region of the laser structure.
=#

module Ex202_Laser_simple
using ChargeTransport
using ExtendableGrids
using GridVisualize
using PyPlot

###########################################################################

numberOfColoumns = Dict(
    "ref1" => [2, 4],
    "ref2" => [4, 8],
    "ref3" => [8, 16],
    "ref4" => [16, 32],
    "ref5" => [32, 64]
)
numberOfRows = Dict(
    "ref1" => [4, 8, 2, 8, 4],
    "ref2" => [8, 16, 4, 16, 8],
    "ref3" => [16, 32, 8, 32, 16],
    "ref4" => [32, 64, 16, 64, 32],
    "ref5" => [64, 128, 32, 128, 64]
)

###################################################################

""" Initializing X and Y coords for the tesca grid"""
function tesca_grid(; refinement = 1, showplot = false, airbox = false)

    @local_unitfactors μm

    ncol = numberOfColoumns["ref$(refinement)"]
    nrow = numberOfRows["ref$(refinement)"]

    widths_columns = [1.0  10.0] * μm
    heights_rows = [1.0  0.5  0.05  0.5  1.0] * μm

    coord_x1 = collect(range(0.0, widths_columns[1], length = ncol[1] + 1))
    coord_x2 = collect(range(widths_columns[1], sum(widths_columns[1:2]), length = ncol[2] + 1))
    X = glue(coord_x1, coord_x2)
    X = glue(reverse(-X), X)

    coord_y1 = collect(range(0.0, heights_rows[1], length = nrow[1] + 1))
    coord_y2 = collect(range(heights_rows[1], sum(heights_rows[1:2]), length = nrow[2] + 1))
    coord_y3 = collect(range(sum(heights_rows[1:2]), sum(heights_rows[1:3]), length = nrow[3] + 1))
    coord_y4 = collect(range(sum(heights_rows[1:3]), sum(heights_rows[1:4]), length = nrow[4] + 1))
    coord_y5 = collect(range(sum(heights_rows[1:4]), sum(heights_rows[1:5]), length = nrow[5] + 1))
    Y = glue(glue(glue(glue(coord_y1, coord_y2), coord_y3), coord_y4), coord_y5)

    grid = simplexgrid(X, Y)

    cellmask!(grid, [X[1], 0.0], [sum(widths_columns), heights_rows[1]], 1)
    cellmask!(grid, [X[1], heights_rows[1]], [sum(widths_columns), sum(heights_rows[1:2])], 2)
    cellmask!(grid, [X[1], sum(heights_rows[1:2])], [sum(widths_columns), sum(heights_rows[1:3])], 3)
    cellmask!(grid, [X[1], sum(heights_rows[1:3])], [sum(widths_columns), sum(heights_rows[1:4])], 4)
    cellmask!(grid, [X[1], sum(heights_rows[1:4])], [sum(widths_columns), sum(heights_rows[1:5])], 5)

    # AIR
    cellmask!(grid, [widths_columns[1], sum(heights_rows[1:4])], [sum(widths_columns), sum(heights_rows[1:5])], 6)
    cellmask!(grid, [X[1], sum(heights_rows[1:4])], [-widths_columns[1], sum(heights_rows[1:5])], 6)

    bregionDonor1 = 1    # bottom boundary
    bregionAcceptor2 = 2    # top boundary
    bregionNoFlux = 3
    bregionAirBox = 4
    bfacemask!(grid, [-widths_columns[1], sum(heights_rows)], [widths_columns[1], sum(heights_rows)], bregionAcceptor2)
    bfacemask!(grid, [X[1], 0.0], [X[end], 0.0], bregionDonor1)

    bfacemask!(grid, [X[1], 0.0], [X[1], sum(heights_rows[1:4])], bregionNoFlux)
    bfacemask!(grid, [X[1], sum(heights_rows[1:4])], [-widths_columns[1], sum(heights_rows[1:4])], bregionNoFlux)
    bfacemask!(grid, [-widths_columns[1], sum(heights_rows[1:4])], [-widths_columns[1], sum(heights_rows)], bregionNoFlux)

    bfacemask!(grid, [X[end], 0.0], [X[end], sum(heights_rows[1:4])], bregionNoFlux)
    bfacemask!(grid, [widths_columns[1], sum(heights_rows[1:4])], [X[end], sum(heights_rows[1:4])], bregionNoFlux)
    bfacemask!(grid, [widths_columns[1], sum(heights_rows[1:4])], [widths_columns[1], sum(heights_rows)], bregionNoFlux)

    bfacemask!(grid, [X[1], sum(heights_rows[1:4])], [X[1], Y[end]], bregionAirBox)
    bfacemask!(grid, [X[1], Y[end]], [-widths_columns[1], Y[end]], bregionAirBox)
    bfacemask!(grid, [widths_columns[1], Y[end]], [X[end], Y[end]], bregionAirBox)
    bfacemask!(grid, [X[end], sum(heights_rows[1:4])], [X[end], Y[end]], bregionAirBox)


    if airbox == false
        grid = subgrid(grid, [1, 2, 3, 4, 5])
    end

    if showplot == true
        GridVisualize.gridplot(
            grid, Plotter = PyPlot, linewidth = 1, fontsize = 35, size = (1200, 900),
            legend = :best, show = true, aspect = 4, colorbar = false, title = "Device Geometry, values in [m]", xlabel = "x-coordinates", ylabel = "y-coordinates"
        )
    end

    return grid
end

function main(;
        refinement = 1, plotting = false, Plotter = PyPlot, verbose = false, test = false,
        unknown_storage = :sparse, numberOfEigenvalues = 1,
        parameter_set = Params_Laser_simple
    ) # choose the parameter set

    # parameter
    p = parameter_set()

    ################################################################################
    if test == false
        println("Set up grid.")
    end

    grid = tesca_grid(refinement = refinement, showplot = plotting, airbox = false)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in data
    data = Data(grid, p.numberOfCarriers, numberOfEigenvalues = numberOfEigenvalues)

    ## Possible choices: Stationary, Transient
    data.modelType = Stationary

    ## Possible choices: Boltzmann, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA,
    ## FermiDiracMinusOne, Blakemore.
    ## Can be a vector with different statistics (for n and p).
    data.F .= FermiDiracOneHalfTeSCA

    data.bulkRecombination = set_bulk_recombination(;
        iphin = p.iphin, iphip = p.iphip,
        bulk_recomb_Auger = true,
        bulk_recomb_radiative = true,
        bulk_recomb_SRH = true
    )

    ## Possible choices: OhmicContact, SchottkyContact (outer boundary) and InterfaceNone,
    ## InterfaceRecombination (inner boundary).
    data.boundaryType[p.bregionAcceptor2] = OhmicContact    # top boundary Dirichlet condition
    data.boundaryType[p.bregionDonor1] = OhmicContact    # bottom boundary Dirichlet condition
    #                                                     # rest is set to Neumann by default

    data.fluxApproximation .= ExcessChemicalPotential
    #data.fluxApproximation           .= ScharfetterGummel           # if F=Boltzmann

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    """ Data from Params_Laser_simple.jl: temperature T, band edge energies EC, EV, effective density of states NC, NV
        mobilities μn, μp, dielectricConstant εs, radiative recombination r0, life times τn, τp,
        Auger recombination coefficients Auger_Cn, Auger_Cp
        doping doping (or vcat(Nd,Na) = doping).
    """
    paramsoptical = ParamsOptical(grid, p.numberOfCarriers, numberOfEigenvalues)
    paramsoptical.laserWavelength = p.λ

    paramsoptical.absorption_0[:] = p.α0
    paramsoptical.gain_0[:] = p.gain0
    paramsoptical.refractiveIndex_0[:] = p.nTilde
    paramsoptical.refractiveIndex_d[:] = p.nTilde_d
    paramsoptical.refractiveIndex_γ[:] = p.γn
    paramsoptical.absorptionFreeCarriers[p.iphin, :] = p.fcnalf
    paramsoptical.absorptionFreeCarriers[p.iphip, :] = p.fcpalf

    paramsoptical.eigenvalues .= 1 + 1 * im   # dummy value for initializing

    data.params = Params(p)
    data.paramsoptical = paramsoptical

    ctsys = System(grid, data, unknown_storage = unknown_storage)


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
    control.abstol = 1.0e-7
    control.reltol = 1.0e-7
    control.tol_round = 1.0e-7
    control.max_round = 3
    control.damp_initial = 0.8   # < 1
    #control.damp_growth  = 1.21 # >= 1

    if test == false
        println("*** done\n")

    end


    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium for Boltzmann")
    end
    ################################################################################

    ## calculate equilibrium solution and set as initial guess

    psi0Vector = electroNeutralSolution(ctsys)

    inival = unknowns(ctsys)
    inival[1, :] = inival[2, :] .= 0.0
    inival[3, :] = psi0Vector


    solution = equilibrium_solve!(ctsys, inival = inival, control = control, nonlinear_steps = 20.0)
    inival = solution


    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Bias loop")
    end
    ################################################################################

    maxBias = p.U[end]  # = 1.81 = topVoltageAcceptor2 # bias goes until the given voltage at acceptor boundary
    biasValues = range(0, stop = maxBias, length = 40)
    IV = zeros(0)

    for Δu in biasValues

        if test == false
            println("bias value: Δu = ", Δu, " V")
        end

        ## set non equilibrium boundary conditions
        set_contact!(ctsys, p.bregionAcceptor2, Δu = Δu)

        solution = solve(ctsys; inival = inival, control = control)
        inival .= solution

    end # bias loop

    if test == false
        println("*** done\n")
    end

    ################################################################################
    ctsys.data.paramsoptical.oldSolution = solution
    currentSolution = solution
    inival = solution


    # eigenvalue and eigenvector of the Helmholtz eigenvalue problem calculated explicitly
    λ1 = [-5.6192359487570194e14 + 4.6950340673509827e11im]
    v1 = [
        -6.44693721037541e-8 + 2.1403391419628273e-8im,
        -0.00016112043789318325 + 5.3345497772375895e-5im,
        -0.04065180212557369 + 0.013461955814899604im,
        -1.625774595841853e-5 + 5.383258567547116e-6im,
        -0.00015713796062169245 + 6.309731311560341e-5im,
        -0.03964667896215338 + 0.015922689828672086im,
        -2.3324630277054376e-5 - 3.9545140054933005e-5im,
        -0.005884084239403538 - 0.00997848428365716im,
        -0.0001809102600978133 + 0.00014171937193898214im,
        -0.04542888332915838 + 0.03459259677971059im,
        0.0005919624102358168 + 0.0007767767309721376im,
        0.14909867925795559 + 0.1961272646164019im,
        0.003509519704937677 + 0.00023236831835518715im,
        0.8850503810764251 + 0.059473291222819455im,
        0.009849688910224585 - 0.00747820637205809im,
        2.485408855146658 - 1.8845167715203734im,
        0.013192949613002478 - 0.03461148153029405im,
        3.3355658648998205 - 8.74169283250292im,
        -0.3699347591173653 - 0.4496821788206425im,
        -93.32201642056931 - 113.46779411302093im,
        -5.037570813855826 + 2.206072044371644im,
        -1271.002622734984 + 556.6956009204584im,
        -5.26426875045551 + 1.532547420703787im,
        -1328.214932850015 + 386.7518276272265im,
        -0.002108593408263752 + 0.0006138388168361884im,
        -0.5311809720795759 + 0.15465513471318287im,
        -0.1873995001277716 + 0.06332832331090844im,
        -7.495085655365893e-5 + 2.532736765454705e-5im,
        -0.18260712076737765 + 0.07481510962027382im,
        -0.026674803247518492 - 0.04652767633233148im,
        -0.22526294915934103 + 0.14803180699535837im,
        0.5575095767350965 + 0.9742390606680793im,
        3.8649440306697325 + 0.7004811462666269im,
        11.588951763025864 - 7.529272956786549im,
        16.804355803465423 - 39.213894439889735im,
        -422.59189173993883 - 528.0537657454727im,
        -5853.433344563868 + 2610.8857632468134im,
        -6124.670145467389 + 1823.000173170905im,
        -2.4492160883354672 + 0.7289325628684069im,
        -0.8239433305568086 + 0.2845918465840918im,
        -0.00032982678049294384 + 0.00011392326428668871im,
        -0.8021001518012466 + 0.3357933076136317im,
        -0.11496413881386891 - 0.20705060118763416im,
        -1.0580169915962223 + 0.589171602111913im,
        1.8196456330246056 + 4.557107300388768im,
        15.860419842117603 + 5.049668769007414im,
        51.386271699647146 - 27.318805261236268im,
        80.81410006786088 - 167.13357161525897im,
        -1821.4316913399055 - 2346.622366949512im,
        -25708.04375973995 + 11695.399375747844im,
        -26937.09783301187 + 8209.613614522255im,
        -10.762756223235524 + 3.2797040120768766im,
        -3.6138102804798353 + 1.2755859110128367im,
        -0.0014623741175362572 + 0.0005163135313304949im,
        -3.5145630376007198 + 1.5032734988162935im,
        -0.4933285127135774 - 0.9187044435324302im,
        -4.894062401090896 + 2.282956861477433im,
        5.177209414344344 + 20.86473547393878im,
        64.07521579720996 + 30.35275983503732im,
        226.20008398762747 - 93.1138049291388im,
        385.34412265699206 - 709.4475177836089im,
        -7827.623921542349 - 10406.96864810031im,
        -112630.72113222913 + 52256.976582767245im,
        -118184.61437710556 + 36871.82426926498im,
        -46.71738226489013 + 14.569128201542856im,
        -5.502085509328842 + 1.9615977177667154im,
        -0.0022038847386733522 + 0.0007857132316632017im,
        -5.348514404225735 + 2.3105034138414866im,
        -0.743025608472823 - 1.405880395918349im,
        -7.9487316574209474 + 3.4142973339001794im,
        5.909391447819046 + 32.14850166767881im,
        93.34651775903117 + 51.58415410110485im,
        344.1331782080538 - 122.01201710647821im,
        581.404392913808 - 1012.0255740920812im,
        -11804.439063787608 - 15928.911486121006im,
        -171392.69424817772 + 80246.44673862122im,
        -179965.44649574556 + 56753.675042061im,
        -71.86136498819324 + 22.659516268995723im,
        -7.421787018608045 + 2.6698838822098505im,
        -0.0029688403973524794 + 0.0010679463471479395im,
        -7.211643408237381 + 3.1432219003998925im,
        -0.9921556840577164 - 1.9051076620073455im,
        -11.12236697804857 + 4.435065398103307im,
        5.584486257303235 + 43.750791167748744im,
        120.69404873942537 + 75.95629988805037im,
        464.02425577267593 - 141.04352156793308im,
        793.0797957198416 - 1310.5565301403608im,
        -15785.932282024054 - 21589.888379501463im,
        -231083.42622129843 + 109083.27467965364im,
        -242790.4498845943 + 77310.0624554971im,
        -97.07507558837803 + 30.90806641804309im,
        -9.383858925467502 + 3.408416729352754im,
        -0.0037529870280622714 + 0.0013631338920585844im,
        -9.114052240537562 + 4.010607801270549im,
        -1.2402367252474724 - 2.4205284629323485im,
        -14.44910647488635 + 5.303408159744719im,
        3.834094150929699 + 55.70440424800665im,
        145.30548474950302 + 104.4361093851943im,
        586.5362605410057 - 146.84860775308178im,
        1025.8795957683074 - 1603.0074869317002im,
        -19773.099793520676 - 27439.895811005757im,
        -292024.44447319285 + 139068.96872458316im,
        -307022.6605554112 + 98781.7890350056im,
        -122.78075650302557 + 39.4999712502306im,
        -11.399446068298042 + 4.18539931425596im,
        -0.004558923554064358 + 0.0016736904376566122im,
        -11.066038214769174 + 4.922063226829554im,
        -1.4866755276094077 - 2.956362893019231im,
        -17.957816340217697 + 5.974091008742144im,
        0.2828026017504354 + 67.99595327872669im,
        166.25604287937432 + 137.9354636101647im,
        712.1662544550392 - 135.74162680302868im,
        1285.571187315357 - 1887.182369064641im,
        -23766.513920437035 - 33530.990588962086im,
        -354543.2177117952 + 170514.43979302194im,
        -373032.61802776245 + 121417.4320054021im,
        -149.1819370487222 + 48.55225275123514im,
        -13.47993738897429 + 5.0093424897800904im,
        -0.005391404217339162 + 0.0020037246479805908im,
        -13.078137918202936 + 5.887282226508628im,
        -1.7307465551301384 - 3.516928810958634im,
        -21.67056238263033 + 6.3981092050115596im,
        -5.450730519272562 + 80.55250782531309im,
        182.4741751127134 + 177.28899341104957im,
        841.1945701844106 - 103.6250494134606im,
        1578.2494390381278 - 2160.6721568509524im,
        -27766.261445523523 - 39917.817751222625im,
        -418974.43635840155 + 203743.07856754382im,
        -441200.14658273524 + 145476.01466779443im,
        -176.43741060709638 + 58.17019124007889im,
        -15.637089651874174 + 5.889212368135343im,
        -0.006254434650354468 + 0.002354508783315053im,
        -15.161160848331418 + 6.9165618468619305im,
        -1.9716067015637926 - 4.106532840931153im,
        -25.60088974521228 + 6.5225176845009im,
        -13.75090019433208 + 93.22731485392029im,
        192.70508138628637 + 223.22403887031845im,
        973.6241651441433 - 45.88918002039038im,
        1910.4126194549474 - 2420.801639878026im,
        -31771.84025553049 - 46658.1940904138im,
        -485661.509853803 + 239093.51114334093im,
        -511916.3154954606 + 171229.32238085562im,
        -204.68359043366377 + 68.45357883104352im,
        -17.883047808546692 + 6.8344704334946655im,
        -0.007160723634580348 + 0.002739405936476204im,
        -17.326246050511795 + 8.02076706098771im,
        -2.2082646801774257 - 4.729539859200686im,
        -29.752130270567722 + 6.290599574782317im,
        -24.99983327398981 + 105.78452934125954im,
        195.4716331836893 + 276.3200814556965im,
        1109.1072359898747 + 42.704980363408914im,
        2289.044110263281 - 2664.573444622618im,
        -35782.05977852464 - 53813.719709690035im,
        -554958.0469773193 + 276922.5582538806im,
        -585585.4157800904 + 198964.41638898733im,
        -233.9727419029617 + 79.46981356679886im,
        -20.230358798349897 + 7.855123260556533im,
        -0.00811845120365489 + 0.0031419174823598328im,
        -19.584916355894528 + 9.211296508747864im,
        -2.439556653579119 - 5.390446802219728im,
        -34.11560810981437 + 5.642448622908839im,
        -39.568748657899576 + 117.88305533809006im,
        189.0324009813342 + 336.9545484606924im,
        1246.8559833984589 + 168.1634310310023im,
        2721.70000114854 - 2888.609964085429im,
        -39794.966236509186 - 61450.410231057984im,
        -627229.2170283698 + 317608.6261915059im,
        -662626.8860967616 + 228986.57997612408im,
        -263.9310823613192 + 91.10639107140054im,
        -20.67812544245915 + 8.051888238380782im,
        -0.008269557785417169 + 0.003213397923115062im,
        -20.013507401253204 + 9.444937635822283im,
        -2.4899919796498184 - 5.514042946446052im,
        -35.170564871956486 + 5.469427165762413im,
        -42.78565458130693 + 120.30311888032864im,
        186.3780382005341 + 350.2320756010373im,
        1273.3316680384241 + 201.10490003338106im,
        2795.4452077473843 - 2905.968804745683im,
        -40393.57847075202 - 63011.49528457448im,
        -640960.8447878574 + 325530.7471538468im,
        -677300.632280452 + 234852.14560145538im,
        -270.88577983596343 + 93.86773281320687im,
        -20.07100482975807 + 7.790355618647301im,
        -0.007999473496351418 + 0.003107746647974658im,
        -19.42730158660435 + 9.143392698642579im,
        -2.4358576737217112 - 5.341545790820189im,
        -34.30519241885582 + 5.620119554638715im,
        -39.34630237719281 + 118.1310656351643im,
        188.5831598481902 + 338.01481475522036im,
        1242.5610450911947 + 177.42780673020712im,
        2675.7409563064234 - 2822.0987494389647im,
        -39195.258269576174 - 61131.828725892934im,
        -622195.5772581355 + 315215.8835689723im,
        -657342.3197575344 + 227266.32319214108im,
        -263.9894521352796 + 91.25234367939245im,
        -17.324263260635075 + 6.639536923535341im,
        -0.0069205343629261794 + 0.002649314654582848im,
        -16.77941398417557 + 7.797613651649034im,
        -2.1444510579824034 - 4.583286391932563im,
        -29.37550509068397 + 6.048796741123717im,
        -24.924506702748513 + 106.34011895712109im,
        193.6076194286342 + 280.1151750876177im,
        1096.174517403649 + 69.34025016339902im,
        2216.5053983231955 - 2517.499998882447im,
        -34299.29619576023 - 52361.38837173969im,
        -537442.2821134315 + 269103.48781699425im,
        -567262.487197688 + 193486.87398722285im,
        -227.0856662806733 + 77.44834557670663im,
        -14.675732186274415 + 5.561809107601641im,
        -0.0058682649491631 + 0.002224342392718851im,
        -14.222129526716342 + 6.535625430607993im,
        -1.8475302611731534 - 3.8617045029109183im,
        -24.654441318519634 + 6.054479078645181im,
        -13.796240035079338 + 94.05522919071355im,
        189.5636311968375 + 229.4739426519581im,
        952.7864214947417 - 3.225655747319569im,
        1809.2340826998511 - 2192.90933981127im,
        -29403.650425568652 - 44060.08343492029im,
        -455569.5136820428 + 225761.32688815377im,
        -480449.18272122485 + 161924.81932649392im,
        -192.1886105216546 + 64.76707737721186im,
        -12.110497517025319 + 4.545529820687713im,
        -0.0048428640240512615 + 0.0018174022345118035im,
        -11.741786510897695 + 5.344037886085459im,
        -1.5462306123713947 - 3.1715799811915417im,
        -20.14379134988877 + 5.691764608228944im,
        -5.583308411402945 + 81.59299496316193im,
        178.27516434162982 + 185.55239046518986im,
        813.3894710863855 - 46.96913525033532im,
        1445.2945889602681 - 1851.9622729210926im,
        -24509.80322941619 - 36152.16389896602im,
        -376143.4864131992 + 184750.71351330582im,
        -396406.3360777614 + 132227.04925797923im,
        -158.5416463001789 + 52.87911938147627im,
        -9.614088439479238 + 3.579685369287116im,
        -0.003844874687944819 + 0.0014315851786295642im,
        -9.325051887876189 + 4.210257149477764im,
        -1.2415664495580676 - 2.5076573685844794im,
        -15.833280441714159 + 5.012263096018838im,
        0.10112577727948024 + 69.16960562377808im,
        161.294453108193 + 147.6254233739202im,
        678.6128203073627 - 67.89423279840094im,
        1116.5793191355583 - 1497.9221442961855im,
        -19618.608171233882 - 28565.99354199371im,
        -298740.74808952503 + 145654.2485485792im,
        -314651.7308015492 + 104058.6206418527im,
        -125.83865712696155 + 41.61241845049713im,
        -7.172387583811412 + 2.6537195987282867im,
        -0.002868375328620592 + 0.0010612085655415265im,
        -6.958884745240946 + 3.1222278155847025im,
        -0.934417057814179 - 1.8646934123503753im,
        -11.702428355786788 + 4.0639902243643435im,
        3.6425545312333867 + 56.91857534738799im,
        139.9415503523682 + 114.8428181096312im,
        548.8366176086113 - 71.46809773274741im,
        815.4047155211617 - 1133.7436766002854im,
        -14730.294802189599 - 21233.370714061366im,
        -222946.9344698138 + 108071.01896684374im,
        -234714.9729809388 + 77098.47425579667im,
        -93.86852896342825 + 30.830868971181964im,
        -4.771508450351986 + 1.7573750667439785im,
        -0.0019082347146658514 + 0.0007027910740619924im,
        -4.630491203570314 + 2.0681366902329708im,
        -0.6255160018699383 - 1.2375366666895296im,
        -7.722303380215889 + 2.8912032097744675im,
        5.418645000736031 + 44.9075132561788im,
        115.3395766264943 + 86.27872847750788im,
        424.2955322109813 - 62.764231749921024im,
        534.4124413967106 - 762.1366831583678im,
        -9844.587145864558 - 14088.808516786683im,
        -148355.00139137285 + 71613.00894129198im,
        -156135.13421044723 + 51036.354531597026im,
        -62.442259499561025 + 20.408837149983928im,
        -2.397794823710983 + 0.8806797896314338im,
        -0.0009589297104529257 + 0.0003521884262832514im,
        -2.3272492493670622 + 1.036546825938581im,
        -0.3155071042961153 - 0.6209894948538655im,
        -3.856962227551925 + 1.5345069357786658im,
        5.795359287916843 + 33.15372921570568im,
        88.44496526325881 + 60.97365550799961im,
        305.17561607816504 - 46.59357464147681im,
        266.4748654584399 - 385.6268247254443im,
        -4960.831713836336 - 7068.834698959753im,
        -74563.38525427971 + 35901.79746803375im,
        -78458.36032403393 + 25570.03924682769im,
        -31.37739153266395 + 10.225147065993388im,
        -0.03775882968648642 + 0.01386776720020158im,
        -1.5206810319048704e-5 + 5.37649191106773e-6im,
        -0.036648058948194034 + 0.016322135177692106im,
        -0.004968679361738421 - 0.009778646751862761im,
        -0.06487652482566891 + 0.03109460186048595im,
        5.122020628146164 + 21.639076871853472im,
        60.066986131709456 + 37.98008766953575im,
        191.7112789219719 - 27.61224504602408im,
        4.606317537246843 - 6.611945314291401im,
        -78.1232090267123 - 111.31248591642158im,
        -1174.1766426287 + 565.3361713321412im,
        -1235.508445360091 + 402.6407202613339im,
        -0.4941068909971907 + 0.16100794673854768im,
        1.7844806040045444 + 4.834163969527212im,
        0.0032462569692771935 + 0.009605170326230683im,
        14.402754824776716 + 7.816992046785759im,
        39.85140167846939 - 4.973005937840169im,
        0.09623049231305934 - 0.03566099478821205im,
        0.5362693866581878 + 1.0719645840189584im,
        0.0010646297147012323 + 0.00214049860079838im,
        3.3860986756636104 + 1.6202267625293463im,
        8.289746507400356 - 0.8773325842090969im,
        0.016637616162072585 - 0.001858633185853087im,
        0.14040384057448585 + 0.22650101157905023im,
        0.0002793087174033717 + 0.00045235122138105984im,
        0.749094528019857 + 0.3250767801935316im,
        1.660137338711462 - 0.14683562584535792im,
        0.003319027978516522 - 0.0002947153138038887im,
        0.0005573516278840598 + 0.0008975016749235774im,
        4.5603758851597413e-7 + 1.2148604691560496e-6im,
        0.002969506128344968 + 0.0012874640305961553im,
        0.006574836820047356 - 0.0005804001481589592im,
        1.3129492382349002e-5 - 1.1631785073414748e-6im,
    ]

    v1 = reshape(v1, length(v1), 1)      # reshaping because in system it must be a 2D array

    ############################
    ctsys.data.paramsoptical.eigenvalues = λ1
    ctsys.data.paramsoptical.eigenvectors = v1
    ctsys.data.paramsoptical.power = p.P[2]

    previousSolution = currentSolution
    solution = solve(ctsys; inival = inival, control = control)
    ctsys.data.paramsoptical.oldSolution = solution
    currentSolution = solution
    inival = solution


    ctsys.data.paramsoptical.eigenvalues = λ1
    ctsys.data.paramsoptical.eigenvectors = v1
    ctsys.data.paramsoptical.power = p.P[2]

    previousSolution = currentSolution
    solution = solve(ctsys; inival = inival, control = control)
    ctsys.data.paramsoptical.oldSolution = solution
    currentSolution = solution
    inival = solution


    if plotting
        vis = GridVisualizer(; Plotter = PyPlot, fignumber = 2, resolution = (1200, 900))

        scalarplot!(
            vis, grid, solution[1, :], Plotter = PyPlot, legend = :best, clear = false, title = "Applied voltage Δu = $maxBias V",
            xlabel = "cross section space along \$x=0\$ [m]", ylabel = "potentials [V]", fontsize = 55, linewidth = 5,
            slice = :x => 0, label = "\$ \\varphi_n \$", color = "blue"
        )

        scalarplot!(
            vis, grid, solution[2, :], Plotter = PyPlot, linewidth = 5,
            slice = :x => 0, label = "\$ \\varphi_p \$", clear = false, color = "mediumvioletred"
        )

        scalarplot!(
            vis, grid, solution[3, :], Plotter = PyPlot, linewidth = 5,
            slice = :x => 0, label = "\$ \\psi \$", clear = false, color = "darkorange"
        )
    end


    testval = sum(solution) / length(solution)
    return testval

end # main

function test()
    testval = 0.5122451923673309
    return main(test = true) ≈ testval
end

if test == false
    println("This message should show when this module is successfully recompiled.")
end

end # module
