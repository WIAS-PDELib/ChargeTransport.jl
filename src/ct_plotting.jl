"""
$(TYPEDSIGNATURES)
Method which can be used to construct the arrays parsed to the plotting routines for labeling.
The description for electrons and holes are predefined. If one wishes to extend by labels for,
e.g. mobile ionic carriers or traps, this can be done within the main file.

"""

function set_plotting_labels(data)

    label_energy = Matrix{LaTeXString}(undef, 2, data.params.numberOfCarriers) # band-edge energies and potential
    label_BEE = Vector{LaTeXString}(undef, data.params.numberOfCarriers)    # band-edge energie parameters
    label_density = Vector{LaTeXString}(undef, data.params.numberOfCarriers)
    label_solution = Vector{LaTeXString}(undef, data.params.numberOfCarriers)

    # indices (∈ IN) of electron and hole quasi Fermi potentials specified by user
    iphin = data.bulkRecombination.iphin # integer index of φ_n
    iphip = data.bulkRecombination.iphip # integer index of φ_p

    ## for electrons
    label_energy[1, iphin] = L"E_c-q\psi"; label_energy[2, iphin] = L"- q \varphi_n"; label_BEE[iphin] = L"E_c"
    label_density[iphin] = L"n_n";              label_solution[iphin] = L"\varphi_n"

    ## for holes
    label_energy[1, iphip] = L"E_v-q\psi"; label_energy[2, iphip] = L"- q \varphi_p"; label_BEE[iphip] = L"E_v"
    label_density[iphip] = L"n_p";              label_solution[iphip] = L"\varphi_p"

    return label_solution, label_density, label_energy, label_BEE
end

"""
$(TYPEDSIGNATURES)
Plotting routine, where the charge carrier densities are depicted
in dependence of space. The case of heterojunctions is tested, but yet
multidimensional plottings are not included.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.

"""
function plot_densities!(visualizer, ctsys, solution, title, label_density, ; plotGridpoints = false)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data
    numberOfRegions = grid[NumCellRegions]

    if dim_space(grid) > 1
        error("plot_densities is so far only tested in 1D")
    end

    if plotGridpoints == true
        marker = :circle
    else
        marker = :none
    end

    params = data.params
    colors = ["green", "red", "gold", "purple", "orange"]
    linestyles = [:solid, :dot, :dash, :dashdot, :solid]

    for icc in 1:params.numberOfCarriers

        label_is_plotted = false

        for ireg in 1:numberOfRegions
            subg = subgrid(grid, [ireg])
            ncc = get_density(solution, ireg, ctsys, icc)

            ## Note that this implies a 1D plot, for multidimensional plots, you may work with
            ## GridVisualize.jl or write your own code.
            scalarplot!(
                visualizer,
                subg,
                1.0e-6 .* ncc;
                clear = false,
                color = colors[icc],
                label = label_is_plotted ? nothing : label_density[icc],
                legend = :cc,
                linestyle = linestyles[icc],
                linewidth = 3,
                title = title,
                xlabel = L"\text{space [m]}",
                ylabel = L"density [$\frac{1}{\text{cm}^3}$]",
                yscale = :log
            )
            label_is_plotted = true

        end
    end

    return nothing
end

function plot_densities(Plotter, ctsys, solution, title, label_density, ; plotGridpoints = false)
    @warn "plot_densities() is deprecated, please use plot_densities!() with a GridVisualizer"

    vis = GridVisualizer(Plotter = Plotter)
    plot_densities!(vis, ctsys, solution, title, label_density, ; plotGridpoints)

    return reveal(vis)
end

"""
$(TYPEDSIGNATURES)

With this method it is possible to plot the energies

``E_\\alpha - q \\psi \\quad \\text{w.r.t. space.}``

The case of heterojunctions is tested, but yet
multidimensional plottings are not included.

One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.

"""

function plot_energies!(visualizer, ctsys, solution, title, label_energy; plotGridpoints = false)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data
    coord = grid[Coordinates]

    (; q) = data.constants

    if length(coord[1]) != 1
        println("plot_energies is so far only implemented in 1D")
    end

    if plotGridpoints
        marker = :circle
    else
        marker = :none
    end

    colors = ["green", "red", "gold", "purple", "orange"]
    linestyles = [:solid, :dot, :dash, :dashdot, :solid]

    for icc in data.electricCarrierList

        #grids = Array{ExtendableGrid, 1}(undef, numberOfRegions)
        #nicc  = Array{Array{Float64, 1}, 1}(undef, numberOfRegions)

        label_is_plotted = false

        for ireg in 1:data.params.numberOfRegions
            subg = subgrid(grid, [ireg])
            Ecc = get_BEE(icc, ireg, ctsys)
            solpsi = view(solution[data.index_psi, :], subg)
            solcc = view(solution[icc, :], subg)

            ## Note that this implies a 1D plot, for multidimensional plots, you may work with
            ## GridVisualize.jl or write your own code.
            scalarplot!(
                visualizer,
                subg,
                Ecc ./ q .- solpsi;
                markershape = marker,
                title = title,
                xlabel = L"\text{space [m]}",
                ylabel = L"\text{energy [eV]}",
                label = label_is_plotted ? nothing : label_energy[1, icc],
                legend = :cc,
                markersize = 8,
                linewidth = 2,
                color = colors[icc],
                linestyle = linestyles[1],
                clear = false
            )

            scalarplot!(
                visualizer,
                subg,
                - solcc,
                markershape = marker,
                label = label_is_plotted ? nothing : label_energy[2, icc],
                legend = :cc,
                markersize = 8,
                linewidth = 2,
                color = colors[icc],
                linestyle = linestyles[2],
                clear = false
            )

            label_is_plotted = true
        end

    end

    for iicc in data.ionicCarrierList
        icc = iicc.ionicCarrier

        label_is_plotted = false

        for ireg in 1:data.params.numberOfRegions
            if ireg ∈ iicc.regions
                subg = subgrid(grid, [ireg])
                Ecc = get_BEE(icc, ireg, ctsys)
                solpsi = view(solution[data.index_psi, :], subg)
                solcc = view(solution[icc, :], subg)

                ## Note that this implies a 1D plot, for multidimensional plots, you may work with
                ## GridVisualize.jl or write your own code.
                scalarplot!(
                    visualizer,
                    subg,
                    data.params.chargeNumbers[icc] .* (solpsi .- Ecc ./ q),
                    label = label_is_plotted ? nothing : label_energy[1, icc],
                    legend = :cc,
                    markershape = marker,
                    linewidth = 2,
                    color = colors[icc],
                    linestyle = linestyles[1]
                )

                scalarplot!(
                    visualizer,
                    subg,
                    data.params.chargeNumbers[icc] .* solcc,
                    label = label_is_plotted ? nothing : label_energy[2, icc],
                    legend = :cc,
                    markershape = marker,
                    linewidth = 2,
                    color = colors[icc],
                    linestyle = linestyles[2]
                )

                label_is_plotted = true
            end
        end
    end


    for iicc in data.trapCarrierList
        icc = iicc.trapCarrier
        count = 0

        for ireg in 1:data.params.numberOfRegions
            if ireg ∈ iicc.regions
                subg = subgrid(grid, [ireg])
                Ecc = get_BEE(icc, ireg, ctsys)
                solpsi = view(solution[data.index_psi, :], subg)
                solcc = view(solution[icc, :], subg)

                if count == 0
                    label1 = label_energy[1, icc]
                    label2 = label_energy[2, icc]
                else
                    label1 = ""
                    label2 = ""
                end
                ## Note that this implies a 1D plot, for multidimensional plots, you may work with
                ## GridVisualize.jl or write your own code.
                scalarplot!(
                    visualizer,
                    subg,
                    data.params.chargeNumbers[icc] .* (solpsi .- Ecc ./ q);
                    label = label1,
                    legend = :best,
                    markershape = marker,
                    linewidth = 2,
                    color = colors[icc],
                    linestyle = linestyles[1]
                )

                scalarplot!(
                    visualizer,
                    subg,
                    data.params.chargeNumbers[icc] .* solcc;
                    label = label2,
                    legend = :best,
                    markershape = marker,
                    linewidth = 2,
                    color = colors[icc],
                    linestyle = linestyles[2]
                )
                count = count + 1
            end
        end
    end

    return nothing
end

function plot_energies(Plotter, ctsys, solution, title, label_energy, ; plotGridpoints = false)
    @warn "plot_energies() is deprecated, please use plot_energies!() with a GridVisualizer"

    vis = GridVisualizer(Plotter = Plotter)
    plot_energies!(vis, ctsys, solution, title, label_energy, ; plotGridpoints)

    return reveal(vis)
end

"""
$(SIGNATURES)
With this method it is possible to depict the band-edge energies ``E_\\alpha ``.
This can be useful for debugging when dealing with heterojunctions.

"""
function plot_energies!(visualizer, ctsys, label_BEE)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data
    params = data.params
    paramsnodal = data.paramsnodal

    q = data.constants.q

    coord = grid[Coordinates]
    cellregions = grid[CellRegions]
    cellnodes = grid[CellNodes]

    if size(coord, 1) != 1
        error("plot_energies is so far only implemented in 1D")
    end

    colors = ["green", "red", "gold", "purple", "orange"]
    linestyles = [:solid, :dot, :dash, :dashdot, :solid]

    #plot different band-edge energies values in interior
    for icc in 1:params.numberOfCarriers

        label_is_plotted = false

        for i in eachindex(cellregions)
            # determine band-edge energy value in cell and number of cell nodes
            cellValue = (params.bandEdgeEnergy[icc, cellregions[i]] + paramsnodal.bandEdgeEnergy[icc, i]) / q
            numberLocalCellNodes = length(cellnodes[:, i])

            # patch together cells
            scalarplot!(
                visualizer,
                coord[cellnodes[:, i]][:],
                fill(cellValue, numberLocalCellNodes);
                title = "Band-edge energies",
                xlabel = L"\text{space [m]}",
                ylabel = L"\text{energy [eV]}",
                label = label_is_plotted ? nothing : label_BEE[icc],
                legend = :rc,
                clear = false,
                markershape = :cross,
                color = colors[icc],
                linewidth = 3,
                linestyle = linestyles[icc],
            )

            label_is_plotted = true
        end

    end

    #plot different band-edge energy values on boundary
    bfaceregions = grid[BFaceRegions]
    bfacenodes = grid[BFaceNodes]

    for icc in 1:params.numberOfCarriers
        for i in eachindex(bfaceregions)
            # determine band-edge energy value in cell and number of cell nodes
            cellValue = (params.bBandEdgeEnergy[icc, bfaceregions[i]] + paramsnodal.bandEdgeEnergy[icc, bfacenodes[i]]) / q
            numberLocalCellNodes = length(bfacenodes[:, i])

            # patch together cells
            scalarplot!(
                visualizer,
                coord[bfacenodes[:, i]],
                fill(cellValue, numberLocalCellNodes);
                clear = false,
                markershape = :cross,
                color = colors[icc]
            )
        end

    end

    return nothing
end

function plot_energies(Plotter, ctsys, label_BEE)
    @warn "plot_energies() is deprecated, please use plot_energies!() with a GridVisualizer"

    vis = GridVisualizer(Plotter = Plotter)
    plot_energies!(vis, ctsys, label_BEE)

    return reveal(vis)
end

"""
$(TYPEDSIGNATURES)
Possibility to plot the considered doping. This is especially useful
for making sure that the interior and the boundary doping agree.

"""
function plot_doping!(visualizer, ctsys, label_density)

    g = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data
    params = data.params
    coord = g[Coordinates]
    cellregions = g[CellRegions]
    cellnodes = g[CellNodes]
    coord = g[Coordinates]

    if length(coord[1]) != 1
        error("plot_doping is so far only implemented in 1D")
    end

    colors = ["green", "red", "gold", "purple", "orange"]
    linestyles = [:solid, :dot, :dash, :dashdot, :solid]

    (ymin, ymax) = 1.0e-6 .* extrema(params.doping)
    ymax *= 10

    # plot different doping values in interior
    for icc in 1:params.numberOfCarriers

        for i in eachindex(cellregions)
            # determine doping value in cell and number of cell nodes
            cellValue = params.doping[icc, cellregions[i]]
            numberLocalCellNodes = length(cellnodes[:, i])


            # patch together cells (multiplying by 1.0e-6 gives us the densities in cm^(-3))
            scalarplot!(
                visualizer,
                coord[cellnodes[:, i]],
                1.0e-6 .* fill(cellValue, numberLocalCellNodes);
                clear = false,
                color = colors[icc],
                label = i == firstindex(cellregions) ? label_density[icc] : nothing,
                legend = :rc,
                linestyle = linestyles[icc],
                linewidth = 3,
                title = "Doping values for charge carriers",
                xlabel = L"\text{space [m]}",
                ylabel = L"doping [$\frac{1}{\text{cm}^3}$]",
                limits = (ymin, ymax),
                yscale = :symlog
            )
        end

    end

    # plot different doping values on boundary
    bfaceregions = g[BFaceRegions]
    bfacenodes = g[BFaceNodes]

    for icc in 1:params.numberOfCarriers

        for i in eachindex(bfaceregions)
            # determine doping value in cell and number of cell nodes
            cellValue = params.bDoping[icc, bfaceregions[i]]
            numberLocalCellNodes = length(bfacenodes[:, i])

            # patch together cells
            scalarplot!(
                visualizer,
                fill(coord[bfacenodes[:, i]]..., 2),
                1.0e-6 .* fill(cellValue, numberLocalCellNodes + 1),
                clear = false,
                markershape = :xcross,
                markersize = 20,
                markevery = 1,
                color = colors[icc],
                limits = (ymin, ymax),
                xladel = L"\text{space [m]}",
                ylabel = L"doping [$\frac{1}{\text{cm}^3}$]",
                yscale = :symlog
            )
        end

    end

    return nothing
end

function plot_doping(Plotter, ctsys, label_density)
    @warn "plot_doping() is deprecated, please use plot_doping!() with a GridVisualizer"

    vis = GridVisualizer(Plotter = Plotter)
    plot_doping!(vis, ctsys, label_density)

    return reveal(vis)
end

"""
Plot doping for nodal dependent doping.
"""

function plot_doping!(visualizer, g::ExtendableGrid, paramsnodal::ParamsNodal)

    coord = g[Coordinates]

    doping = 1.0e-6 .* paramsnodal.doping[:]

    doping[doping .<= 0] .= NaN

    scalarplot!(
        visualizer,
        coord[:],
        doping,
        title = "Doping values for charge carriers",
        xlabel = L"\text{space [m]}",
        ylabel = L"doping [$\frac{1}{\text{cm}^3}$]",
        yscale = :log,
        color = :green,
        markershape = :cross,
        markersize = 8
    )

    return nothing

end

function plot_doping(Plotter, g::ExtendableGrid, paramsnodal::ParamsNodal)
    @warn "plot_doping() is deprecated, please use plot_doping!() with a GridVisualizer"

    vis = GridVisualizer(Plotter = Plotter)
    plot_doping!(vis, g::ExtendableGrid, paramsnodal::ParamsNodal)

    return reveal(vis)
end

"""
$(TYPEDSIGNATURES)
Plotting routine for depicting the electroneutral potential.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plot_electroNeutralSolutionBoltzmann!(visualizer, grid, psi0; plotGridpoints = false)

    if plotGridpoints
        marker = :circle
    else
        marker = :none
    end

    coord = grid[Coordinates]

    scalarplot!(
        visualizer,
        coord[:],
        psi0,
        title = "Electroneutral potential",
        xlabel = L"\text{space [m]}",
        ylabel = L"\text{potential [V]}",
        color = :blue,
        markershape = marker,
        markersize = 8
    )

    return nothing
end

function plot_electroNeutralSolutionBoltzmann(Plotter, grid, psi0; plotGridpoints = false)
    @warn "plot_electroNeutralSolutionBoltzmann() is deprecated, please use plot_electroNeutralSolutionBoltzmann!() with a GridVisualizer"

    vis = GridVisualizer(Plotter = Plotter)
    plot_electroNeutralSolutionBoltzmann!(vis, grid, psi0; plotGridpoints)

    return reveal(vis)
end

"""
$(TYPEDSIGNATURES)
Method for plotting the solution vectors: the electrostatic potential ``\\psi``
as well as the charge carriers.
The case of heterojunctions is tested, but yet
multidimensional plottings are not included.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plot_solution!(visualizer, ctsys, solution, title, label_solution; plotGridpoints = false)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    if dim_space(grid) > 1
        error("plot_solution is so far only tested in 1D")
    end

    if plotGridpoints == true
        marker = :circle
    else
        marker = :none
    end

    coord = grid[Coordinates]'
    ipsi = data.index_psi

    colors = ["green", "red", "gold", "purple", "orange"]
    linestyles = [:solid, :dot, :dash, :dashdot, :solid]

    scalarplot!(
        visualizer,
        grid,
        solution[ipsi, :];
        clear = false,
        color = :blue,
        label = L"\psi",
        legend = :rc,
        linewidth = 3,
        markershape = marker,
        markersize = 8,
        title = title,
        xlabel = L"\text{space [m]}",
        ylabel = L"\text{potential [V]}"
    )

    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOn
        ipsiStandard = data.barrierLoweringInfo.ipsiStandard

        scalarplot!(
            visualizer,
            grid,
            solution[ipsiStandard, :];
            clear = false,
            color = :black,
            label = L"$\psi$ (Schottky contacts)",
            legend = :cc,
            linestyle = :dot,
            linewidth = 3,
            markershape = marker,
            markersize = 8
        )
    end

    # electrons and holes
    for icc in data.electricCarrierList

        scalarplot!(
            visualizer,
            grid,
            solution[icc, :];
            clear = false,
            color = colors[icc],
            label = label_solution[icc],
            legend = :cc,
            linestyle = linestyles[1],
            linewidth = 3,
            markershape = marker,
            markersize = 8
        )
    end

    for icc in data.ionicCarrierList
        # DA: could be modified by using subgrids from ExtendableGrids.
        # this is needed to only plot present ionic charge carrier in respective defined regions
        regions = grid[CellRegions]

        subregions = zeros(Int64, 0)
        for ix in 1:length(icc.regions)
            subreg = findall(x -> x == icc.regions[ix], regions)
            append!(subregions, subreg)
            push!(subregions, subregions[end] + 1)
        end
        subgrid = coord[subregions]

        icc = icc.ionicCarrier

        scalarplot!(
            visualizer,
            subgrid ./ 1,
            solution[icc, subregions];
            clear = false,
            color = colors[icc],
            label = label_solution[icc],
            legend = :cc,
            linestyle = linestyles[1],
            linewidth = 3,
            markershape = marker,
            markersize = 8
        )
    end
    for icc in data.trapCarrierList
        # DA: could be modified by using subgrids from ExtendableGrids.
        # this is needed to only plot present ionic charge carrier in respective defined regions
        regions = grid[CellRegions]

        subregions = zeros(Int64, 0)
        for it in 1:length(icc.regions)
            subreg = findall(x -> x == icc.regions[it], regions)
            append!(subregions, subreg)
            push!(subregions, subregions[end] + 1)
        end
        subgrid = coord[subregions]

        icc = icc.trapCarrier

        Plotter.plot(subgrid ./ 1, solution[icc, subregions], label = label_solution[icc], marker = marker, color = colors[icc], linestyle = linestyles[1], linewidth = 3)
    end

    return nothing

end

function plot_solution(Plotter, ctsys, solution, title, label_solution; plotGridpoints = false)
    @warn "plot_solution() is deprecated, please use plot_solution!() with a GridVisualizer"

    vis = GridVisualizer(Plotter = Plotter)
    plot_solution!(vis, ctsys, solution, title, label_solution; plotGridpoints)

    return reveal(vis)
end

"""
$(TYPEDSIGNATURES)
Method for showing the total current.
One input parameter is the boolean plotGridpoints which makes it possible to plot markers,
which indicate where the nodes are located.
"""
function plot_IV!(visualizer, biasValues, IV, title, ; plotGridpoints = false)

    if plotGridpoints == true
        marker = :circle
    else
        marker = :none
    end

    scalarplot!(
        visualizer,
        biasValues[1:length(IV)],
        IV;
        color = :blue,
        markershape = marker,
        markersize = 8,
        title = title,
        xlabel = L"\text{bias [V]}",
        ylabel = L"\text{total current [A]}"
    )

    return nothing
end

function plot_IV(Plotter, biasValues, IV, title, ; plotGridpoints = false)
    @warn "plot_IV() is deprecated, please use plot_IV!() with a GridVisualizer"

    vis = GridVisualizer(Plotter = Plotter)
    plot_IV!(visualizer, biasValues, IV, title, ; plotGridpoints)

    return reveal(vis)
end
