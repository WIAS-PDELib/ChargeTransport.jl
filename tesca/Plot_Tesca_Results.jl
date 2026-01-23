#=
In this file we read in the results from Tesca and
use julia plotting routines to plot the results.
=#

module Plot_Tesca_Results

using ChargeTransport
using ExtendableGrids

function main(; Plotter = nothing)

    Plotter.pygui(true) # Plots as Pop-Ups

    # Read in grid from dom file
    grid = simplexgrid("tesca.dom", format = "dom")
    grid[Coordinates] *= 1.0e-6
    trim!(grid)


    # Read in data from diodat file
    data_eq = read_diodat("tesca-mosS-eq_dio.dat")
    data = read_diodat("tesca-mosS_dio.dat")

    ################################################################################
    # Solution psi with PythonPlot
    ################################################################################
    solution_psi_eq = data_eq["ElectrostaticPotential"]

    XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

    fig = Plotter.figure()
    ax = fig.add_subplot(111, projection = "3d")
    ax.plot_trisurf(XX[:], YY[:], solution_psi_eq)

    Plotter.title("Electrostatic Potential Equilibrium", fontsize = 16)
    Plotter.xlabel("Width [m]", fontsize = 12)
    Plotter.ylabel("Height [m]", fontsize = 12)
    Plotter.zlabel("potential [V]", fontsize = 12)

    ax = Plotter.gca()
    ax.tick_params(axis = "x", labelsize = 10)
    ax.tick_params(axis = "y", labelsize = 10)
    ax.tick_params(axis = "z", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    ax.zaxis.get_offset_text().set_fontsize(10)

    Plotter.tight_layout()
    current_figure = Plotter.gcf()
    display(current_figure)


    ################################################################################
    # Solution psi with PythonPlot
    ################################################################################
    solution_psi = data["ElectrostaticPotential"]

    XX = grid[Coordinates][1, :]; YY = grid[Coordinates][2, :]

    fig = Plotter.figure()
    ax = fig.add_subplot(111, projection = "3d")
    ax.plot_trisurf(XX[:], YY[:], solution_psi)

    Plotter.title("Electrostatic Potential", fontsize = 16)
    Plotter.xlabel("Width [m]", fontsize = 12)
    Plotter.ylabel("Height [m]", fontsize = 12)
    Plotter.zlabel("potential [V]", fontsize = 12)

    ax = Plotter.gca()
    ax.tick_params(axis = "x", labelsize = 10)
    ax.tick_params(axis = "y", labelsize = 10)
    ax.tick_params(axis = "z", labelsize = 10)
    ax.xaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontsize(10)
    ax.zaxis.get_offset_text().set_fontsize(10)

    Plotter.tight_layout()
    current_figure = Plotter.gcf()
    display(current_figure)

    return nothing
end

end
