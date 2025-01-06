#=
# Example code for a 2D non rectangular grid.
([source code](SOURCE_URL))

This code provides an unstructured grid for a non rectangular two-dimensional domain.
The grid is produced with Triangulate.jl.
=#
ENV["LC_NUMERIC"] = "C"

module Non_RectangularGrid_2D

    using ChargeTransport
    using ExtendableGrids
    using PyPlot

    ## For using this example, one additionally needs to add Triangulate. SimplexGridFactory is a wrapper for using this meshgenerator.
    using SimplexGridFactory
    using Triangulate

    function main(; Plotter = PyPlot, plotting = false)

        ## region numbers
        regionDonor = 1                           # n doped region
        regionIntrinsic = 2                           # intrinsic region
        regionAcceptor = 3                           # p doped region

        ## boundary region numbers
        bregionDonor = 1
        bregionAcceptor = 2
        bregionJunction1 = 3
        bregionJunction2 = 4
        bregionNoFlux = 5

        ## grid
        h_ndoping = 9.9e-6 * cm
        h_intrinsic = 4.0e-5 * cm + 2.0e-7 * cm
        h_pdoping = 1.99e-5 * cm
        height = 3.0e-5 * cm

        function unsuitable(x1, y1, x2, y2, x3, y3, area)
            bary_x = (x1 + x2 + x3) / 3.0
            bary_y = (y1 + y2 + y3) / 3.0
            dx = bary_x - refinement_center[1]
            dy = bary_y - refinement_center[2]
            qdist = dx^2 + dy^2
            return area > 0.1 * max(8.0e-16, qdist)
        end

        b = SimplexGridBuilder(Generator = Triangulate)

        ## specify boundary nodes
        length_0 = point!(b, 0.0, 0.0)
        length_n = point!(b, h_ndoping, 0.0)
        length_ni = point!(b, h_ndoping + h_intrinsic, 0.0)
        length_nip = point!(b, h_ndoping + h_intrinsic + h_pdoping, 0.0)
        height_0 = point!(b, 0.0, height)
        height_n = point!(b, h_ndoping, height)

        ## for L shape
        height_ni12 = point!(b, h_ndoping + h_intrinsic / 2, height)
        height_ni2 = point!(b, h_ndoping + h_intrinsic / 2, height / 2)
        height_ni = point!(b, h_ndoping + h_intrinsic, height / 2)
        height_nip = point!(b, h_ndoping + h_intrinsic + h_pdoping, height / 2)

        ## specify boundary regions
        ## metal interface
        facetregion!(b, bregionDonor)
        facet!(b, length_0, height_0)
        facetregion!(b, bregionAcceptor)
        facet!(b, length_nip, height_nip)

        ## no flux
        facetregion!(b, bregionNoFlux)
        facet!(b, length_0, length_nip)
        facetregion!(b, bregionNoFlux)
        facet!(b, height_0, height_n)
        facetregion!(b, bregionNoFlux)
        facet!(b, height_0, height_ni12)
        facetregion!(b, bregionNoFlux)
        facet!(b, height_ni12, height_ni2)
        facetregion!(b, bregionNoFlux)
        facet!(b, height_ni2, height_nip)

        ## inner interface
        facetregion!(b, bregionJunction1)
        facet!(b, length_n, height_n)
        facetregion!(b, bregionJunction2)
        facet!(b, length_ni, height_ni)

        refinement_center = [h_ndoping + h_intrinsic / 2, height / 2]
        ## Activate unsuitable callback
        options!(b, unsuitable = unsuitable)

        ## cell regions
        cellregion!(b, regionDonor)
        regionpoint!(b, h_ndoping - 1.0e-6 * cm, height / 2 - 1.0e-6 * cm)
        cellregion!(b, regionIntrinsic)
        regionpoint!(b, h_ndoping + h_intrinsic - 1.0e-6 * cm, height / 2 - 1.0e-6 * cm)
        cellregion!(b, regionAcceptor)
        regionpoint!(b, h_ndoping + h_intrinsic + h_pdoping - 1.0e-6 * cm, height / 2 - 1.0e-6 * cm)

        options!(b, maxvolume = 7.0e-16)

        grid = simplexgrid(b)

        if plotting
            gridplot(grid, Plotter = Plotter, resolution = (600, 400), linewidth = 0.6)
            Plotter.xlabel("length [m]")
            Plotter.ylabel("height [m]")
            Plotter.tight_layout()
            # builderplot(b,Plotter=Plotter,resolution=(750,700))
        end

        return nothing

    end # main

end # module
