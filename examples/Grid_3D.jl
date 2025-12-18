#=
# Example code for a 3D grid.
([source code](SOURCE_URL))

This code provides an unstructured grid for a three-dimensional domain with an additional
quadratic hole within the device. The grid is produced with TetGen.jl.
=#

module Grid_3D

using ExtendableGrids

using ChargeTransport
## For using this example, one additionally needs to add TetGen. SimplexGridFactory is a wrapper for using this meshgenerator.
using SimplexGridFactory
using TetGen

function main(; Plotter = nothing) # plotting is currently only tested with GLMakie and PythonPlot

    @local_unitfactors μm cm s ns V K ps Hz W

    builder3d = let

        b = SimplexGridBuilder(Generator = TetGen)

        ## region numbers
        regionDonor = 1                           # n doped region
        regionIntrinsic = 2                           # intrinsic region
        regionAcceptor = 3                           # p doped region

        ## boundary region numbers
        bregionDonor = 1
        bregionAcceptor = 2
        bregionJ1 = 3
        bregionJ2 = 4
        bregionNoFlux = 5

        ## grid
        h_ndoping = 9.9e-6 * cm
        h_intrinsic = 4.0e-5 * cm + 2.0e-7 * cm
        h_pdoping = 1.99e-5 * cm
        height = 3.0e-5 * cm
        width = 3.0e-5 * cm

        ## lower area (nregion)
        Area1length_0 = point!(b, 0.0, 0.0, 0.0)
        Area1length_n = point!(b, h_ndoping, 0.0, 0.0)
        Area1width_n = point!(b, h_ndoping, width, 0)
        Area1width_0 = point!(b, 0, width, 0)

        ## upper area (nregion)
        Area2length_0 = point!(b, 0.0, 0.0, height)
        Area2length_n = point!(b, h_ndoping, 0.0, height)
        Area2height_n = point!(b, h_ndoping, width, height)
        Area2height_0 = point!(b, 0, width, height)

        ## lower area (iregion)
        Area1length_ni = point!(b, h_ndoping + h_intrinsic, 0.0, 0.0)
        Area1width_ni = point!(b, h_ndoping + h_intrinsic, width, 0)
        ## upper area (iregion)
        Area2length_ni = point!(b, h_ndoping + h_intrinsic, 0.0, height)
        Area2height_ni = point!(b, h_ndoping + h_intrinsic, width, height)

        ## lower area (pregion)
        Area1length_nip = point!(b, h_ndoping + h_intrinsic + h_pdoping, 0.0, 0.0)
        Area1width_nip = point!(b, h_ndoping + h_intrinsic + h_pdoping, width, 0)
        ## upper area (pregion)
        Area2length_nip = point!(b, h_ndoping + h_intrinsic + h_pdoping, 0.0, height)
        Area2height_nip = point!(b, h_ndoping + h_intrinsic + h_pdoping, width, height)

        ## n-region
        facetregion!(b, bregionNoFlux) # surface below
        facet!(b, Area1length_0, Area1length_n, Area1width_n, Area1width_0)
        facetregion!(b, bregionNoFlux) # surface up
        facet!(b, Area2length_0, Area2length_n, Area2height_n, Area2height_0)
        facetregion!(b, bregionNoFlux) # surface front
        facet!(b, Area1length_0, Area1length_n, Area2length_n, Area2length_0)
        facetregion!(b, bregionNoFlux) # surface back
        facet!(b, Area1width_n, Area1width_0, Area2height_0, Area2height_n)

        ## i-region
        facetregion!(b, bregionNoFlux) # surface below
        facet!(b, Area1length_n, Area1length_ni, Area1width_ni, Area1width_n)
        facetregion!(b, bregionNoFlux) # surface up
        facet!(b, Area2length_n, Area2length_ni, Area2height_ni, Area2height_n)
        facetregion!(b, bregionNoFlux) # surface front
        facet!(b, Area1length_n, Area1length_ni, Area2length_ni, Area2length_n)
        facetregion!(b, bregionNoFlux) # surface back
        facet!(b, Area1width_ni, Area1width_n, Area2height_n, Area2height_ni)

        ## p-region
        facetregion!(b, bregionNoFlux) # untere Oberfläche
        facet!(b, Area1length_ni, Area1length_nip, Area1width_nip, Area1width_ni)
        facetregion!(b, bregionNoFlux) # surface up
        facet!(b, Area2length_ni, Area2length_nip, Area2height_nip, Area2height_ni)
        facetregion!(b, bregionNoFlux) # surface front
        facet!(b, Area1length_ni, Area1length_nip, Area2length_nip, Area2length_ni)
        facetregion!(b, bregionNoFlux) # surface back
        facet!(b, Area1width_nip, Area1width_ni, Area2height_ni, Area2height_nip)

        ## inner interfaces
        facetregion!(b, bregionJ1) # inner interface n/i
        facet!(b, Area1length_n, Area1width_n, Area2height_n, Area2length_n)
        facetregion!(b, bregionJ2) # inner interface i/n
        facet!(b, Area1length_ni, Area1width_ni, Area2height_ni, Area2length_ni)

        facetregion!(b, bregionDonor) # metalinterface left
        facet!(b, Area1width_0, Area1length_0, Area2length_0, Area2height_0)
        facetregion!(b, bregionAcceptor) # metalinterface right
        facet!(b, Area1width_nip, Area1length_nip, Area2length_nip, Area2height_nip)

        distance = 8.0e-6 * cm

        hp1 = point!(b, h_ndoping + h_intrinsic / 2 - distance, width / 2 - distance, height / 2 - distance)
        hp2 = point!(b, h_ndoping + h_intrinsic / 2 + distance, width / 2 - distance, height / 2 - distance)
        hp3 = point!(b, h_ndoping + h_intrinsic / 2 + distance, width / 2 + distance, height / 2 - distance)
        hp4 = point!(b, h_ndoping + h_intrinsic / 2 - distance, width / 2 + distance, height / 2 - distance)
        hp5 = point!(b, h_ndoping + h_intrinsic / 2 - distance, width / 2 - distance, height / 2 + distance)
        hp6 = point!(b, h_ndoping + h_intrinsic / 2 + distance, width / 2 - distance, height / 2 + distance)
        hp7 = point!(b, h_ndoping + h_intrinsic / 2 + distance, width / 2 + distance, height / 2 + distance)
        hp8 = point!(b, h_ndoping + h_intrinsic / 2 - distance, width / 2 + distance, height / 2 + distance)

        facetregion!(b, 6)
        facet!(b, hp1, hp2, hp3, hp4)
        facet!(b, hp5, hp6, hp7, hp8)
        facet!(b, hp1, hp2, hp6, hp5)
        facet!(b, hp2, hp3, hp7, hp6)
        facet!(b, hp3, hp4, hp8, hp7)
        facet!(b, hp4, hp1, hp5, hp8)
        holepoint!(b, h_ndoping + h_intrinsic / 2, width / 2, height / 2)

        cellregion!(b, regionDonor)
        regionpoint!(b, h_ndoping - 2.0e-6 * cm, width / 2 - 2.0e-6 * cm, height / 2 - 2.0e-6 * cm)
        cellregion!(b, regionIntrinsic)
        regionpoint!(b, h_ndoping + h_intrinsic - 6.0e-6 * cm, width / 2 - 6.0e-6 * cm, height / 2 - 6.0e-6 * cm)
        cellregion!(b, regionAcceptor)
        regionpoint!(b, h_ndoping + h_intrinsic + h_pdoping - 6.0e-6 * cm, width / 2 - 6.0e-6 * cm, height / 2 - 6.0e-6 * cm)

        options!(b, maxvolume = 1.0e-24)

        b

    end

    grid = simplexgrid(builder3d)

    if Plotter !== nothing # plotting is currently only tested with GLMakie and PythonPlot
        gridplot(grid; Plotter, zplane = 1.5e-7, azim = 20, elev = 60, linewidth = 0.5, legend = :lt)
    end

    return grid

end # main

end # module
