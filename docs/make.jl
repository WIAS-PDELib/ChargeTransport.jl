using Documenter, ExampleJuggler, Literate, ChargeTransport, PlutoStaticHTML, PlutoSliderServer, VoronoiFVM

DocMeta.setdocmeta!(ExampleJuggler, :DocTestSetup, :(using ExampleJuggler); recursive = true)

exampledir = joinpath(@__DIR__, "..", "examples")
notebookdir = joinpath(@__DIR__, "..", "notebooks")

modules = modules = [
    "Ex101_PIN.jl",
    "Ex102_PIN_nodal_doping.jl",
    "Ex103_PSC_IVMeasurement.jl",
    "Ex104_PSC_Photogeneration.jl",
    "Ex105_PSC_gradedFlux.jl",
    "Ex106_PSC_SurfaceRecombination.jl",
    "Ex107_MoS2_withIons_BarrierLowering.jl",
    "Ex108_CIGS.jl",
    "Ex109_PSC_NonDimensional.jl",
    "Ex201_PSC_tensorGrid.jl",
    "Ex202_Laser_simple.jl",
]


cleanexamples()

module_examples = @docmodules(exampledir, modules, use_module_titles = true)

#############################################################################

notebooks = [
    "PSC example" => "PSC_example.jl",
]

#pluto_examples = @docplutonotebooks(notebookdir, notebooks, iframe = false)

#############################################################################

makedocs(
    sitename = "ChargeTransport.jl",
    modules = [ChargeTransport],
    clean = false,
    doctest = true,
    authors = "D. Abdel, Z. Amer, P. Farrell, J. Fuhrmann, P. Jaap",
    repo = "https://github.com/WIAS-PDELib/ChargeTransport.jl",
    pages = [
        "Home" => "general.md",
        "Changelog" => "changes.md",
        "Mathematical drift-diffusion models" => "backgroundinfo.md",
        "How to get started" => [
            #"GeneralInformation.md",
            "GaAs.md",
            "PSC.md",
        ],
        "Types, Constructors and Methods" => "allindex.md",
        #"Pluto Notebooks" => pluto_examples,
        "Examples" => module_examples,
    ]
)

cleanexamples()

deploydocs(
    repo = "github.com/WIAS-PDELib/ChargeTransport.jl.git"
)
