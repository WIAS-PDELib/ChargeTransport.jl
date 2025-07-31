using Documenter, ExampleJuggler, Literate, ChargeTransport, PlutoSliderServer, VoronoiFVM

DocMeta.setdocmeta!(ExampleJuggler, :DocTestSetup, :(using ExampleJuggler); recursive = true)

exampledir = joinpath(@__DIR__, "..", "examples")


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

notebookjl = last.(notebooks)
notebookmd = []

# Use sliderserver to generate html
notebook_html_dir = joinpath(@__DIR__, "src", "nbhtml")
export_directory(
    joinpath(@__DIR__, "..", "pluto-examples"),
    notebook_paths = notebookjl,
    Export_output_dir = joinpath(notebook_html_dir),
    Export_offer_binder = false,
)
# generate frame markdown for each notebook
for notebook in notebookjl
    base = split(notebook, ".")[1]
    mdstring = """
    ##### [$(base).jl](@id $(base))
    [Download](https://github.com/WIAS-PDELib/ChargeTransport.jl/blob/master/pluto-examples/$(notebook))
    this [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook.
    ```@raw html
    <iframe style="height:20000px" width="100%" src="../$(base).html"> </iframe>
    ```
    """
    mdname = base * ".md"
    push!(notebookmd, joinpath("nbhtml", mdname))
    io = open(joinpath(notebook_html_dir, mdname), "w")
    write(io, mdstring)
    close(io)
end

notebooks = first.(notebooks) .=> notebookmd

#############################################################################

makedocs(
    sitename = "ChargeTransport.jl",
    modules = [ChargeTransport],
    clean = false,
    doctest = false,
    authors = "D. Abdel, Z. Amer, P. Farrell, J. Fuhrmann, P. Jaap",
    repo = "https://github.com/WIAS-PDELib/ChargeTransport.jl",
    pages = Any[
        "ChargeTransport.jl -- Simulating charge transport in semiconductors" => "general.md",
        "Changelog" => "changes.md",
        "Mathematical drift-diffusion models" => "backgroundinfo.md",
        "How to get started" => [
            "GeneralInformation.md",
            "GaAs.md",
            "PSC.md",
        ],
        "Types, Constructors and Methods" => "allindex.md",
        #"Pluto Notebooks" => notebooks,
        "Examples" => module_examples,
    ]
)

cleanexamples()

deploydocs(
    repo = "github.com/WIAS-PDELib/ChargeTransport.jl.git"
)
