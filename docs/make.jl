using Documenter, Literate, ChargeTransport, PlutoSliderServer

# Turn block comments starting in the first column into "normal" hash comments
# as block comments currently are not handled by Literate.jl.
# from https://github.com/chmerdon/GradientRobustMultiPhysics.jl/.
function hashify_block_comments(input)
    lines_in = collect(eachline(IOBuffer(input)))
    lines_out = IOBuffer()
    line_number = 0
    in_block_comment_region = false
    for line in lines_in
        line_number += 1
        if occursin(r"^#=", line)
            if in_block_comment_region
                error("line $(line_number): already in block comment region\n$(line)")
            end
            println(lines_out, replace(line, r"^#=" => "#"))
            in_block_comment_region = true
        elseif occursin(r"^=#", line)
            if !in_block_comment_region
                error("line $(line_number): not in block comment region\n$(line)")
            end
            println(lines_out, replace(line, r"^=#" => "#"))
            in_block_comment_region = false
        else
            if in_block_comment_region
                println(lines_out, "# " * line)
            else
                println(lines_out, line)
            end
        end
    end
    return String(take!(lines_out))
end

#
# Replace SOURCE_URL marker with github url of source
#
function replace_source_url(input, source_url)
    lines_in = collect(eachline(IOBuffer(input)))
    lines_out = IOBuffer()
    for line in lines_in
        println(lines_out, replace(line, "SOURCE_URL" => source_url))
    end
    return String(take!(lines_out))
end

#
# Generate Markdown pages from examples
#
example_jl_dir = joinpath(@__DIR__, "..", "examples")
example_md_dir = joinpath(@__DIR__, "src", "examples")

for example_source in readdir(example_jl_dir)
    base, ext = splitext(example_source)
    if ext == ".jl"
        source_url = "https://github.com/PatricioFarrell/ChargeTransport.jl/tree/master/examples/" * example_source
        preprocess(buffer) = replace_source_url(buffer, source_url) |> hashify_block_comments
        Literate.markdown(
            joinpath(@__DIR__, "..", "examples", example_source),
            example_md_dir,
            documenter = false,
            execute = false,
            info = false,
            preprocess = preprocess
        )
    end
end

generated_examples = joinpath.("examples", readdir(example_md_dir))

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
    [Download](https://github.com/PatricioFarrell/ChargeTransport.jl/blob/master/pluto-examples/$(notebook))
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
    clean = true,
    doctest = false,
    authors = "D. Abdel, P. Farrell, J. Fuhrmann",
    repo = "https://github.com/PatricioFarrell/ChargeTransport.jl",
    pages = Any[
        "ChargeTransport.jl -- Simulating charge transport in semiconductors" => "general.md",
        "Mathematical drift-diffusion models" => "backgroundinfo.md",
        "How to get started" => [
            "GeneralInformation.md",
            "GaAs.md",
            "PSC.md",
        ],
        "Types, Constructors and Methods" => "allindex.md",
        "Pluto Notebooks" => notebooks,
        "Examples" => generated_examples,
    ]
)

deploydocs(
    repo = "github.com/PatricioFarrell/ChargeTransport.jl.git"
)
