"""
    read_diodat(filename)

Read scalar vertex functions from  WIAS-TeSCA dios "*.dat" output files.
Return a dictionary containing the functions found.
"""
function read_diodat(fname)
    tks = TokenStream(fname)
    expecttoken(tks, "DF-ISE")
    expecttoken(tks, "text")
    data = Dict{String, Vector{Float64}}()
    while !eof(tks)
        token = gettoken(tks)
        while token != "function" && !eof(tks)
            token = gettoken(tks)
        end
        if eof(tks)
            continue
        end
        expecttoken(tks, "=")
        funcname = gettoken(tks)

        expecttoken(tks, "type")
        expecttoken(tks, "=")
        expecttoken(tks, "scalar")

        expecttoken(tks, "dimension")
        expecttoken(tks, "=")
        expecttoken(tks, "1")

        expecttoken(tks, "location")
        expecttoken(tks, "=")
        expecttoken(tks, "vertex")

        expecttoken(tks, "validity")
        expecttoken(tks, "=")
        expecttoken(tks, "[")
        validity = gettoken(tks)
        expecttoken(tks, "]")

        expecttoken(tks, "Values")
        expecttoken(tks, "(")
        nval = parse(Int, gettoken(tks))
        expecttoken(tks, ")")

        expecttoken(tks, "{")
        func = zeros(nval)
        for i in 1:nval
            func[i] = parse(Float64, gettoken(tks))
        end
        expecttoken(tks, "}")
        @info "Parsed function $(funcname) valid on $(validity)"
        data[funcname] = func
    end
    return data
end
