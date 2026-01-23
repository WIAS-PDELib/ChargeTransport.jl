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
        typestr = gettoken(tks)

        expecttoken(tks, "dimension")
        expecttoken(tks, "=")
        dim = parse(Int, gettoken(tks))

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

        if typestr == "scalar"
            # read scalar functions
            func = zeros(nval)
            for i in 1:nval
                func[i] = parse(Float64, gettoken(tks))
            end
            data[funcname] = func
            @info "Parsed scalar function $(funcname) valid on $(validity)"

        else
            # skip vector function for now (not sure for what we need them)
            for _ in 1:(nval * dim)
                gettoken(tks)
            end
            @info "Skipped $(typestr) function $(funcname)"
        end

        expecttoken(tks, "}")

    end
    return data
end
