"""
    read_diodat(filename)

Read scalar and vector vertex functions from  WIAS-TeSCA dios "*.dat" output files.
Return a dictionary containing the functions found.
"""
function read_diodat(fname)
    tks = TokenStream(fname)
    expecttoken(tks, "DF-ISE")
    expecttoken(tks, "text")
    data = Dict{String, Union{Vector{Float64}, Matrix{Float64}}}()
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

        function add_entry!(data, funcname, func)
            key = funcname
            if haskey(data, key)
                i = 2
                key = "$(funcname)_$i"
                while haskey(data, key)
                    i = i + 1
                    key = "$(funcname)_$i"
                end
            end
            data[key] = func
            return key
        end

        function parse_or_nan(tks)
            x = tryparse(Float64, tks)
            return x === nothing ? NaN : x
        end

        if typestr == "scalar"
            func = zeros(nval)
            for ival in 1:nval
                func[ival] = parse_or_nan(gettoken(tks))
            end
            key = add_entry!(data, funcname, func)
            @info """Scalar $(key) on $(validity) ∈ $(extrema(func))"""
        else
            func = zeros(dim, nval)
            for ival in 1:nval
                for idim in 1:dim
                    func[idim, ival] = parse_or_nan(gettoken(tks))
                end
            end
            key = add_entry!(data, funcname, func)
            @info """Vector $(key) on $(validity) ∈ $(extrema(func, dims = (2,)))"""
        end

        expecttoken(tks, "}")

    end
    return data
end
