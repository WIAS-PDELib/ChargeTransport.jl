using Test, VoronoiFVM
using Aqua
using ChargeTransport

# Mitigate https://github.com/JuliaLang/julia/issues/58634
function treshape(X::AbstractArray, n, m)
    Y = reshape(X, n, m)
    Y .= 1.0
    return X
end

function talloc(; n = 10, m = 20)
    X = rand(n, m)
    treshape(X, n, m)
    return @allocated treshape(X, n, m)
end

if talloc() > 0
    VoronoiFVM.check_allocs!(false)
    @warn "Disabling allocation checks due to julia issue #58634"
end

modname(fname) = splitext(basename(fname))[1]

#
# Include all Julia files in `testdir` whose name starts with `prefix`,
# Each file `prefixModName.jl` must contain a module named
# `prefixModName` which has a method test() returning true
# or false depending on success.
#
function run_tests_from_directory(testdir, prefix)
    println("Directory $(testdir):")
    examples = modname.(readdir(testdir))
    for example in examples
        if length(example) >= length(prefix) && example[1:length(prefix)] == prefix
            println("  $(example):")
            path = joinpath(testdir, "$(example).jl")
            @eval begin
                include($path)
                @test eval(Meta.parse("$($example).test()"))
            end
        end
    end
    return
end

function run_all_tests()
    return @time begin
        # @testset "Basictest" begin
        #     run_tests_from_directory(@__DIR__,"test_")
        # end
        @testset "Examples" begin
            run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Ex")
        end
        @testset "PSC" begin
            run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "PSC")
        end
    end
end

@testset "Aqua.jl" begin
    Aqua.test_all(ChargeTransport)
end

run_all_tests()
