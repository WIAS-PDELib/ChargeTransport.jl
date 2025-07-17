using Test

# Activate assembly loop allocation checking
# as default.
ENV["VORONOIFVM_CHECK_ALLOCS"] = "false"

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
        if length(example) >= length(prefix) &&example[1:length(prefix)] == prefix
            println("  $(example):")
            path = joinpath(testdir, "$(example).jl")
            @eval begin
                include($path)
                # Compile + run test
                print("   compile:")
                @time @test eval(Meta.parse("$($example).test()"))
                # Second run: pure execution time.
                print("       run:")
                @time eval(Meta.parse("$($example).test()"))
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
    end
end

run_all_tests()
