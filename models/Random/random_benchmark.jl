using BenchmarkTools
using BenchmarkTools: minimum, median

directory = joinpath(@__DIR__, "..", "..", "rand")
if !isdir(directory)
    throw(ArgumentError("no `rand` folder was found in the base directory"))
end

include("random_parser.jl")
include("random_algorithm.jl")

SUITE = BenchmarkGroup()
model = "Random"
cases = String[]
SUITE[model] = BenchmarkGroup()
validation = Dict{String,Int}()

LazySets.deactivate_assertions()
const silent = false

for (i, file) in enumerate(readdir(directory))
    if endswith(file, ".json")
        case, _ = split(file, ".")
        push!(cases, case)
        silent || @info "processing random instance $case"

        path = abspath(joinpath(directory, file))

        ivp, spec, T = parse_problem(path)

        res = solve_random(ivp, spec, T; silent=silent)
        validation[case] = Int(res)
        SUITE[model][case] = @benchmarkable solve_random($ivp, $predicate_safe,
            $predicate_unsafe, $T)
    else
        @warn "ignoring $file in $directory"
    end
end

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# tune parameters
tune!(SUITE)

# run the benchmarks
results = run(SUITE, verbose=true)

# return the sample with the smallest time value in each test
println("minimum time for each benchmark:\n", minimum(results))

# return the median for each test
println("median time for each benchmark:\n", median(results))

if !@isdefined io
    io = stdout
end

for c in cases
    local t = median(results[model][c]).time * 1e-9
    runtime = round(t, digits=4)
    print(io, "$model,$c,$(validation[c]),$(runtime)\n")
end
