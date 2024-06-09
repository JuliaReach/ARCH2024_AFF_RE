# deactivate plot GUI, which is not available in Docker
ENV["GKSwstype"] = "100"

# instantiate project
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

const TEST_LONG = true   # if true, the longer test suite is run; may take > 1 hour
                         # and requires at least 16gb RAM
global const TARGET_FOLDER = "results"
const RESULTS_FILE = "results.csv"

function main()
    if !isdir(TARGET_FOLDER)
        mkdir(TARGET_FOLDER)
    end
    global io = open(joinpath(TARGET_FOLDER, RESULTS_FILE), "w")
    print(io, "benchmark,instance,result,time\n")

    println("Running AFF benchmarks...")

    # Heat 3D benchmark
    println("###\nRunning Heat 3D benchmark\n###")
    include("models/Heat3D/heat3d_benchmark.jl")

    # Clamped Beam benchmark
    println("###\nRunning Clamped Beam benchmark\n###")
    include("models/Clamped/clamped_benchmark.jl")

    # Spacecraft benchmark
    if TEST_LONG
        println("###\nRunning Spacecraft benchmark\n###")
        include("models/Spacecraft/Spacecraft_benchmark.jl")
    end

    # Powertrain benchmark
    println("###\nRunning Powertrain benchmark\n###")
    include("models/Powertrain/Powertrain_benchmark.jl")

    # Platoon benchmark
    println("###\nRunning Platoon benchmark\n###")
    include("models/Platoon/Platoon_benchmark.jl")

    # Gearbox benchmark
    println("###\nRunning Gearbox benchmark\n###")
    include("models/Gearbox/gearbox_benchmark.jl")

    # Brake benchmark
    if TEST_LONG
        println("###\nRunning Electromechanic Brake benchmark\n###")
        include("models/EMBrake/embrake_benchmark.jl")
    end

    # Random benchmark
    println("###\nRunning Random benchmark\n###")
    try
        include("models/Random/random_benchmark.jl")
    catch e
        showerror(stdout, e)
        println()
    end

    print(io, "\n")
    println("Finished running benchmarks.")
    close(io)
    nothing
end

main()
