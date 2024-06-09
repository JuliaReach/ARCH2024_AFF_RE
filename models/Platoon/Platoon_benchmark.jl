using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "Platoon"
cases = ["PLAD01-BND42", "PLAD01-BND42-discrete",
         "PLAD01-BND30", "PLAD01-BND30-discrete"]

SUITE[model] = BenchmarkGroup()

include("Platoon.jl")
validation = Dict(c => 0 for c in cases)

LazySets.deactivate_assertions()

boxdirs = BoxDirections{Float64, Vector{Float64}}(10)
octdirs = CustomDirections([Vector(vi) for vi in OctDirections(10)])

# ----------------------------------------
#  PLAD01-BND42 (dense time)
# ----------------------------------------
case = "PLAD01-BND42"
prob_PLAD01_BND42 = platoon(; deterministic_switching=true)
sol_PLAD01_BND42 = solve(prob_PLAD01_BND42, alg=BOX(δ=0.01),
                         clustering_method=BoxClustering(1),
                         intersection_method=TemplateHullIntersection(boxdirs),
                         intersect_source_invariant=false,
                         tspan = (0.0 .. 20.0))
property = dmin_specification(sol_PLAD01_BND42, -42.0)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_PLAD01_BND42, alg=BOX(δ=0.01),
                               clustering_method=BoxClustering(1),
                               intersection_method=TemplateHullIntersection($boxdirs),
                               intersect_source_invariant=false,
                               tspan = (0.0 .. 20.0))

# ----------------------------------------
#  PLAD01-BND42 (discrete time)
# ----------------------------------------
case = "PLAD01-BND42-discrete"
prob_PLAD01_BND42 = platoon(; deterministic_switching=true)
sol_PLAD01_BND42_d = solve(prob_PLAD01_BND42,
                           alg=BOX(δ=0.1, approx_model=NoBloating()),
                           clustering_method=BoxClustering(1, [3,1,1,1,1,1,1,1,1,1]),
                           intersection_method=TemplateHullIntersection(boxdirs),
                           intersect_source_invariant=false,
                           tspan = (0.0 .. 20.0))
property = dmin_specification(sol_PLAD01_BND42_d, -42.0)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_PLAD01_BND42,
                               alg=BOX(δ=0.1, approx_model=NoBloating()),
                               clustering_method=BoxClustering(1, [3,1,1,1,1,1,1,1,1,1]),
                               intersection_method=TemplateHullIntersection($boxdirs),
                               intersect_source_invariant=false,
                               tspan = (0.0 .. 20.0))

# ----------------------------------------
#  PLAD01-BND30 (dense time)
# ----------------------------------------
case = "PLAD01-BND30"
prob_PLAD01_BND30 = platoon(; deterministic_switching=true)
imethod = TemplateHullIntersection(octdirs)
cmethod = LazyClustering(1)
alg = LGG09(δ=0.03, template=octdirs, approx_model=Forward(setops=octdirs))
sol_PLAD01_BND30 = solve(prob_PLAD01_BND30,
                         alg=alg,
                         clustering_method=cmethod,
                         intersection_method=imethod,
                         intersect_source_invariant=false,
                         tspan = (0.0 .. 20.0))
property = dmin_specification(sol_PLAD01_BND30, -30.0)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_PLAD01_BND30,
                            alg=$alg,
                            clustering_method=$cmethod,
                            intersection_method=$imethod,
                            intersect_source_invariant=false,
                            tspan = (0.0 .. 20.0))

# ----------------------------------------
#  PLAD01-BND30 (discrete time)
# ----------------------------------------
case = "PLAD01-BND30-discrete"
prob_PLAD01_BND30 = platoon(; deterministic_switching=true)
imethod = TemplateHullIntersection(octdirs)
cmethod = LazyClustering(1)
alg = LGG09(δ=0.1, template=octdirs, approx_model=NoBloating())
sol_PLAD01_BND30_d = solve(prob_PLAD01_BND30,
                           alg=alg,
                           clustering_method=cmethod,
                           intersection_method=imethod,
                           intersect_source_invariant=false,
                           tspan = (0.0 .. 20.0))
property = dmin_specification(sol_PLAD01_BND30_d, -30.0)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_PLAD01_BND30,
                            alg=$alg,
                            clustering_method=$cmethod,
                            intersection_method=$imethod,
                            intersect_source_invariant=false,
                            tspan = (0.0 .. 20.0))

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

# ==============================================================================
# Plot
# ==============================================================================

if !@isdefined TARGET_FOLDER
    TARGET_FOLDER = @__DIR__
end

fig = Plots.plot()
Plots.plot!(fig, sol_PLAD01_BND30, vars=(0, 1), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
    tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"t",
    ylab=L"x_{1}",
    xtick=[0, 5, 10, 15, 20.], ytick=[-30, -20, -10, 0],
    xlims=(0., 20.), ylims=(-31, 7),
    bottom_margin=0mm, left_margin=-2mm, right_margin=4mm, top_margin=0mm,
    size=(1000, 1000))
Plots.hline!(fig, [-30.0], lc=:red, ls=:dash, lw=2, lab="")

savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP24-JuliaReach-$model-PLAD01-BND30"))

# release memory
sol_PLAD01_BND42 = nothing
sol_PLAD01_BND42_d = nothing
sol_PLAD01_BND30 = nothing
sol_PLAD01_BND30_d = nothing
GC.gc()
