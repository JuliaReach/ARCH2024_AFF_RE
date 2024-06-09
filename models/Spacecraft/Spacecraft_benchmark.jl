using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "Rendezvous"
cases = ["NA01", "NA01-discrete",
         "A01", "A01-discrete",
         "A02", "A02-discrete",
         "A03", "A03-discrete",
         "A04", # "A04-discrete",
         #"A05", "A05-discrete",
         #"A06", "A06-discrete",
         #"A07", "A07-discrete",
         #"A08", "A08-discrete",
         "U01", # "U01-discrete",
         "U02", "U02-discrete"]

SUITE[model] = BenchmarkGroup()

include("Spacecraft.jl")
validation = Dict(c => 0 for c in cases)

LazySets.deactivate_assertions()

boxdirs = BoxDirections{Float64, Vector{Float64}}(5)

# ----------------------------------------
#  NA01 (dense time)
# ----------------------------------------
case = "NA01"
prob_NA01 = spacecraft(abort_time=-1.)

sol_NA01 = solve(prob_NA01, alg=BOX(δ=0.04),
                 clustering_method=LazyClustering(),
                 intersection_method=TemplateHullIntersection(boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))
property = SR02_specification(sol_NA01)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_NA01, alg=BOX(δ=0.04),
                 clustering_method=LazyClustering(),
                 intersection_method=TemplateHullIntersection($boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))

# ----------------------------------------
#  NA01 (discrete time)
# ----------------------------------------
case = "NA01-discrete"
sol_NA01 = solve(prob_NA01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                 clustering_method=LazyClustering(),
                 intersection_method=TemplateHullIntersection(boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))
property = SR02_specification(sol_NA01)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_NA01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                 clustering_method=LazyClustering(),
                 intersection_method=TemplateHullIntersection($boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))

 # ----------------------------------------
 #  A01 (dense time)
 # ----------------------------------------
case = "A01"
prob_A01 = spacecraft(abort_time=120.)

sol_A01 = solve(prob_A01, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(3),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A01)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_A01, alg=BOX(δ=0.04),
                 clustering_method=LazyClustering(3),
                 intersection_method=TemplateHullIntersection($boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))

# ----------------------------------------
#  A01 (discrete time)
# ----------------------------------------
case = "A01-discrete"
sol_A01 = solve(prob_A01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                clustering_method=LazyClustering(3),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A01)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_A01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                  clustering_method=LazyClustering(3),
                  intersection_method=TemplateHullIntersection($boxdirs),
                  intersect_source_invariant=false,
                  tspan = (0.0 .. 300.0))


# ----------------------------------------
#  A02 (dense time)
# ----------------------------------------
case = "A02"
prob_A02 = spacecraft(abort_time=[120., 125.])

sol_A02 = solve(prob_A02, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(3),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A02)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_A02, alg=BOX(δ=0.04),
               clustering_method=LazyClustering(3),
               intersection_method=TemplateHullIntersection($boxdirs),
               intersect_source_invariant=false,
               tspan = (0.0 .. 300.0))

# ----------------------------------------
#  A02 (discrete time)
# ----------------------------------------
case = "A02-discrete"
sol_A02 = solve(prob_A02, alg=BOX(δ=0.1, approx_model=NoBloating()),
              clustering_method=LazyClustering(3),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A02)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_A02, alg=BOX(δ=0.1, approx_model=NoBloating()),
               clustering_method=LazyClustering(3),
               intersection_method=TemplateHullIntersection($boxdirs),
               intersect_source_invariant=false,
               tspan = (0.0 .. 300.0))


# ----------------------------------------
#  A03 (dense time)
# ----------------------------------------
case = "A03"
prob_A03 = spacecraft(abort_time=[120., 145.])

sol_A03 = solve(prob_A03, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(4),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A03)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_A03, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(4),
              intersection_method=TemplateHullIntersection($boxdirs),
              intersect_source_invariant=false,
              intersect_source_invariant_method=TemplateHullIntersection($boxdirs),
              tspan = (0.0 .. 300.0))

# ----------------------------------------
#  A03 (discrete time)
# ----------------------------------------
case = "A03-discrete"
sol_A03 = solve(prob_A03, alg=BOX(δ=0.1, approx_model=NoBloating()),
             clustering_method=LazyClustering(4),
             intersection_method=TemplateHullIntersection(boxdirs),
             intersect_source_invariant=false,
             intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
             tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A03)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_A03, alg=BOX(δ=0.1, approx_model=NoBloating()),
              clustering_method=LazyClustering(4),
              intersection_method=TemplateHullIntersection($boxdirs),
              intersect_source_invariant=false,
              intersect_source_invariant_method=TemplateHullIntersection($boxdirs),
              tspan = (0.0 .. 300.0))


# ----------------------------------------
#  A04 (dense time)
# ----------------------------------------
case = "A04"
prob_A04 = spacecraft(abort_time=240.)

sol_A04 = solve(prob_A04, alg=BOX(δ=0.01),
              clustering_method=LazyClustering(40),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A04)
validation[case] = Int(property)
SUITE[model][case] = @benchmarkable solve($prob_A04, alg=BOX(δ=0.01),
            clustering_method=LazyClustering(40),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            intersect_source_invariant_method=TemplateHullIntersection($boxdirs),
            tspan = (0.0 .. 300.0))

# ----------------------------------------
#  A04 (discrete time)
# ----------------------------------------
#=
case = "A04-discrete"
sol_A04 = solve(prob_A04, alg=BOX(δ=0.1, approx_model=NoBloating()),
           clustering_method=LazyClustering(16),
           intersection_method=TemplateHullIntersection(boxdirs),
           intersect_source_invariant=false,
           intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
           tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A04)
validation[case] = Int(property)
SUITE[model][cases[10]] = @benchmarkable solve($prob_A04, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            intersect_source_invariant_method=TemplateHullIntersection($boxdirs),
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A05 (dense time)
# ----------------------------------------
#=
case = "A05"
prob_A05 = spacecraft(abort_time=[235., 240.])

sol_A05 = solve(prob_A05, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(13),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A05)
validation[case] = Int(property)
SUITE[model][cases[11]] = @benchmarkable solve($prob_A05, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(13),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A05 (discrete time)
# ----------------------------------------
#=
case = "A05-discrete"
sol_A05 = solve(prob_A05, alg=BOX(δ=0.1, approx_model=NoBloating()),
           clustering_method=LazyClustering(13),
           intersection_method=TemplateHullIntersection(boxdirs),
           intersect_source_invariant=false,
           tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A05)
validation[case] = Int(property)
SUITE[model][cases[12]] = @benchmarkable solve($prob_A05, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(13),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A06 (dense time)
# ----------------------------------------
#=
case = "A06"
prob_A06 = spacecraft(abort_time=[230., 240.])

sol_A06 = solve(prob_A06, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(16),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A06)
validation[case] = Int(property)
SUITE[model][cases[13]] = @benchmarkable solve($prob_A06, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A06 (discrete time)
# ----------------------------------------
#=
case = "A06-discrete"
sol_A06 = solve(prob_A06, alg=BOX(δ=0.1, approx_model=NoBloating()),
           clustering_method=LazyClustering(16),
           intersection_method=TemplateHullIntersection(boxdirs),
           intersect_source_invariant=false,
           tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A06)
validation[case] = Int(property)
SUITE[model][cases[14]] = @benchmarkable solve($prob_A06, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A07 (dense time)
# ----------------------------------------
#=
case = "A07"
prob_A07 = spacecraft(abort_time=[50., 150.])

sol_A07 = solve(prob_A07, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(50),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A07)
validation[case] = Int(property)
SUITE[model][cases[15]] = @benchmarkable solve($prob_A07, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(50),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A07 (discrete time)
# ----------------------------------------
#=
case = "A07-discrete"
sol_A07 = solve(prob_A07, alg=BOX(δ=0.1, approx_model=NoBloating()),
           clustering_method=LazyClustering(50),
           intersection_method=TemplateHullIntersection(boxdirs),
           intersect_source_invariant=false,
           tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A07)
validation[case] = Int(property)
SUITE[model][cases[16]] = @benchmarkable solve($prob_A07, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(50),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A08 (dense time)
# ----------------------------------------
#=
case = "A08"
prob_A08 = spacecraft(abort_time=[0., 240.])

sol_A08 = solve(prob_A08, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(800),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A08)
validation[case] = Int(property)
SUITE[model][cases[17]] = @benchmarkable solve($prob_A08, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(800),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A08 (discrete time)
# ----------------------------------------
#=
case = "A08-discrete"
sol_A08 = solve(prob_A08, alg=BOX(δ=0.1, approx_model=NoBloating()),
                clustering_method=LazyClustering(900),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A08)
validation[case] = Int(property)
SUITE[model][cases[18]] = @benchmarkable solve($prob_A08, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(900),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  U01 (dense time)
# ----------------------------------------
case = "U01"
prob_U01 = spacecraft(abort_time=260.)

sol_U01 = solve(prob_U01, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(16),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_U01)
validation[case] = Int(property)
SUITE[model][cases[10]] = @benchmarkable solve($prob_U01, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))

# ----------------------------------------
#  U01 (discrete time)
# ----------------------------------------
#=
case = "U01-discrete"
sol_U01 = solve(prob_U01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                clustering_method=LazyClustering(16),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_U01)
validation[case] = Int(property)
SUITE[model][cases[12]] = @benchmarkable solve($prob_U01, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  U02 (dense time)
# ----------------------------------------
case = "U02"
prob_U02 = spacecraft(abort_time=[0., 260.])

sol_U02 = solve(prob_U02, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(16),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_U02)
validation[case] = Int(property)
SUITE[model][cases[11]] = @benchmarkable solve($prob_U02, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))

# ----------------------------------------
#  U02 (discrete time)
# ----------------------------------------
case = "U02-discrete"
sol_U02 = solve(prob_U02, alg=BOX(δ=0.1, approx_model=NoBloating()),
                clustering_method=LazyClustering(16),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_U02)
validation[case] = Int(property)
SUITE[model][cases[12]] = @benchmarkable solve($prob_U02, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))



sol_NA01 = nothing
sol_A01 = nothing
sol_A02 = nothing
sol_A03 = nothing
sol_A04 = nothing
#sol_A05 = nothing
#sol_A06 = nothing
#sol_A07 = nothing
#sol_A08 = nothing
sol_U01 = nothing
sol_U02 = nothing
GC.gc()

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

function plot_SRNA01()
    prob_NA01 = spacecraft(abort_time=-1.)
    sol = solve(prob_NA01, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))

    idx_approaching = findall(x -> x == 1, location.(sol))
    idx_attempt = findall(x -> x == 2, location.(sol))
    idx_aborting = findall(x -> x == 3, location.(sol))

    fig = plot(legend=:bottomright, tickfont=font(30, "Times"), guidefontsize=45,
               xlab=L"x", ylab=L"y",
               xtick=[-800, -600, -400, -200.0, 0.0], ytick=[-400, -300, -200, -100, 0.],
               xlims=(-1000.0, 0.0), ylims=(-450.0, 0.0),
               bottom_margin=-6mm, left_margin=-3mm, right_margin=1mm, top_margin=3mm,
               size=(1000, 1000))

    for idx in idx_approaching
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:lightgreen, alpha=1.)
    end
    for idx in idx_attempt
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:red, alpha=1.)
    end
    for idx in idx_aborting
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:cyan, alpha=1.)
    end
    fig
end

fig = plot_SRNA01()
savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP24-JuliaReach-$model-SRNA01"))

# ----

function plot_SRA01()
    prob_A01 = spacecraft(abort_time=120.)
    @time sol = solve(prob_A01, alg=BOX(δ=0.04),
                      clustering_method=LazyClustering(3),
                      intersection_method=TemplateHullIntersection(boxdirs),
                      intersect_source_invariant=false,
                      tspan = (0.0 .. 300.0))

    idx_approaching = findall(x -> x == 1, location.(sol))
    idx_attempt = findall(x -> x == 2, location.(sol))
    idx_aborting = findall(x -> x == 3, location.(sol))

    fig = plot(legend=:bottomright, tickfont=font(30, "Times"), guidefontsize=45,
               xlab=L"x", ylab=L"y", lw=0.0,
               xtick=[-750, -500, -250, 0, 250.], ytick=[-400, -300, -200, -100, 0.],
               xlims=(-1000.0, 400.0), ylims=(-450.0, 0.0),
               bottom_margin=-6mm, left_margin=-3mm, right_margin=1mm, top_margin=3mm,
               size=(1000, 1000))

    for idx in idx_approaching
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:lightgreen, alpha=1.)
    end
    for idx in idx_attempt
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:red, alpha=1.)
    end
    for idx in idx_aborting
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:cyan, alpha=1.)
    end
    fig
end

fig = plot_SRA01()
savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP24-JuliaReach-$model-SRA01"))

GC.gc()
