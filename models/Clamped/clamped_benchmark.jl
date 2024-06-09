using Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets

# We don't use BenchmarkTools for this model

model = "Beam"
cases = ["CBC01", # constant input set
         "CBC01-discrete",
         "CBF01", # bounded but arbitrarily-varying input set
         "CBF01-discrete",
         "CBC02",
         "CBC02-discrete",
         "CBF02",
         "CBF02-discrete",
         "CBC03",
#          "CBC03-discrete",
#          "CBF03",
#          "CBF03-discrete"
        ]
if TEST_LONG
    push!(cases, "CBC03-discrete")
    push!(cases, "CBF03-discrete")
end

results = Dict(model => Dict(c => -1.0 for c in cases))
validation = Dict(c => 1 for c in cases)  # nothing to verify
sol = Dict{String,Any}()
# measure = [] # maximum value of the velocity at node 70

include("clamped.jl")
include("clamped_krylov.jl")

LazySets.deactivate_assertions()

# ----------------------------------------
#  CBC01
# ----------------------------------------

case = "CBC01"
m, N, Z = parse_clamped_beam("C", 100)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-6, vars=[70, 170], n=201)
res = @timed solve_krylov_discr(prob, NSTEPS=10_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CBC01-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[70, 170], n=201, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CBF01
# ----------------------------------------

case = "CBF01"
m, N, Z = parse_clamped_beam("F", 100)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-7, vars=[70, 170], n=200)
res = @timed solve(prob, NSTEPS=100_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CBF01-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[70, 170], n=200, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CBC02
# ----------------------------------------

case = "CBC02"
m, N, Z = parse_clamped_beam("C", 500)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-6, vars=[350, 850], n=1001)
res = @timed solve_krylov_discr(prob, NSTEPS=10_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CBC02-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[350, 850], n=1001, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CBF02
# ----------------------------------------

case = "CBF02"
m, N, Z = parse_clamped_beam("F", 500)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-7, vars=[350, 850], n=1000)
res = @timed solve(prob, NSTEPS=100_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CBF02-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[350, 850], n=1000, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CBC03
# ----------------------------------------

case = "CBC03"
m, N, Z = parse_clamped_beam("C", 1000)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-6, vars=[700, 1700], n=2001)
res = @timed solve_krylov_discr(prob, NSTEPS=10_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CBC03-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[700, 1700], n=2001, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CBF03
# ----------------------------------------

case = "CBF03"
m, N, Z = parse_clamped_beam("F", 1000)
#=
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-7, vars=[700, 1700], n=2000)
res = @timed solve(prob, NSTEPS=100_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time
=#

case = "CBF03-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[700, 1700], n=2000, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

if !@isdefined io
    io = stdout
end

for c in cases
   local t = results[model][c] # in seconds
   runtime = round(t, digits=4)
   print(io, "$model,$c,$(validation[c]),$(runtime)\n")
end

# ==============================================================================
# Plot
# ==============================================================================

if !@isdefined TARGET_FOLDER
    TARGET_FOLDER = @__DIR__
end

# Constant force case
fig = Plots.plot()
Plots.plot!(fig, sol["CBC01"], vars=(0, 170),
           color=:blue, alpha=0.5, lw=1.0, linecolor=:blue,
           tickfont=font(30, "Times"), guidefontsize=45,
           xlab=L"t",
           ylab=L"v_{70}",
           xtick=([0, 0.0025, 0.005, 0.0075, 0.01],
                  [L"0.0", L"0.0025", L"0.0050", L"0.0075", L"0.0100"]),
           ytick=([-80, -60, -40, -20, 0, 20, 40, 60, 80],
                  [L"-80", L"-60", L"-40", L"-20", L"0", L"20", L"40", L"60", L"80"]),
           xlims=(0, 0.01), ylims=(-80, 80),
           bottom_margin=0mm, left_margin=0mm, right_margin=15mm, top_margin=3mm,
           size=(1000, 1000))

# savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP24-JuliaReach-$model-CBC01.pdf"))
savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP24-JuliaReach-$model-CBC01.png"))

# Time-varying force case
fig = Plots.plot()
Plots.plot!(fig, sol["CBF01"], vars=(0, 170),
           color=:blue, alpha=0.5, lw=1.0, linecolor=:blue,
           tickfont=font(30, "Times"), guidefontsize=45,
           xlab=L"t",
           ylab=L"v_{70}",
           xtick=([0, 0.0025, 0.005, 0.0075, 0.01],
                  [L"0.0", L"0.0025", L"0.0050", L"0.0075", L"0.0100"]),
           ytick=([-80, -60, -40, -20, 0, 20, 40, 60, 80],
                  [L"-80", L"-60", L"-40", L"-20", L"0", L"20", L"40", L"60", L"80"]),
           xlims=(0, 0.01), ylims=(-80, 80),
           bottom_margin=0mm, left_margin=0mm, right_margin=15mm, top_margin=3mm,
           size=(1000, 1000))

# savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP24-JuliaReach-$model-CBF01.pdf"))
savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP24-JuliaReach-$model-CBF01.png"))

sol = nothing
GC.gc()
