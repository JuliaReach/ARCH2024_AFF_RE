using ReachabilityAnalysis

const δ0 = 1e-2
const δ_min = 1e-5

function solve_random(ivp, specification, T; silent::Bool=true)
    predicate_safe, predicate_unsafe, unsafe_states, C = specification

    # simple check whether initial states are unsafe
    if !isdisjoint(C * ivp.x0, unsafe_states)
        silent || @info "proven unsafe in the initial states"
        return false
    end

    δ = δ0
    while δ >= δ_min
        silent || @info "δ = $δ"
        alg = GLGM06(; δ=δ, approx_model=Forward())
        sol = solve(ivp, alg; T=T);

        # check predicates
        if predicate_safe(sol; silent=silent)
            silent || @info "proven safe"
            return true
        elseif predicate_unsafe(sol; silent=silent)
            silent || @info "proven unsafe"
            return false
        end

        # refinement
        δ /= 2
    end
    silent || @info "giving up"
    return -1
end
