using ReachabilityAnalysis
import JSON

const known_keys = ("version", "A", "B", "C", "U", "X0", "tend", "unsafeSet")
const compatible_version = "1.0"

function parse_problem(path)
    dict = JSON.parsefile(path)

    # validate input
    if length(dict) != 8 || !all(in(keys(dict)), known_keys)
        @warn "file $path has unknown structure; parsing may fail or result " *
              "in unexpected output"
    end
    if get(dict, "version", "-1") != compatible_version
        @warn "incompatible version; only $compatible_version is supported"
    end

    # parse IVP
    ivp = parse_ivp(dict)

    # parse unsafe states
    unsafe_states = parse_unsafe_states(dict)

    # parse C matrix
    C = parse_matrix(dict, "C")

    # build predicate for guaranteed safety
    predicate_safe_set(R) = isdisjoint(C * R, unsafe_states)

    function predicate_safe(sol; silent::Bool=false)
        for R in sol
            if !predicate_safe_set(R)
                silent || println("  Potential violation for time range $(tspan(R)).")
                return false
            end
        end
        return true
    end

    # build predicate for guaranteed violation
    if unsafe_states isa UnionSet || unsafe_states isa UnionSetArray
        predicate_unsafe_set = R -> any(C * R ⊆ U for U in unsafe_states)
    else
        predicate_unsafe_set = R -> C * R ⊆ unsafe_states
    end

    function predicate_unsafe(sol; silent::Bool=false)
        for R in sol
            if predicate_unsafe_set(R)
                silent || println("  Guaranteed violation for time range $(tspan(R)).")
                return true
            end
        end
        return false
    end

    # parse time horizon
    T = parse_T(dict)

    return ivp, (predicate_safe, predicate_unsafe, unsafe_states, C), T
end

function parse_ivp(dict)
    A = parse_matrix(dict, "A")
    B = parse_matrix(dict, "B")
    U = parse_set(dict, "U")
    X0 = parse_set(dict, "X0")

    n = size(A, 2)
    X = Universe(n)

    S = @system(x' = A * x + B * u, x ∈ X, u ∈ U)
    return @ivp(S, x(0) ∈ X0)
end

# convert a vector of `Any`-vectors (rows) to a matrix
function parse_matrix(dict, name)
    rows = dict[name]
    return Float64.(reduce(hcat, rows)')
end

function parse_set(dict, name)
    return parse_set(dict[name])
end

function parse_set(dict)
    type = get(dict, "type", nothing)
    if isnothing(type)
        throw(ArgumentError("expected `type` key in set description"))
    end
    if type == "interval"
        return parse_box(dict)
    end
    throw(ArgumentError("unsupported set type `$type`"))
end

function parse_box(dict)
    l = Float64.(dict["lowerbound"])
    u = Float64.(dict["upperbound"])
    return Hyperrectangle(low=l, high=u)
end

function parse_unsafe_states(dict)
    list = dict["unsafeSet"]
    sets = [parse_set(e) for e in list]
    if length(sets) == 1
        return sets[1]
    elseif length(sets) == 2
        return UnionSet(sets[1], sets[2])
    end
    return UnionSetArray(sets)
end

function parse_T(dict)
    return Float64(dict["tend"])
end
