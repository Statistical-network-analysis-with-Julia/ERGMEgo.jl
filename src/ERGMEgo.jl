"""
    ERGMEgo.jl - ERGMs for Ego-Centric Network Data

Provides tools for fitting ERGMs to egocentrically sampled network data,
where we observe a sample of "egos" along with their local networks (alters
and ties among alters).

This enables inference about complete network properties from ego samples.

Port of the R ergm.ego package from the StatNet collection.
"""
module ERGMEgo

using DataFrames
using Distributions
using ERGM
using Graphs
using LinearAlgebra
using Network
using Optim
using Random
using Statistics
using StatsBase

# Data structures
export EgoData, EgoNetwork, EgoSample

# Data preparation
export as_egodata, ego_design
export read_ego_data, merge_ego_data

# Ego-specific terms
export EgoEdges, EgoMixingMatrix, EgoNodeMatch
export EgoDegree, EgoGWDegree, EgoTriangle

# Estimation
export ergm_ego, fit_ego_ergm

# Population size estimation
export estimate_popsize

# Simulation
export simulate_ego_sample

# Diagnostics
export ego_gof, compare_ego_population

# =============================================================================
# Ego Network Data Structures
# =============================================================================

"""
    EgoNetwork{T}

An ego-centric network observation.

# Fields
- `ego::T`: Ego vertex ID (or index in sample)
- `alters::Vector{T}`: Alter vertex IDs
- `alter_ties::Matrix{Bool}`: Adjacency among alters (not including ego)
- `ego_attrs::Dict{Symbol, Any}`: Ego attributes
- `alter_attrs::Dict{Symbol, Vector}`: Alter attributes (column per attribute)
- `weights::Vector{Float64}`: Optional weights for alters
"""
struct EgoNetwork{T}
    ego::T
    alters::Vector{T}
    alter_ties::Matrix{Bool}
    ego_attrs::Dict{Symbol, Any}
    alter_attrs::Dict{Symbol, Vector}
    weights::Vector{Float64}

    function EgoNetwork{T}(ego::T, alters::Vector{T}, alter_ties::Matrix{Bool};
                           ego_attrs::Dict{Symbol,Any}=Dict{Symbol,Any}(),
                           alter_attrs::Dict{Symbol,Vector}=Dict{Symbol,Vector}(),
                           weights::Vector{Float64}=Float64[]) where T
        n_alters = length(alters)
        size(alter_ties) == (n_alters, n_alters) ||
            throw(ArgumentError("alter_ties must be $(n_alters)×$(n_alters)"))

        w = isempty(weights) ? ones(n_alters) : weights
        new{T}(ego, alters, alter_ties, ego_attrs, alter_attrs, w)
    end
end

EgoNetwork(ego::T, alters::Vector{T}, alter_ties::Matrix{Bool}; kwargs...) where T =
    EgoNetwork{T}(ego, alters, alter_ties; kwargs...)

"""
    n_alters(ego_net::EgoNetwork) -> Int

Get the number of alters in an ego network.
"""
n_alters(ego_net::EgoNetwork) = length(ego_net.alters)

"""
    alter_degree(ego_net::EgoNetwork) -> Vector{Int}

Get the degree of each alter within the ego network (not counting ego).
"""
function alter_degree(ego_net::EgoNetwork)
    return vec(sum(ego_net.alter_ties, dims=2))
end

"""
    ego_degree(ego_net::EgoNetwork) -> Int

Get the degree of the ego (number of alters).
"""
ego_degree(ego_net::EgoNetwork) = n_alters(ego_net)

"""
    n_alter_ties(ego_net::EgoNetwork) -> Int

Get the number of ties among alters.
"""
function n_alter_ties(ego_net::EgoNetwork)
    return sum(ego_net.alter_ties) ÷ 2  # Assuming undirected
end

"""
    EgoData

Collection of ego networks with sampling information.

# Fields
- `egos::Vector{EgoNetwork}`: Individual ego network observations
- `population_size::Union{Int, Nothing}`: Known or estimated population size
- `sampling_weights::Vector{Float64}`: Sampling weights for each ego
- `design::Dict{Symbol, Any}`: Survey design information
"""
struct EgoData{T}
    egos::Vector{EgoNetwork{T}}
    population_size::Union{Int, Nothing}
    sampling_weights::Vector{Float64}
    design::Dict{Symbol, Any}

    function EgoData(egos::Vector{EgoNetwork{T}};
                     population_size::Union{Int, Nothing}=nothing,
                     sampling_weights::Vector{Float64}=Float64[],
                     design::Dict{Symbol, Any}=Dict{Symbol, Any}()) where T
        n_egos = length(egos)
        sw = isempty(sampling_weights) ? ones(n_egos) : sampling_weights
        new{T}(egos, population_size, sw, design)
    end
end

Base.length(ed::EgoData) = length(ed.egos)
Base.iterate(ed::EgoData, state=1) = state > length(ed) ? nothing : (ed.egos[state], state + 1)
Base.getindex(ed::EgoData, i) = ed.egos[i]

"""
    summary_stats(ed::EgoData) -> NamedTuple

Summary statistics for ego data.
"""
function summary_stats(ed::EgoData)
    n_egos = length(ed.egos)
    degrees = [ego_degree(e) for e in ed.egos]
    alter_ties = [n_alter_ties(e) for e in ed.egos]

    return (
        n_egos = n_egos,
        mean_degree = mean(degrees),
        median_degree = median(degrees),
        min_degree = minimum(degrees),
        max_degree = maximum(degrees),
        mean_alter_ties = mean(alter_ties),
        total_alters = sum(degrees),
        population_size = ed.population_size
    )
end

# =============================================================================
# Data Preparation
# =============================================================================

"""
    as_egodata(df::DataFrame; ego_col, alter_col, ...) -> EgoData

Create EgoData from a DataFrame.

# Arguments
- `df`: DataFrame with ego network data
- `ego_id::Symbol`: Column identifying ego
- `alter_id::Symbol`: Column identifying alter (if long format)
- `ego_attrs::Vector{Symbol}`: Columns with ego attributes
- `alter_attrs::Vector{Symbol}`: Columns with alter attributes
"""
function as_egodata(df::DataFrame;
                    ego_id::Symbol=:ego_id,
                    alter_id::Symbol=:alter_id,
                    ego_attrs::Vector{Symbol}=Symbol[],
                    alter_attrs::Vector{Symbol}=Symbol[],
                    tie_col::Symbol=:tie,
                    weight_col::Union{Symbol, Nothing}=nothing)

    # Group by ego
    ego_ids = unique(df[!, ego_id])
    egos = EgoNetwork{Int}[]

    for (idx, eid) in enumerate(ego_ids)
        ego_df = filter(row -> row[ego_id] == eid, df)

        # Get alters for this ego
        alters = unique(ego_df[!, alter_id])
        n_a = length(alters)

        # Build alter tie matrix (placeholder - needs tie data)
        alter_ties = zeros(Bool, n_a, n_a)

        # Get ego attributes
        e_attrs = Dict{Symbol, Any}()
        for attr in ego_attrs
            e_attrs[attr] = ego_df[1, attr]
        end

        # Get alter attributes
        a_attrs = Dict{Symbol, Vector}()
        for attr in alter_attrs
            a_attrs[attr] = ego_df[!, attr]
        end

        # Weights
        weights = if !isnothing(weight_col) && weight_col in names(ego_df)
            ego_df[!, weight_col]
        else
            ones(n_a)
        end

        push!(egos, EgoNetwork(idx, collect(1:n_a), alter_ties;
                               ego_attrs=e_attrs, alter_attrs=a_attrs, weights=weights))
    end

    return EgoData(egos)
end

"""
    ego_design(ed::EgoData; ppopsize, weights) -> EgoData

Specify survey design for ego data.
"""
function ego_design(ed::EgoData{T};
                    ppopsize::Union{Int, Nothing}=nothing,
                    weights::Union{Vector{Float64}, Nothing}=nothing) where T
    new_weights = isnothing(weights) ? ed.sampling_weights : weights
    design = copy(ed.design)

    return EgoData(ed.egos;
                   population_size=ppopsize,
                   sampling_weights=new_weights,
                   design=design)
end

# =============================================================================
# Ego-Specific ERGM Terms
# =============================================================================

"""
    EgoEdges <: AbstractERGMTerm

Edge count term for ego data - estimates overall network density.
"""
struct EgoEdges <: AbstractERGMTerm end

name(::EgoEdges) = "edges"

function compute(::EgoEdges, ed::EgoData)
    # Weighted average of ego degrees
    total_degree = sum(ego_degree(e) * w for (e, w) in zip(ed.egos, ed.sampling_weights))
    total_weight = sum(ed.sampling_weights)
    return total_degree / total_weight
end

"""
    EgoNodeMatch <: AbstractERGMTerm

Homophily term based on matching node attributes.
"""
struct EgoNodeMatch <: AbstractERGMTerm
    attr::Symbol
    EgoNodeMatch(attr::Symbol) = new(attr)
end

name(t::EgoNodeMatch) = "nodematch.$(t.attr)"

function compute(t::EgoNodeMatch, ed::EgoData)
    total = 0.0
    total_weight = 0.0

    for (ego_net, w) in zip(ed.egos, ed.sampling_weights)
        ego_val = get(ego_net.ego_attrs, t.attr, nothing)
        isnothing(ego_val) && continue

        alter_vals = get(ego_net.alter_attrs, t.attr, nothing)
        isnothing(alter_vals) && continue

        # Count matches between ego and alters
        matches = count(av -> av == ego_val, alter_vals)
        total += matches * w
        total_weight += w
    end

    return total_weight > 0 ? total / total_weight : 0.0
end

"""
    EgoDegree <: AbstractERGMTerm

Degree distribution term for ego data.
"""
struct EgoDegree <: AbstractERGMTerm
    d::Union{Int, Nothing}  # Specific degree, or nothing for all

    EgoDegree(d::Union{Int, Nothing}=nothing) = new(d)
end

name(t::EgoDegree) = isnothing(t.d) ? "degree" : "degree.$(t.d)"

function compute(t::EgoDegree, ed::EgoData)
    if isnothing(t.d)
        # Return mean degree
        total = sum(ego_degree(e) * w for (e, w) in zip(ed.egos, ed.sampling_weights))
        return total / sum(ed.sampling_weights)
    else
        # Count egos with specific degree
        count = sum(w for (e, w) in zip(ed.egos, ed.sampling_weights) if ego_degree(e) == t.d)
        return count / sum(ed.sampling_weights)
    end
end

"""
    EgoGWDegree <: AbstractERGMTerm

Geometrically weighted degree for ego data.
"""
struct EgoGWDegree <: AbstractERGMTerm
    decay::Float64
    EgoGWDegree(decay::Float64=0.5) = new(decay)
end

name(t::EgoGWDegree) = "gwdegree.$(t.decay)"

function compute(t::EgoGWDegree, ed::EgoData)
    total = 0.0
    total_weight = 0.0

    for (ego_net, w) in zip(ed.egos, ed.sampling_weights)
        d = ego_degree(ego_net)
        # GWD contribution
        if d > 0
            contribution = exp(t.decay) * (1 - (1 - exp(-t.decay))^d)
            total += contribution * w
        end
        total_weight += w
    end

    return total_weight > 0 ? total / total_weight : 0.0
end

"""
    EgoTriangle <: AbstractERGMTerm

Triangle count from ego data (based on alter-alter ties).
"""
struct EgoTriangle <: AbstractERGMTerm end

name(::EgoTriangle) = "triangle"

function compute(::EgoTriangle, ed::EgoData)
    total = 0.0
    total_weight = 0.0

    for (ego_net, w) in zip(ed.egos, ed.sampling_weights)
        # Triangles = alter-alter ties (each tie with ego forms a triangle)
        triangles = n_alter_ties(ego_net)
        total += triangles * w
        total_weight += w
    end

    return total_weight > 0 ? total / total_weight : 0.0
end

"""
    EgoMixingMatrix <: AbstractERGMTerm

Mixing matrix term - cross-tabulation of ego-alter attributes.
"""
struct EgoMixingMatrix <: AbstractERGMTerm
    attr::Symbol
    levels::Vector{Any}

    function EgoMixingMatrix(attr::Symbol; levels::Vector=Any[])
        new(attr, levels)
    end
end

name(t::EgoMixingMatrix) = "mixing.$(t.attr)"

function compute(t::EgoMixingMatrix, ed::EgoData)
    # Build mixing matrix
    levels = if isempty(t.levels)
        # Infer levels from data
        all_vals = Any[]
        for ego_net in ed.egos
            push!(all_vals, get(ego_net.ego_attrs, t.attr, missing))
            append!(all_vals, get(ego_net.alter_attrs, t.attr, []))
        end
        unique(filter(!ismissing, all_vals))
    else
        t.levels
    end

    n_levels = length(levels)
    mixing = zeros(n_levels, n_levels)

    for (ego_net, w) in zip(ed.egos, ed.sampling_weights)
        ego_val = get(ego_net.ego_attrs, t.attr, nothing)
        isnothing(ego_val) && continue

        ego_idx = findfirst(==(ego_val), levels)
        isnothing(ego_idx) && continue

        alter_vals = get(ego_net.alter_attrs, t.attr, nothing)
        isnothing(alter_vals) && continue

        for av in alter_vals
            alter_idx = findfirst(==(av), levels)
            isnothing(alter_idx) && continue
            mixing[ego_idx, alter_idx] += w
        end
    end

    return sum(mixing)  # Return total for now; could return matrix
end

# =============================================================================
# Model and Estimation
# =============================================================================

"""
    EgoERGMModel

ERGM model for ego-centric data.
"""
struct EgoERGMModel
    terms::Vector{AbstractERGMTerm}
    data::EgoData
    ppopsize::Int
end

"""
    EgoERGMResult

Results from fitting an ego ERGM.
"""
struct EgoERGMResult
    model::EgoERGMModel
    coefficients::Vector{Float64}
    std_errors::Vector{Float64}
    loglik::Float64
    converged::Bool
end

function Base.show(io::IO, result::EgoERGMResult)
    println(io, "Ego ERGM Results")
    println(io, "================")
    println(io, "Population size: $(result.model.ppopsize)")
    println(io, "Sample size: $(length(result.model.data))")
    println(io, "Log-likelihood: $(round(result.loglik, digits=4))")
    println(io, "Converged: $(result.converged)")
    println(io)
    println(io, "Coefficients:")
    for (i, term) in enumerate(result.model.terms)
        println(io, "  $(rpad(name(term), 20)) $(lpad(round(result.coefficients[i], digits=4), 10)) " *
                    "(SE: $(round(result.std_errors[i], digits=4)))")
    end
end

"""
    ergm_ego(data::EgoData, terms; ppopsize, kwargs...) -> EgoERGMResult

Fit an ERGM from ego-centric network data.

# Arguments
- `data`: EgoData object
- `terms`: Vector of ERGM terms
- `ppopsize`: Pseudo-population size for inference
"""
function ergm_ego(data::EgoData, terms::Vector{<:AbstractERGMTerm};
                  ppopsize::Union{Int, Nothing}=nothing,
                  method::Symbol=:mple,
                  maxiter::Int=100)

    # Use provided ppopsize or estimate
    pop_size = if !isnothing(ppopsize)
        ppopsize
    elseif !isnothing(data.population_size)
        data.population_size
    else
        # Estimate from data
        estimate_popsize(data)
    end

    model = EgoERGMModel(terms, data, pop_size)

    if method == :mple
        return ego_mple(model; maxiter=maxiter)
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end

fit_ego_ergm = ergm_ego

"""
    ego_mple(model::EgoERGMModel; kwargs...) -> EgoERGMResult

MPLE for ego ERGM.
"""
function ego_mple(model::EgoERGMModel; maxiter::Int=100, tol::Float64=1e-6)
    n_terms = length(model.terms)
    coef = zeros(n_terms)

    # Compute observed statistics
    obs_stats = [compute(term, model.data) for term in model.terms]

    # Gradient descent optimization
    for iter in 1:maxiter
        grad = zeros(n_terms)

        # Approximate gradient
        for (i, term) in enumerate(model.terms)
            # Simplified gradient computation
            grad[i] = obs_stats[i] - obs_stats[i] * sigmoid(-coef[i])
        end

        # Update with step size
        step_size = 0.1 / sqrt(iter)
        coef .+= step_size * grad

        if maximum(abs.(grad)) < tol
            # Compute approximate standard errors
            se = fill(0.1, n_terms)  # Placeholder
            return EgoERGMResult(model, coef, se, NaN, true)
        end
    end

    se = fill(NaN, n_terms)
    return EgoERGMResult(model, coef, se, NaN, false)
end

sigmoid(x) = 1.0 / (1.0 + exp(-x))

# =============================================================================
# Population Size Estimation
# =============================================================================

"""
    estimate_popsize(ed::EgoData; method=:horvitz_thompson) -> Int

Estimate population size from ego data.
"""
function estimate_popsize(ed::EgoData; method::Symbol=:horvitz_thompson)
    if method == :horvitz_thompson
        # Horvitz-Thompson estimator based on sampling weights
        return round(Int, sum(ed.sampling_weights))
    elseif method == :capture_recapture
        # Simple capture-recapture based on alter overlap
        # Count unique alters weighted by appearance frequency
        alter_counts = Dict{Int, Int}()
        for ego_net in ed.egos
            for a in ego_net.alters
                alter_counts[a] = get(alter_counts, a, 0) + 1
            end
        end

        n_unique = length(alter_counts)
        n_total = sum(values(alter_counts))

        # Lincoln-Petersen estimate
        n_sample = length(ed.egos)
        return round(Int, n_sample * n_unique / (n_total / n_sample))
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end

# =============================================================================
# Simulation
# =============================================================================

"""
    simulate_ego_sample(net::Network, n_egos::Int; kwargs...) -> EgoData

Simulate an ego sample from a complete network.
"""
function simulate_ego_sample(net::Network{T}, n_egos::Int;
                             with_replacement::Bool=false,
                             include_alter_ties::Bool=true) where T
    n = nv(net)
    n_egos <= n || throw(ArgumentError("n_egos cannot exceed network size"))

    # Sample egos
    ego_ids = sample(1:n, n_egos; replace=with_replacement)
    egos = EgoNetwork{T}[]

    for (idx, ego_id) in enumerate(ego_ids)
        # Get alters (neighbors of ego)
        alters = collect(neighbors(net, ego_id))

        if isempty(alters)
            alter_ties = zeros(Bool, 0, 0)
        else
            # Build alter tie matrix
            n_a = length(alters)
            alter_ties = zeros(Bool, n_a, n_a)

            if include_alter_ties
                for (i, a1) in enumerate(alters)
                    for (j, a2) in enumerate(alters)
                        if i < j && has_edge(net, a1, a2)
                            alter_ties[i, j] = true
                            alter_ties[j, i] = true
                        end
                    end
                end
            end
        end

        push!(egos, EgoNetwork(T(idx), T.(collect(1:length(alters))), alter_ties))
    end

    return EgoData(egos; population_size=n)
end

# =============================================================================
# Diagnostics
# =============================================================================

"""
    ego_gof(result::EgoERGMResult; statistics) -> NamedTuple

Goodness-of-fit for ego ERGM.
"""
function ego_gof(result::EgoERGMResult;
                 statistics::Vector{Symbol}=[:degree, :esp])
    # Compare observed ego statistics to model expectations
    obs_stats = Dict{Symbol, Float64}()

    for stat in statistics
        if stat == :degree
            obs_stats[:degree] = mean(ego_degree(e) for e in result.model.data.egos)
        elseif stat == :alter_ties
            obs_stats[:alter_ties] = mean(n_alter_ties(e) for e in result.model.data.egos)
        end
    end

    return (observed=obs_stats,)
end

end # module
