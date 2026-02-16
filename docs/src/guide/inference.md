# Population Inference

ERGMEgo.jl estimates population-level network properties from ego-centric samples. This page covers the estimation procedure, population size estimation, survey design considerations, diagnostics, and best practices.

## Overview

The inference pipeline follows these steps:

1. **Specify survey design**: Set population size and sampling weights
2. **Construct pseudo-population**: Scale ego data to represent the target population
3. **Compute statistics**: Calculate ego-level sufficient statistics for each ERGM term
4. **Estimate parameters**: Maximize the pseudo-likelihood to obtain ERGM coefficients
5. **Validate**: Assess model fit and sensitivity

## Pseudo-Likelihood Estimation

### Why Pseudo-Likelihood?

For ego-centric data, we cannot compute the full ERGM likelihood because we do not observe the complete network. Instead, ERGMEgo.jl uses pseudo-likelihood:

- Each ego provides a local view of the network
- Ego-level statistics estimate population-level ERGM sufficient statistics
- Pseudo-likelihood combines these local estimates into a coherent model

This approach is computationally efficient and yields consistent estimates under standard survey sampling assumptions.

### How It Works

For each ego observation, the conditional distribution of ties given the local network structure is:

$$P(\text{tie}_{ij} = 1 \mid \text{rest}) = \text{logistic}\left(\sum_k \theta_k \, \delta_k(y, i, j)\right)$$

Where $\delta_k(y, i, j)$ is the change statistic for term $k$. ERGMEgo.jl estimates these change statistics from ego data and aggregates across the sample using sampling weights.

## Fitting Models

### Basic Usage

```julia
result = ergm_ego(data, terms; ppopsize=10000)
```

### Full Options

```julia
result = ergm_ego(data, terms;
    ppopsize = 10000,   # Pseudo-population size
    method = :mple,     # Estimation method
    maxiter = 100       # Maximum iterations
)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `ppopsize` | `Union{Int, Nothing}` | Pseudo-population size | From data or estimated |
| `method` | `Symbol` | Estimation method (currently `:mple`) | `:mple` |
| `maxiter` | `Int` | Maximum optimization iterations | `100` |

### Alternative Syntax

The `fit_ego_ergm` function is an alias for `ergm_ego`:

```julia
# These are equivalent
result = ergm_ego(data, terms; ppopsize=10000)
result = fit_ego_ergm(data, terms; ppopsize=10000)
```

### Population Size Priority

When multiple sources provide population size, ERGMEgo.jl uses this priority order:

1. `ppopsize` keyword argument (highest priority)
2. `data.population_size` field
3. Estimated from sampling weights via `estimate_popsize` (fallback)

## Population Size Estimation

When the true population size is unknown, ERGMEgo.jl provides two estimation methods.

### Horvitz-Thompson Estimator

The Horvitz-Thompson estimator uses sampling weights to estimate the population size:

```julia
N = estimate_popsize(data; method=:horvitz_thompson)
```

**How it works**: Sums the sampling weights across all egos. If each ego $i$ has sampling weight $w_i$ (typically the inverse of the sampling probability), then:

$$\hat{N} = \sum_{i=1}^{n} w_i$$

**When to use**: When reliable sampling weights are available (e.g., from a probability sample).

### Capture-Recapture Estimator

The capture-recapture estimator uses alter overlap between egos to estimate the population:

```julia
N = estimate_popsize(data; method=:capture_recapture)
```

**How it works**: Uses the Lincoln-Petersen method. If $n$ egos name a total of $n_{\text{unique}}$ unique alters with $n_{\text{total}}$ total alter mentions:

$$\hat{N} = \frac{n \times n_{\text{unique}}}{n_{\text{total}} / n}$$

**When to use**: When alter identifiers allow matching across egos (i.e., you can detect when two egos name the same alter).

### Choosing a Method

| Scenario | Recommended Method |
|----------|--------------------|
| Probability sample with known weights | Horvitz-Thompson |
| Alter IDs available and matchable | Capture-recapture |
| No weights, no alter matching | Use domain knowledge or census data |

### Using Estimated Population Size

```julia
# Estimate population size
N = estimate_popsize(data; method=:horvitz_thompson)

# Use in model fitting
result = ergm_ego(data, terms; ppopsize=N)

# Or set on data object
data = ego_design(data; ppopsize=N)
result = ergm_ego(data, terms)
```

## Understanding Results

The `EgoERGMResult` object contains:

| Field | Type | Description |
|-------|------|-------------|
| `model` | `EgoERGMModel` | The fitted model (terms, data, population size) |
| `coefficients` | `Vector{Float64}` | Estimated coefficients |
| `std_errors` | `Vector{Float64}` | Standard errors |
| `loglik` | `Float64` | Log-pseudo-likelihood at convergence |
| `converged` | `Bool` | Whether optimization converged |

### Displaying Results

```julia
println(result)
```

Output:

```text
Ego ERGM Results
================
Population size: 10000
Sample size: 50
Log-likelihood: -234.5678
Converged: true

Coefficients:
  edges                    -2.1234 (SE: 0.0812)
  triangle                  0.3456 (SE: 0.0923)
  nodematch.gender          0.5678 (SE: 0.0567)
```

## Interpreting Coefficients

Coefficients have the same interpretation as standard ERGM coefficients -- they describe the log-odds effect on tie formation in the population network:

| Term | Coefficient | exp(coef) | Interpretation |
|------|-------------|-----------|----------------|
| EgoEdges | -2.0 | 0.14 | Baseline tie probability ~12% |
| EgoTriangle | 0.5 | 1.65 | Each shared partner increases tie odds by 65% |
| EgoNodeMatch(:race) | 0.8 | 2.23 | Same-race ties are 2.2x more likely |
| EgoGWDegree(0.5) | -0.3 | 0.74 | Slight preference for uniform degree |

### Confidence Intervals

```julia
using Distributions

alpha = 0.05
z = quantile(Normal(), 1 - alpha/2)

lower = result.coefficients .- z .* result.std_errors
upper = result.coefficients .+ z .* result.std_errors

for (i, term) in enumerate(result.model.terms)
    println("$(name(term)): [$(round(lower[i], digits=3)), $(round(upper[i], digits=3))]")
end
```

## Sensitivity Analysis

### Sensitivity to Population Size

The pseudo-population size affects estimates. Check sensitivity:

```julia
pop_sizes = [1000, 5000, 10000, 50000, 100000]

for N in pop_sizes
    result = ergm_ego(data, terms; ppopsize=N)
    println("N=$N: ", round.(result.coefficients, digits=3))
end
```

If results are stable across reasonable population sizes, the estimates are robust. The edges coefficient is most sensitive to population size, while other coefficients (homophily, clustering) are typically more stable.

### Sensitivity to Sampling Weights

```julia
# Equal weights
data_equal = ego_design(data; weights=ones(length(data)))
result_equal = ergm_ego(data_equal, terms; ppopsize=10000)

# Original weights
result_weighted = ergm_ego(data, terms; ppopsize=10000)

println("Equal weights: ", round.(result_equal.coefficients, digits=3))
println("Survey weights: ", round.(result_weighted.coefficients, digits=3))
```

## Model Comparison

### Comparing Models

```julia
# Model 1: Density only
result1 = ergm_ego(data, [EgoEdges()]; ppopsize=10000)

# Model 2: Add clustering
result2 = ergm_ego(data, [EgoEdges(), EgoTriangle()]; ppopsize=10000)

# Model 3: Add homophily
result3 = ergm_ego(data, [EgoEdges(), EgoTriangle(), EgoNodeMatch(:race)]; ppopsize=10000)

for (i, r) in enumerate([result1, result2, result3])
    println("Model $i: $(length(r.coefficients)) terms, converged=$(r.converged)")
end
```

## Goodness-of-Fit

Use `ego_gof` to assess model fit by comparing observed ego statistics to model expectations:

```julia
gof = ego_gof(result; statistics=[:degree, :alter_ties])

println("Observed mean degree: ", gof.observed[:degree])
println("Observed mean alter ties: ", gof.observed[:alter_ties])
```

## Convergence

### Checking Convergence

```julia
if result.converged
    println("Model converged")
else
    println("WARNING: Model did not converge")
end
```

### Common Issues

| Issue | Symptom | Solution |
|-------|---------|----------|
| Too few egos | Non-convergence | Collect more ego observations |
| Sparse alter-alter ties | Triangle term unstable | Remove EgoTriangle term |
| Rare attribute levels | Large standard errors | Collapse attribute categories |
| Too many terms | Slow convergence | Reduce model complexity |
| Incorrect ppopsize | Implausible coefficients | Verify population size |

### Handling Non-Convergence

```julia
# Increase iterations
result = ergm_ego(data, terms; ppopsize=10000, maxiter=500)

# Simplify model
terms_simple = [EgoEdges(), EgoNodeMatch(:race)]
result_simple = ergm_ego(data, terms_simple; ppopsize=10000)
```

## Simulation for Validation

Simulate ego samples from a complete network and verify that the method recovers known parameters:

```julia
using Network

# Create a network with known properties
net = Network{Int}(; n=100)
for i in 1:100, j in (i+1):100
    if rand() < 0.1
        add_edge!(net, i, j)
    end
end

# Draw ego sample
ego_data = simulate_ego_sample(net, 30;
    with_replacement=false,
    include_alter_ties=true
)

# Fit model
terms = [EgoEdges(), EgoTriangle()]
result = ergm_ego(ego_data, terms; ppopsize=100)

println("Estimated coefficients: ", round.(result.coefficients, digits=3))
```

## Best Practices

1. **Know your population size**: Results are sensitive to `ppopsize`; use the best estimate available
2. **Check convergence**: Always verify `result.converged == true`
3. **Start simple**: Begin with EgoEdges and add terms incrementally
4. **Use sampling weights**: If the ego sample is not a simple random sample
5. **Report population size**: Always report the `ppopsize` used in your analysis
6. **Run sensitivity analysis**: Check how results change with different population size estimates
7. **Validate with simulation**: Simulate ego samples from known networks to verify the method
8. **Sufficient sample**: More ego observations yield more reliable estimates
