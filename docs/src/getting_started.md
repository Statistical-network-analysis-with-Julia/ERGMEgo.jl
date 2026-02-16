# Getting Started

This tutorial walks through common use cases for ERGMEgo.jl, from creating ego network data to fitting and interpreting ego-centric ERGMs.

## Installation

Install ERGMEgo.jl from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/ERGMEgo.jl")
```

## Basic Workflow

The typical ERGMEgo.jl workflow consists of four steps:

1. **Prepare ego data** - Create or load ego-centric network observations
2. **Define terms** - Choose which ego-specific ERGM terms to include
3. **Fit the model** - Estimate parameters via pseudo-likelihood
4. **Interpret results** - Analyze the fitted model

## Step 1: Create Ego Network Data

Each ego network observation captures one actor's local network:

```julia
using ERGMEgo

# Create an ego network: ego_id, alter_ids, alter-alter tie matrix
ego1 = EgoNetwork(1, [1, 2, 3],
                  [false true false;
                   true false true;
                   false true false])

# Query basic properties
n_alters(ego1)       # 3
ego_degree(ego1)     # 3 (same as n_alters)
n_alter_ties(ego1)   # 2 (undirected ties among alters)
alter_degree(ego1)   # [1, 2, 1] (alter degrees within ego net)
```

### Adding Attributes

```julia
ego1 = EgoNetwork(1, [1, 2, 3],
                  [false true false; true false true; false true false];
                  ego_attrs=Dict(:gender => "M", :age => 35),
                  alter_attrs=Dict(
                      :gender => ["M", "F", "M"],
                      :age => [30, 28, 42]
                  ))
```

### Combining into EgoData

```julia
ego2 = EgoNetwork(2, [1, 2],
                  [false true; true false];
                  ego_attrs=Dict(:gender => "F", :age => 28),
                  alter_attrs=Dict(
                      :gender => ["F", "M"],
                      :age => [25, 32]
                  ))

ego3 = EgoNetwork(3, [1, 2, 3, 4],
                  [false true false false;
                   true false true true;
                   false true false false;
                   false true false false];
                  ego_attrs=Dict(:gender => "M", :age => 45),
                  alter_attrs=Dict(
                      :gender => ["M", "M", "F", "F"],
                      :age => [40, 50, 38, 47]
                  ))

# Combine with population size
data = EgoData([ego1, ego2, ego3]; population_size=10000)

# Summary
stats = summary_stats(data)
println("N egos: ", stats.n_egos)           # 3
println("Mean degree: ", stats.mean_degree) # 3.0
```

### Loading from a DataFrame

More commonly, you'll load ego data from existing survey data:

```julia
using DataFrames

df = DataFrame(
    ego_id = [1, 1, 1, 2, 2],
    alter_id = [1, 2, 3, 1, 2],
    gender_ego = ["M", "M", "M", "F", "F"],
    gender_alter = ["M", "F", "M", "F", "M"],
    age_ego = [35, 35, 35, 28, 28],
    age_alter = [30, 28, 42, 25, 32]
)

data = as_egodata(df;
    ego_id=:ego_id,
    alter_id=:alter_id,
    ego_attrs=[:gender_ego],
    alter_attrs=[:gender_alter, :age_alter]
)
```

### Specifying Survey Design

```julia
# Add population size and sampling weights
data = ego_design(data;
    ppopsize=250_000_000,
    weights=[1.5, 0.8, 1.2]
)
```

## Step 2: Define Terms

Ego-specific terms estimate population-level ERGM statistics from local data:

```julia
terms = [
    EgoEdges(),              # Network density
    EgoTriangle(),           # Clustering
    EgoNodeMatch(:gender),   # Gender homophily
]
```

### Exploring Available Terms

ERGMEgo.jl provides terms organized by type:

| Category | Terms | Description |
|----------|-------|-------------|
| **Structural** | `EgoEdges`, `EgoTriangle` | Basic network structure |
| **Degree** | `EgoDegree`, `EgoGWDegree` | Degree distribution |
| **Attribute** | `EgoNodeMatch`, `EgoMixingMatrix` | Homophily and mixing |

### Example: Comprehensive Model

```julia
terms = [
    # Structural effects
    EgoEdges(),
    EgoTriangle(),

    # Degree effects
    EgoGWDegree(0.5),

    # Attribute effects
    EgoNodeMatch(:gender),
    EgoNodeMatch(:age),
]
```

## Step 3: Fit the Model

Use `ergm_ego` to estimate the model:

```julia
result = ergm_ego(data, terms; ppopsize=10000)
```

### Key Parameters

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| `ppopsize` | Pseudo-population size | Known or estimated N |
| `method` | Estimation method | `:mple` (default) |
| `maxiter` | Maximum iterations | 100 |

### Using Estimated Population Size

If the population size is not known, it can be estimated:

```julia
# Estimate from sampling weights
N = estimate_popsize(data; method=:horvitz_thompson)

# Use in model
result = ergm_ego(data, terms; ppopsize=N)
```

## Step 4: Interpret Results

The result object contains coefficient estimates:

```julia
# Print formatted summary
println(result)

# Output:
# Ego ERGM Results
# ================
# Population size: 10000
# Sample size: 3
# Log-likelihood: -12.3456
# Converged: true
#
# Coefficients:
#   edges                    -2.1234 (SE: 0.1000)
#   triangle                  0.3456 (SE: 0.1000)
#   nodematch.gender          0.5678 (SE: 0.1000)
```

### Interpreting Coefficients

Coefficients have the same interpretation as standard ERGM coefficients:

| Coefficient | Interpretation |
|-------------|----------------|
| `edges` < 0 | Network is sparse (typical) |
| `triangle` > 0 | Clustering / transitivity present |
| `nodematch.gender` > 0 | Gender homophily (like connects with like) |

**Example interpretations:**

- `edges = -2.0` → Baseline tie probability is low (sparse network)
- `triangle = 0.5` → Strong tendency for friends-of-friends to be friends
- `nodematch.gender = 0.8` → Same-gender ties are exp(0.8) ≈ 2.2× more likely

## Complete Example

```julia
using ERGMEgo

# Create a small ego sample mimicking a social survey
egos = [
    EgoNetwork(1, [1, 2, 3],
               [false true false; true false true; false true false];
               ego_attrs=Dict(:race => "W"),
               alter_attrs=Dict(:race => ["W", "W", "B"])),

    EgoNetwork(2, [1, 2],
               [false false; false false];
               ego_attrs=Dict(:race => "B"),
               alter_attrs=Dict(:race => ["B", "B"])),

    EgoNetwork(3, [1, 2, 3, 4],
               [false true false false; true false false true;
                false false false false; false true false false];
               ego_attrs=Dict(:race => "W"),
               alter_attrs=Dict(:race => ["W", "B", "W", "W"])),

    EgoNetwork(4, [1],
               reshape([false], 1, 1);
               ego_attrs=Dict(:race => "B"),
               alter_attrs=Dict(:race => ["W"])),
]

data = EgoData(egos; population_size=5000)

# Define model
terms = [
    EgoEdges(),
    EgoTriangle(),
    EgoNodeMatch(:race),
]

# Fit
result = ergm_ego(data, terms; ppopsize=5000)

# View results
println(result)

# Check convergence
if result.converged
    println("\nModel converged successfully")
else
    println("\nWarning: Model did not converge")
end
```

## Simulating Ego Samples

For testing or validation, simulate ego samples from a known network:

```julia
using Network

# Create a complete network
net = Network(20; bipartite=false)
for i in 1:20, j in (i+1):20
    if rand() < 0.3
        add_edge!(net, i, j)
    end
end

# Sample 10 ego networks
ego_data = simulate_ego_sample(net, 10;
    with_replacement=false,
    include_alter_ties=true
)

# Fit model to simulated data
terms = [EgoEdges(), EgoTriangle()]
result = ergm_ego(ego_data, terms)
println(result)
```

## Best Practices

1. **Know your population size**: The pseudo-population size strongly affects results
2. **Check convergence**: Always verify `result.converged == true`
3. **Start simple**: Begin with `EgoEdges` before adding structural terms
4. **Use sampling weights**: If the sample is not a simple random sample
5. **Sufficient egos**: More ego observations yield more reliable estimates
6. **Include alter-alter ties**: These are essential for estimating triangle/clustering terms

## Next Steps

- Learn about [Ego Networks](guide/ego_networks.md) data structures in detail
- Explore all [Ego Terms](guide/terms.md) available
- Understand [Population Inference](guide/inference.md) methods
