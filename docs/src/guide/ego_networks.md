# Ego Networks and Data

This guide covers how to work with ego-centric network data in ERGMEgo.jl.

## Ego Networks

An `EgoNetwork` represents a single ego-centric observation: one focal actor (ego) along with their direct contacts (alters) and the ties among those alters.

```julia
using ERGMEgo

# Basic ego network: ego_id, alter_ids, alter-alter tie matrix
ego = EgoNetwork(1, [1, 2, 3],
                 [false true false;
                  true false true;
                  false true false])
```

### EgoNetwork Fields

| Field | Type | Description | Default |
|-------|------|-------------|---------|
| `ego` | `T` | Ego identifier | Required |
| `alters` | `Vector{T}` | Alter identifiers | Required |
| `alter_ties` | `Matrix{Bool}` | Adjacency among alters | Required |
| `ego_attrs` | `Dict{Symbol, Any}` | Ego attributes | Empty |
| `alter_attrs` | `Dict{Symbol, Vector}` | Alter attributes | Empty |
| `weights` | `Vector{Float64}` | Alter weights | Ones |

### Accessing Ego Network Data

```julia
ego = EgoNetwork(1, [1, 2, 3],
                 [false true false; true false true; false true false];
                 ego_attrs=Dict(:gender => "M", :age => 35),
                 alter_attrs=Dict(:gender => ["M", "F", "M"]))

ego.ego          # 1
ego.alters       # [1, 2, 3]
ego.alter_ties   # 3×3 Bool matrix
ego.ego_attrs    # Dict(:gender => "M", :age => 35)
ego.alter_attrs  # Dict(:gender => ["M", "F", "M"])
```

### Query Functions

```julia
n_alters(ego)      # 3 — number of alters
ego_degree(ego)    # 3 — same as n_alters (ego degree)
n_alter_ties(ego)  # 2 — number of undirected ties among alters
alter_degree(ego)  # [1, 2, 1] — degree of each alter within the ego net
```

### Creating Ego Networks with Attributes

```julia
# Categorical attributes
ego = EgoNetwork(1, [1, 2, 3, 4],
                 zeros(Bool, 4, 4);
                 ego_attrs=Dict(:race => "White", :education => "College"),
                 alter_attrs=Dict(
                     :race => ["White", "Black", "White", "Hispanic"],
                     :education => ["College", "HS", "College", "College"]
                 ))

# Numeric attributes
ego = EgoNetwork(1, [1, 2],
                 [false true; true false];
                 ego_attrs=Dict(:income => 55000.0),
                 alter_attrs=Dict(:income => [48000.0, 62000.0]))
```

### Alter-Alter Tie Matrix

The `alter_ties` matrix represents ties among alters (not including the ego). For undirected networks, it should be symmetric:

```julia
# 3 alters: alter 1↔2 connected, alter 2↔3 connected
ties = [false true  false;
        true  false true;
        false true  false]

ego = EgoNetwork(1, [1, 2, 3], ties)
n_alter_ties(ego)  # 2
```

For ego networks with no alter-alter tie information:

```julia
# No ties observed
ego = EgoNetwork(1, [1, 2, 3], zeros(Bool, 3, 3))
n_alter_ties(ego)  # 0
```

For isolated egos (no alters):

```julia
ego = EgoNetwork(1, Int[], zeros(Bool, 0, 0))
n_alters(ego)  # 0
```

## Ego Data

An `EgoData` object collects multiple ego networks with sampling metadata:

```julia
ego1 = EgoNetwork(1, [1, 2, 3],
                  [false true false; true false true; false true false])
ego2 = EgoNetwork(2, [1, 2],
                  [false true; true false])
ego3 = EgoNetwork(3, [1, 2, 3, 4],
                  zeros(Bool, 4, 4))

data = EgoData([ego1, ego2, ego3]; population_size=5000)
```

### EgoData Fields

| Field | Type | Description | Default |
|-------|------|-------------|---------|
| `egos` | `Vector{EgoNetwork{T}}` | Ego network observations | Required |
| `population_size` | `Union{Int, Nothing}` | Known/estimated population size | `nothing` |
| `sampling_weights` | `Vector{Float64}` | Sampling weights per ego | Ones |
| `design` | `Dict{Symbol, Any}` | Survey design info | Empty |

### Accessing EgoData

```julia
# Basic access
data[1]           # First ego network
length(data)      # Number of egos

# Iteration
for ego_net in data
    println("Ego $(ego_net.ego): $(n_alters(ego_net)) alters")
end
```

### Summary Statistics

```julia
stats = summary_stats(data)

stats.n_egos          # Number of ego observations
stats.mean_degree     # Mean ego degree
stats.median_degree   # Median ego degree
stats.min_degree      # Minimum ego degree
stats.max_degree      # Maximum ego degree
stats.mean_alter_ties # Mean alter-alter ties per ego
stats.total_alters    # Total alters across all egos
stats.population_size # Population size (if set)
```

## Loading Data

### From DataFrame

The most common way to create ego data from survey responses:

```julia
using DataFrames

df = DataFrame(
    ego_id = [1, 1, 1, 2, 2, 3, 3, 3, 3],
    alter_id = [1, 2, 3, 1, 2, 1, 2, 3, 4],
    gender_ego = ["M","M","M", "F","F", "M","M","M","M"],
    gender_alter = ["M","F","M", "F","M", "F","M","F","M"],
)

data = as_egodata(df;
    ego_id=:ego_id,
    alter_id=:alter_id,
    ego_attrs=[:gender_ego],
    alter_attrs=[:gender_alter]
)
```

### Custom Column Names

When your DataFrame has different column names:

```julia
df = DataFrame(
    respondent = [1, 1, 2, 2, 2],
    contact = [10, 20, 10, 30, 40],
    resp_race = ["W", "W", "B", "B", "B"],
    contact_race = ["W", "B", "B", "B", "W"]
)

data = as_egodata(df;
    ego_id=:respondent,
    alter_id=:contact,
    ego_attrs=[:resp_race],
    alter_attrs=[:contact_race]
)
```

### Adding Sampling Weights

When egos have unequal sampling probabilities (e.g., from a complex survey):

```julia
data = as_egodata(df; ego_id=:ego_id, alter_id=:alter_id)

# Add weights and population size
data = ego_design(data;
    ppopsize=250_000_000,
    weights=[1.5, 0.8, 1.2]
)
```

## Survey Design

### Specifying Population Size

The population size is important for scaling ego statistics to the population level:

```julia
# Known population size
data = EgoData(egos; population_size=10000)

# Or add later via ego_design
data = ego_design(data; ppopsize=10000)
```

### Sampling Weights

When egos are not sampled with equal probability:

```julia
# Inverse probability weights
weights = [1/p for p in sampling_probabilities]

data = ego_design(data;
    ppopsize=N,
    weights=weights
)
```

### Design Effects

The `design` dictionary stores additional survey design information:

```julia
data = EgoData(egos;
    population_size=10000,
    design=Dict{Symbol,Any}(
        :sampling_frame => "GSS 2018",
        :strata => :region
    ))
```

## Data Validation

### Checking Alter Tie Matrix

The alter tie matrix must be square with dimensions matching the number of alters:

```julia
# This is correct: 3 alters, 3×3 matrix
ego = EgoNetwork(1, [1, 2, 3], zeros(Bool, 3, 3))

# This throws an error: dimensions don't match
try
    ego = EgoNetwork(1, [1, 2, 3], zeros(Bool, 2, 2))
catch e
    println(e)  # ArgumentError: alter_ties must be 3×3
end
```

### Common Data Issues

| Issue | Symptom | Solution |
|-------|---------|----------|
| Asymmetric alter ties | Different results from i→j and j→i | Symmetrize: `ties = ties .| ties'` |
| Missing alter-alter ties | Triangles cannot be estimated | Use only degree-based terms |
| Egos with zero alters | May cause division issues | Filter or handle gracefully |
| Duplicate alters | Inflated statistics | De-duplicate before creating ego network |

### Filtering Egos

```julia
# Keep only egos with at least 2 alters
filtered_egos = [e for e in data.egos if n_alters(e) >= 2]
data_filtered = EgoData(filtered_egos; population_size=data.population_size)
```

## Simulating Ego Samples

For testing and validation, simulate ego samples from a complete network:

```julia
using Network

# Create a known network
net = Network(50; bipartite=false)
for i in 1:50, j in (i+1):50
    if rand() < 0.15
        add_edge!(net, i, j)
    end
end

# Draw ego sample
ego_data = simulate_ego_sample(net, 20;
    with_replacement=false,
    include_alter_ties=true
)

println("Sampled $(length(ego_data)) egos from network of $(nv(net)) nodes")
```

### Simulation Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `with_replacement` | Sample egos with replacement | `false` |
| `include_alter_ties` | Observe ties among alters | `true` |

Setting `include_alter_ties=false` simulates a scenario where only ego-alter ties are observed (no alter-alter information).

## Working with Multiple Surveys

To combine ego data from different sources:

```julia
# Two separate surveys
data1 = EgoData(egos_survey1; population_size=5000)
data2 = EgoData(egos_survey2; population_size=5000)

# Combine manually
combined_egos = vcat(data1.egos, data2.egos)
combined_weights = vcat(data1.sampling_weights, data2.sampling_weights)
combined = EgoData(combined_egos;
    population_size=5000,
    sampling_weights=combined_weights
)
```
