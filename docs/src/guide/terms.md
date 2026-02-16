# Ego-Specific Terms

Terms in ERGMEgo.jl estimate population-level ERGM statistics from ego-centric data. All terms subtype `AbstractERGMTerm` and can be combined freely in models.

## Terms Interface

All terms implement two methods:

```julia
compute(term, ed::EgoData) -> Float64
name(term) -> String
```

The `compute` function calculates the weighted aggregate statistic across all ego networks in the `EgoData`.

## Term Categories

ERGMEgo.jl organizes terms into three categories:

| Type | Description | Examples |
|------|-------------|----------|
| Structural | Basic network structure | EgoEdges, EgoTriangle |
| Degree | Degree distribution | EgoDegree, EgoGWDegree |
| Attribute | Homophily and mixing | EgoNodeMatch, EgoMixingMatrix |

## Structural Terms

These capture basic structural properties of the network estimated from ego data.

### EgoEdges

Estimates overall network density from mean ego degree:

```julia
EgoEdges()
```

**Computation**: Weighted average of ego degrees across all ego networks.

**Interpretation**: The edges coefficient represents the baseline log-odds of a tie. A negative coefficient (typical) indicates a sparse network.

### EgoTriangle

Estimates the transitivity (clustering) of the network from alter-alter ties:

```julia
EgoTriangle()
```

**Computation**: Weighted average of alter-alter tie counts. Each alter-alter tie forms a triangle with the ego.

**Interpretation**: A positive coefficient indicates clustering — friends of friends tend to be friends.

**Requirement**: Alter-alter tie data must be available. If `alter_ties` is all zeros, this term contributes nothing.

Visual representation:

```text
    Alter₁
   / |
Ego  |     ← The alter₁—alter₂ tie completes the triangle
   \ |
    Alter₂
```

## Degree Terms

These capture heterogeneity in the degree distribution.

### EgoDegree

Degree distribution term:

```julia
# Mean ego degree
EgoDegree()

# Proportion of egos with specific degree d
EgoDegree(3)   # Proportion with exactly 3 alters
```

**Computation**:

- `EgoDegree()`: Weighted mean of ego degrees
- `EgoDegree(d)`: Weighted proportion of egos with degree `d`

**Interpretation**:

- `EgoDegree()` captures average connectivity
- `EgoDegree(d)` captures whether degree `d` is over- or under-represented

### EgoGWDegree

Geometrically weighted degree — a smooth summary of the degree distribution:

```julia
# With decay parameter
EgoGWDegree(0.5)   # decay = 0.5
EgoGWDegree(1.0)   # decay = 1.0
EgoGWDegree()      # decay = 0.5 (default)
```

**Parameter decay (α)**: Controls how quickly additional alters are down-weighted. Lower α = stronger down-weighting of high-degree egos.

**Computation**: For each ego with degree $d$:

$$\text{GWD contribution} = e^\alpha \left(1 - (1 - e^{-\alpha})^d\right)$$

**Interpretation**: A positive coefficient indicates a preference for degree heterogeneity. Unlike `EgoDegree(d)`, GWDegree captures the entire degree distribution in a single parameter.

### Choosing Between Degree Terms

| Scenario | Recommended Term |
|----------|-----------------|
| Simple model, one degree parameter | `EgoDegree()` |
| Testing for specific degree concentration | `EgoDegree(d)` |
| Smooth degree distribution modeling | `EgoGWDegree(α)` |
| Exploratory analysis | Start with `EgoDegree()`, refine with `EgoGWDegree` |

## Attribute Terms

These incorporate actor-level attributes for homophily and mixing effects.

### EgoNodeMatch

Homophily term — do egos tend to connect with alters who share the same attribute?

```julia
# Homophily on a categorical attribute
EgoNodeMatch(:gender)
EgoNodeMatch(:race)
EgoNodeMatch(:education)
```

**Computation**: Weighted average of the number of ego-alter attribute matches across all egos.

**Interpretation**: A positive coefficient indicates homophily — actors prefer ties with similar others.

**Example**:

```julia
ego = EgoNetwork(1, [1, 2, 3],
                 zeros(Bool, 3, 3);
                 ego_attrs=Dict(:gender => "M"),
                 alter_attrs=Dict(:gender => ["M", "F", "M"]))

# EgoNodeMatch(:gender) for this ego:
# Ego is "M", alters are ["M", "F", "M"]
# Matches: alter 1 (M=M) and alter 3 (M=M) → 2 matches out of 3
```

### EgoMixingMatrix

Full mixing pattern — captures who connects with whom based on attributes:

```julia
# Cross-tabulation of ego-alter attributes
EgoMixingMatrix(:race)
EgoMixingMatrix(:race; levels=["White", "Black", "Hispanic"])
```

**Computation**: Builds a weighted cross-tabulation of ego-alter attribute combinations.

**Parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `attr` | Attribute symbol | Required |
| `levels` | Specific levels to include | Auto-detected |

**Interpretation**: Captures asymmetric mixing patterns that `EgoNodeMatch` cannot. For example, whether White egos disproportionately nominate Black alters (or vice versa).

**Use case**: When you need to distinguish between different types of cross-group ties.

## Using Terms in Practice

### Building a Model

```julia
# Minimal model
terms_basic = [
    EgoEdges(),
]

# Standard model with clustering
terms_standard = [
    EgoEdges(),
    EgoTriangle(),
]

# Comprehensive model with attributes
terms_full = [
    EgoEdges(),
    EgoTriangle(),
    EgoGWDegree(0.5),
    EgoNodeMatch(:race),
    EgoNodeMatch(:education),
]
```

### Computing Terms Manually

```julia
# Compute a single term statistic
edges = EgoEdges()
value = compute(edges, data)
println("Mean ego degree: ", value)

# Compute all terms
for term in terms
    val = compute(term, data)
    println("$(name(term)): $val")
end
```

### Custom Term Names

All terms have automatic names based on their type:

```julia
name(EgoEdges())              # "edges"
name(EgoTriangle())           # "triangle"
name(EgoDegree())             # "degree"
name(EgoDegree(3))            # "degree.3"
name(EgoGWDegree(0.5))        # "gwdegree.0.5"
name(EgoNodeMatch(:gender))   # "nodematch.gender"
name(EgoMixingMatrix(:race))  # "mixing.race"
```

## Effect of Sampling Weights

All terms use sampling weights when computing statistics. Egos with higher weights contribute more to the aggregate statistic:

```julia
# Equal weights (default)
data_equal = EgoData(egos)
val_equal = compute(EgoEdges(), data_equal)

# Unequal weights
data_weighted = ego_design(data_equal;
    weights=[2.0, 1.0, 0.5]
)
val_weighted = compute(EgoEdges(), data_weighted)

# Results differ because egos contribute proportionally to their weights
```

## Choosing Terms

### By Research Question

| Question | Terms |
|----------|-------|
| What is the overall network density? | `EgoEdges` |
| Is the network clustered? | `EgoTriangle` |
| Is the degree distribution heterogeneous? | `EgoGWDegree`, `EgoDegree` |
| Is there racial homophily? | `EgoNodeMatch(:race)` |
| Are mixing patterns asymmetric? | `EgoMixingMatrix(:race)` |
| What is the degree distribution shape? | Multiple `EgoDegree(d)` terms |

### Best Practices

1. **Always include `EgoEdges`**: This controls for baseline density
2. **Use `EgoTriangle` only with alter-alter data**: Without alter-alter ties, this term is uninformative
3. **Prefer `EgoGWDegree` over multiple `EgoDegree(d)`**: One parameter captures the entire distribution
4. **Check attribute coverage**: `EgoNodeMatch` requires the attribute to be present in both ego and alter data
5. **Start simple**: Add terms incrementally and check convergence at each step
6. **Avoid redundant terms**: `EgoEdges` and `EgoDegree()` are closely related — typically include only one

## Relationship to Standard ERGM Terms

| Ego Term | Standard ERGM Term | What It Estimates |
|----------|--------------------|-------------------|
| `EgoEdges` | `Edges` | Network density |
| `EgoTriangle` | `Triangle` | Transitivity |
| `EgoDegree(d)` | `Degree(d)` | Degree-d proportion |
| `EgoGWDegree(α)` | `GWDegree(α)` | Degree distribution shape |
| `EgoNodeMatch(:x)` | `NodeMatch(:x)` | Attribute homophily |
| `EgoMixingMatrix(:x)` | `NodeMix(:x)` | Mixing patterns |
