# ERGMEgo.jl

*ERGMs for Ego-Centric Network Data in Julia*

A Julia package for fitting Exponential Random Graph Models to egocentrically sampled network data.

## Overview

Ego-centric network sampling is one of the most common methods for collecting network data in social science. Rather than observing a complete network, a sample of "egos" is drawn from the population, and each ego reports on their local network: their direct contacts (alters) and the ties among those alters.

ERGMEgo.jl provides tools for fitting ERGMs to such ego-centric data, enabling inference about population-level network properties from ego samples.

ERGMEgo.jl is a port of [ergm.ego](https://github.com/statnet/ergm.ego) from the StatNet collection.

### What is an Ego-Centric Network?

An ego-centric network is a local view of a network centered on a focal actor (the "ego"):

```text
       Alter₁
      / |
Ego --- Alter₂
      \ |
       Alter₃
```

The data collected for each ego includes:

- **Ego-alter ties**: Which alters the ego is connected to
- **Alter-alter ties**: Which alters are connected to each other
- **Attributes**: Characteristics of the ego and their alters

### Key Concepts

| Concept | Description |
|---------|-------------|
| **Ego** | A sampled actor from the population |
| **Alter** | A direct contact of the ego |
| **Ego Network** | The ego, their alters, and ties among alters |
| **Ego Data** | A collection of ego networks with sampling information |
| **Pseudo-Population** | A synthetic population used for ERGM inference |
| **Sampling Weights** | Weights accounting for unequal sampling probabilities |

### Applications

Ego-centric ERGMs are widely used in:

- **Survey research**: Analyzing General Social Survey (GSS) network modules
- **Public health**: Studying disease transmission networks from contact surveys
- **Organizational studies**: Understanding professional networks from employee surveys
- **Community research**: Measuring social capital from ego-centric data
- **Hard-to-reach populations**: Network analysis when full network observation is infeasible

## Features

- **Ego network data structures**: `EgoNetwork` and `EgoData` types for managing ego-centric observations
- **Ego-specific ERGM terms**: Adapted terms that estimate population-level network statistics from ego data
- **Pseudo-likelihood estimation**: MPLE estimation adapted for ego-centric data
- **Population size estimation**: Horvitz-Thompson and capture-recapture methods
- **Simulation**: Generate ego samples from complete networks for testing and validation
- **Survey design support**: Incorporate sampling weights and population size information

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/ERGMEgo.jl")
```

Or for development:

```julia
using Pkg
Pkg.develop(path="/path/to/ERGMEgo.jl")
```

## Quick Start

```julia
using ERGMEgo

# Create ego network observations
ego1 = EgoNetwork(1, [1, 2, 3],
                  [false true false; true false true; false true false];
                  ego_attrs=Dict(:gender => "M"),
                  alter_attrs=Dict(:gender => ["M", "F", "M"]))

ego2 = EgoNetwork(2, [1, 2],
                  [false true; true false];
                  ego_attrs=Dict(:gender => "F"),
                  alter_attrs=Dict(:gender => ["F", "M"]))

# Combine into EgoData with population size
data = EgoData([ego1, ego2]; population_size=1000)

# Define ego-specific ERGM terms
terms = [
    EgoEdges(),              # Network density
    EgoTriangle(),           # Clustering (transitivity)
    EgoNodeMatch(:gender),   # Gender homophily
]

# Fit ego ERGM
result = ergm_ego(data, terms; ppopsize=1000)

# View results
println(result)
```

## Choosing Terms

| Use Case | Recommended Terms |
|----------|------------------|
| Network density | [`EgoEdges`](@ref) |
| Degree heterogeneity | [`EgoDegree`](@ref), [`EgoGWDegree`](@ref) |
| Transitivity / clustering | [`EgoTriangle`](@ref) |
| Homophily (categorical) | [`EgoNodeMatch`](@ref) |
| Mixing patterns | [`EgoMixingMatrix`](@ref) |

## Documentation

```@contents
Pages = [
    "getting_started.md",
    "guide/ego_networks.md",
    "guide/terms.md",
    "guide/inference.md",
    "api/types.md",
    "api/terms.md",
    "api/estimation.md",
]
Depth = 2
```

## Theoretical Background

### Ego-Centric ERGMs

ERGMs model the probability of a network as:

$$P(Y = y) = \frac{1}{\kappa(\theta)} \exp\left(\sum_k \theta_k g_k(y)\right)$$

Where $g_k(y)$ are sufficient statistics and $\theta_k$ are parameters. The key challenge with ego-centric data is that we observe only local neighborhoods rather than the full network.

### Pseudo-Population Approach

ERGMEgo.jl uses the pseudo-population approach from Krivitsky & Morris (2017):

1. **Construct a pseudo-population** from the ego sample, scaled to the estimated population size
2. **Compute ego-level sufficient statistics** that estimate the population-level ERGM statistics
3. **Maximize the pseudo-likelihood** over the pseudo-population to obtain parameter estimates

This approach yields consistent estimates under standard survey sampling assumptions.

### Estimable Statistics

Not all ERGM statistics can be estimated from ego-centric data. Statistics that depend only on local structure are estimable:

| Estimable | Not Estimable |
|-----------|---------------|
| Edges (density) | Geodesic distance |
| Degree distribution | Betweenness centrality |
| Triangles (from alter-alter ties) | k-cycles (k > 3) |
| Attribute mixing | Global clustering coefficient |

## References

1. Krivitsky, P.N. & Morris, M. (2017). Inference for social network models from egocentrically sampled data, with application to understanding persistent racial disparities in HIV prevalence in the US. *Annals of Applied Statistics*, 11(1), 427-455.

2. Handcock, M.S. & Gile, K.J. (2010). Modeling social networks from sampled data. *Annals of Applied Statistics*, 4(1), 5-25.

3. Morris, M. (2004). *Network Epidemiology: A Handbook for Survey Design and Data Collection*. Oxford University Press.

4. Krivitsky, P.N., Hunter, D.R., Morris, M., & Klumb, C. (2023). ergm 4: New features for analyzing exponential-family random graph models. *Journal of Statistical Software*, 105(6), 1-44.

5. Smith, J.A. (2012). Macrostructure from microstructure: Generating whole systems from ego networks. *Sociological Methodology*, 42(1), 155-205.
