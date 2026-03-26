# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ERGMEgo.jl is a Julia port of the R `ergm.ego` package (StatNet collection) for fitting Exponential Random Graph Models (ERGMs) to egocentrically sampled network data. It enables inference about population-level network properties from ego samples where each sampled "ego" reports their local network (alters and alter-alter ties).

## Development Commands

- **Run tests**: `julia --project -e 'using Pkg; Pkg.test()'`
- **Install dependencies**: `julia --project -e 'using Pkg; Pkg.instantiate()'`
- **Load package locally**: `julia --project -e 'using ERGMEgo'`
- **Build docs**: `julia --project=docs docs/make.jl`

## Architecture

The entire package lives in a single module file: `src/ERGMEgo.jl`. It is organized into these sections (in order):

1. **Data Structures** — `EgoNetwork{T}` (single ego observation with alters, alter ties, attributes, weights) and `EgoData{T}` (collection of ego networks with population size, sampling weights, survey design info). Helper functions: `n_alters`, `alter_degree`, `ego_degree`, `n_alter_ties`, `summary_stats`.
2. **Data Preparation** — `as_egodata` (DataFrame to EgoData conversion) and `ego_design` (attach survey design/weights).
3. **Ego-Specific ERGM Terms** — All subtypes of `AbstractERGMTerm` (from the ERGM package), each with `name()` and `compute()` methods: `EgoEdges`, `EgoNodeMatch`, `EgoDegree`, `EgoGWDegree`, `EgoTriangle`, `EgoMixingMatrix`.
4. **Model and Estimation** — `EgoERGMModel`, `EgoERGMResult`, `ergm_ego` (main entry point, alias `fit_ego_ergm`), and `ego_mple` (pseudo-likelihood estimation via gradient descent).
5. **Population Size Estimation** — `estimate_popsize` with Horvitz-Thompson and capture-recapture methods.
6. **Simulation** — `simulate_ego_sample` generates ego samples from a complete `Network`.
7. **Diagnostics** — `ego_gof` for goodness-of-fit checks.

## Key Dependencies

- **ERGM.jl** — Provides `AbstractERGMTerm` base type that all ego terms extend
- **Network.jl** — Network data structure used in simulation (`simulate_ego_sample`)
- **DataFrames.jl** — Used in `as_egodata` for DataFrame ingestion
- **Graphs.jl** — Graph primitives (`nv`, `neighbors`, `has_edge`) used in simulation
- **Optim.jl**, **Distributions.jl**, **StatsBase.jl** — Statistical computing support

## Conventions

- All code resides in a single file (`src/ERGMEgo.jl`); there are no `include` statements.
- Types are parametric on vertex ID type `T` (e.g., `EgoNetwork{T}`, `EgoData{T}`).
- Each ERGM term is a struct subtyping `AbstractERGMTerm` with two required methods: `name(t)::String` and `compute(t, ed::EgoData)::Float64`.
- Weighted statistics use `ed.sampling_weights` paired with `ed.egos` via `zip`.
- Networks are assumed undirected (alter ties divided by 2 in `n_alter_ties`).
- Julia docstrings are used on all exported functions and types.
- Requires Julia >= 1.9.
