using Documenter
using ERGMEgo

DocMeta.setdocmeta!(ERGMEgo, :DocTestSetup, :(using ERGMEgo); recursive=true)

makedocs(
    sitename = "ERGMEgo.jl",
    modules = [ERGMEgo],
    authors = "Statistical Network Analysis with Julia",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://Statistical-network-analysis-with-Julia.github.io/ERGMEgo.jl",
        edit_link = "main",
    ),
    repo = "https://github.com/Statistical-network-analysis-with-Julia/ERGMEgo.jl/blob/{commit}{path}#{line}",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => [
            "Ego Networks" => "guide/ego_networks.md",
            "Ego Terms" => "guide/terms.md",
            "Population Inference" => "guide/inference.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Terms" => "api/terms.md",
            "Estimation" => "api/estimation.md",
        ],
    ],
    warnonly = [:missing_docs, :docs_block],
)

deploydocs(
    repo = "github.com/Statistical-network-analysis-with-Julia/ERGMEgo.jl.git",
    devbranch = "main",
    versions = [
        "stable" => "dev",
        "dev" => "dev",
    ],
    push_preview = true,
)
