using Documenter
using MATLABDiffEq

DocMeta.setdocmeta!(MATLABDiffEq, :DocTestSetup, :(using MATLABDiffEq); recursive = true)

makedocs(;
    modules = [MATLABDiffEq],
    sitename = "MATLABDiffEq.jl",
    doctest = true,
    checkdocs = :exports,
    format = Documenter.HTML(;
        canonical = "https://docs.sciml.ai/MATLABDiffEq/stable/",
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = ["Home" => "index.md", "API" => "api.md"],
)

deploydocs(; repo = "github.com/SciML/MATLABDiffEq.jl.git", push_preview = true)
