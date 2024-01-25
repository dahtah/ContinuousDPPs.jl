# see documentation at https://juliadocs.github.io/Documenter.jl/stable/

using Documenter, ContinuousDPPs

makedocs(
    modules = [ContinuousDPPs],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Simon Barthelmé",
    sitename = "ContinuousDPPs.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "github.com/simon.barthelme/ContinuousDPPs.jl.git",
    push_preview = true
)
