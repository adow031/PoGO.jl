using Documenter, PoGO

makedocs(
    sitename = "PoGO.jl: Piecewise approximation for Global Optimization",
    modules = [PoGO],
    clean = true,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["PoGO.jl" => "index.md", "API Reference" => "api.md"],
)

deploydocs(repo = "github.com/adow031/PoGO.jl.git", devurl = "docs")
