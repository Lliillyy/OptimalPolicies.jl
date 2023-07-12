import Pkg
Pkg.add("Documenter")

using Documenter
using OptimalPolicies

makedocs(
    sitename = "OptimalPolicies.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = ["Introduction" => "index.md", "API" => "api.md"],
)

deploydocs(repo = "github.com/Lliillyy/OptimalPolicies.jl.git", devbranch = "main")
