using Documenter, AmplNLReader, NLPModels

makedocs(
  modules = [AmplNLReader],
  doctest = true,
  linkcheck = false,
  strict = true,
  format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    assets = ["assets/style.css"],
  ),
  sitename = "AmplNLReader.jl",
  pages = Any["Home" => "index.md", "API" => "api.md", "Reference" => "reference.md"],
)

deploydocs(
  repo = "github.com/JuliaSmoothOptimizers/AmplNLReader.jl.git",
  push_preview = true,
  devbranch = "main",
)
