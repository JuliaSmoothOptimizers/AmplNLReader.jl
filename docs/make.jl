using Documenter, AmplNLReader, NLPModels

makedocs(
  modules = [AmplNLReader],
  doctest = true,
  linkcheck = true,
  format = Documenter.HTML(assets = ["assets/style.css"],
                           ansicolor = true,
                           prettyurls = get(ENV, "CI", nothing) == "true"),
  sitename = "AmplNLReader.jl",
  pages = Any["Home" => "index.md",
              "API" => "api.md",
              "Reference" => "reference.md"],
)

deploydocs(
  repo = "github.com/JuliaSmoothOptimizers/AmplNLReader.jl.git",
  push_preview = true,
  devbranch = "main",
)
