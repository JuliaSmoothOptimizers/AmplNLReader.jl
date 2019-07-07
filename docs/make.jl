using Documenter, AmplNLReader, NLPModels

makedocs(
  modules = [AmplNLReader],
  doctest = true,
  linkcheck = true,
  strict = true,
  assets = ["assets/style.css"],
  format = Documenter.HTML(
              prettyurls = get(ENV, "CI", nothing) == "true"
            ),
  sitename = "AmplNLReader.jl",
  pages = Any["Home" => "index.md",
              "API" => "api.md",
              "Reference" => "reference.md"]
)

deploydocs(deps = nothing, make = nothing,
  repo = "github.com/JuliaSmoothOptimizers/AmplNLReader.jl.git",
  target = "build",
  devbranch = "master"
)
