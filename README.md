# AmplNLReader.jl: A [Julia](http://julialang.org) interface to [AMPL](http://www.ampl.com)

## How to Cite

If you use AmplNLReader.jl in your work, please cite using the format given in [CITATION.bib](https://github.com/JuliaSmoothOptimizers/AmplNLReader.jl/blob/main/CITATION.bib).

### Stable release [![Github release](https://img.shields.io/github/release/JuliaSmoothOptimizers/AmplNLReader.jl.svg)](https://github.com/JuliaSmoothOptimizers/AmplNLReader.jl/releases/latest)

- Documentation: [![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSmoothOptimizers.github.io/AmplNLReader.jl/stable)

### Development version

- Documentation: [![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaSmoothOptimizers.github.io/AmplNLReader.jl/latest)
- Tests:
    - ![CI](https://github.com/JuliaSmoothOptimizers/AmplNLReader.jl/workflows/CI/badge.svg?branch=main) (Linux, macOS, Windows)
    - [![Build Status](https://api.cirrus-ci.com/github/JuliaSmoothOptimizers/AmplNLReader.jl.svg)](https://cirrus-ci.com/github/JuliaSmoothOptimizers/AmplNLReader.jl) (FreeBSD)
    - [![codecov](https://codecov.io/gh/JuliaSmoothOptimizers/AmplNLReader.jl/branch/main/graph/badge.svg?token=KwtSr5vCBr)](https://codecov.io/gh/JuliaSmoothOptimizers/AmplNLReader.jl)

## How to Install

At the Julia prompt,

````JULIA
pkg> add AmplNLReader
````

## Testing

````JULIA
pkg> test AmplNLReader
````

## Creating a Model

For an introduction to the AMPL modeling language, see

* R. Fourer, D. M. Gay, and B. W. Kernighan, [AMPL: A Mathematical Programming Language](http://ampl.com/REFS/amplmod.pdf), Management Science 36, pp. 519-554, 1990.
* R. Fourer, D. M. Gay, and B. W. Kernighan, [AMPL: A Modeling Language for Mathematical Programming](http://ampl.com/BOOK/download.html), Duxbury Press / Brooks/Cole Publishing Company, 2003.
* D. Orban, [The Lightning AMPL Tutorial. A Guide for Nonlinear Optimization Users](https://gerad.ca/en/papers/G-2009-66), [GERAD](http://www.gerad.ca) Technical Report G-2009-66, 2009.

Suppose you have an AMPL model represented by the model and data files `mymodel.mod` and `mymodel.dat`. Decode this model as a so-called `nl` file using

    ampl -ogmymodel mymodel.mod mymodel.dat

For example:

````Julia
julia> using AmplNLReader

julia> hs33 = AmplModel("hs033.nl")
Minimization problem hs033.nl
nvar = 3, ncon = 2 (0 linear)

julia> print(hs33)
Minimization problem hs033.nl
nvar = 3, ncon = 2 (0 linear)
lvar = 1x3 Array{Float64,2}:
 0.0  0.0  0.0
uvar = 1x3 Array{Float64,2}:
 Inf  Inf  5.0
lcon = 1x2 Array{Float64,2}:
 -Inf  4.0
ucon = 1x2 Array{Float64,2}:
 0.0  Inf
x0 = 1x3 Array{Float64,2}:
 0.0  0.0  3.0
y0 = 1x2 Array{Float64,2}:
 -0.0  -0.0
````

There is support for holding multiple models in memory simultaneously. This should be transparent to the user.

## Optimization Problems

`AmplNLReader.jl` currently focuses on continuous problems conforming to [`NLPModels.jl`](https://github.com/JuliaSmoothOptimizers/NLPModels.jl).

`AmplModel` objects support all methods associated to `NLPModel` objects.
Please see the [`NLPModels.jl` documentation](https://juliasmoothoptimizers.github.io/NLPModels.jl/latest) for more information.
The following table lists extra methods associated to an `AmplModel`.
See [Hooking your Solver to AMPL](http://ampl.com/REFS/hooking2.pdf) for background.

Method                          | Notes
--------------------------------|--------------------------------
`write_sol(nlp, msg, x, y)`     | Write primal and dual solutions to file

## Missing Methods

* methods for LPs (sparse cost, sparse constraint matrix)
* methods to check optimality conditions.
