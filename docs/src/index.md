# AmplNLReader.jl documentation

This package provides an interface to optimization problems modeled in the
[AMPL](https://ampl.com) modeling language.

Currently, only smooth problems are supported.

This package implements the
[NLPModels.jl](https://github.com/JuliaSmoothOptimizers/NLPModels.jl/) API.

Interaction with AMPL models is realized by way of the AMPL Solver Library (ASL)
as [implemented](http://netlib.org/ampl/) by David Gay.

Most, but not all, AMPL models come in the form of a model file (`problem.mod`)
and a data file (`problem.dat`).
For this reason, AmplNLReader assumes that an `nl` file is available.
An `nl` file is the result of decoding a model and, possibly, a data file to
instantiate a concrete problem.
One can be generated from the command line using
```
ampl -ogproblem problem.mod problem.dat
```
where `problem` should be replaced with your problem's name.

For an introduction to the AMPL modeling language and background on the ASL, see

* R. Fourer, D. M. Gay, and B. W. Kernighan, [AMPL: A Mathematical Programming Language](https://ampl.com/REFS/amplmod.pdf), Management Science 36, pp. 519-554, 1990.
* R. Fourer, D. M. Gay, and B. W. Kernighan, [AMPL: A Modeling Language for Mathematical Programming](https://ampl.com/resources/the-ampl-book/chapter-downloads/), Duxbury Press / Brooks/Cole Publishing Company, 2003
* D. M. Gay, [Hooking your Solver to AMPL](https://ampl.com/REFS/hooking2.pdf), Technical Report 97-4-06, Computing Sciences Research Center, Bell Laboratories, Murray Hill, NJ, 1997.
* D. Orban, [The Lightning AMPL Tutorial. A Guide for Nonlinear Optimization Users](https://www.gerad.ca/en/papers/G-2009-66), [GERAD](https://www.gerad.ca/en) Technical Report G-2009-66, 2009.

## Installing

The following command should automatically download and install AmplNLReader.jl and its
dependencies (Use `]` to enter `pkg>` mode):
````julia
pkg> add AmplNLReader
pkg> build AmplNLReader
pkg> test AmplNLReader
````

## Usage

````julia
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

Check the [NLPModels
API](https://JuliaSmoothOptimizers.github.io/NLPModels.jl/stable/api/) for details on the complete API.

## Contents

```@contents
```
