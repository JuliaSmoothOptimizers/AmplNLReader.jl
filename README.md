# AmplNLReader.jl: A [Julia](http://julialang.org) interface to [AMPL](http://www.ampl.com)

[![Build Status](https://travis-ci.org/dpo/AmplNLReader.jl.svg?branch=master)](https://travis-ci.org/dpo/AmplNLReader.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/tvi3dng5dio174fi?svg=true)](https://ci.appveyor.com/project/dpo/amplnlreader-jl)
[![Coverage Status](https://coveralls.io/repos/dpo/AmplNLReader.jl/badge.svg?branch=master)](https://coveralls.io/r/dpo/AmplNLReader.jl?branch=master)

## How to Install

At the Julia prompt, clone this repository and build:

````JULIA
julia> Pkg.clone("https://github.com/dpo/AmplNLReader.jl.git")
julia> Pkg.build("AmplNLReader")
````

This will automatically use BinDeps to download and install the AMPL interface library.

## Testing

````JULIA
julia> Pkg.test("AmplNLReader")
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

There is preliminary support for holding multiple models in memory simultaneously. This should be transparent to the user.

## Optimization Problems

`AmplNLReader.jl` currently focuses on continuous problems written in the form

    optimize f(x)  subject to l ≤ x ≤ u,  L ≤ c(x) ≤ U,

where `f` is the objective function, `c` is the (vector-valued) constraint function, `l` and `u` are vectors of lower and upper bounds on the variables, and `L` and `U` are vectors of lower and upper bounds on the general constraints.

## Attributes

`AmplModel` objects have the following attributes:

Attribute   | Type               | Notes
------------|--------------------|------------------------------------
`meta`      | `NLPModelMeta`     | See [NLP.jl](https://github.com/optimizers/NLP.jl)
`__asl`     | `Ptr{Void}`        | A private pointer to the internal ASL data structure

## Methods

The following table lists the methods associated to an `AmplModel`. See [Hooking your Solver to AMPL](http://ampl.com/REFS/hooking2.pdf) for background.

Method                          | Notes
--------------------------------|--------------------------------
`varscale(nlp, s)`              | Scale the vector of variables by the vector `s`
`obj(nlp, x)`                   | Evaluate the objective function at `x`
`grad(nlp, x)`                  | Evaluate the objective function gradient at `x`
`lagscale(nlp, s)`              | Set the scaling factor in the Lagrangian
`conscale(nlp, s)`              | Scale the vector of constraints by the vector `s`
`cons(nlp, x)`                  | Evaluate the vector of constraints at `x`
`jth_con(nlp, x, j)`            | Evaluate the `j`-th constraint at `x`
`jth_congrad(nlp, x, j)`        | Evaluate the `j`-th constraint gradient at `x`
`jth_sparse_congrad(nlp, x, j)` | Evaluate the `j`-th constraint sparse gradient at `x`
`jac_coord(nlp, x)`             | Evaluate the sparse Jacobian of the constraints at `x` in coordinate format
`jac(nlp, x)`                   | Evaluate the sparse Jacobian of the constraints at `x`
`hprod(nlp, x, v, y=y0, w=1)`   | Evaluate the product of the Hessian of the Lagrangian at (`x`,`y`) with `v` using the objective weight `w`
`jth_hprod(nlp, x, v, j)`       | Compute the product of the Hessian of the `j`-th constraint at `x` with `v`
`ghjvprod(nlp, x, g, v)`        | Compute the vector of dot products (`g`, `Hj*v`)
`hess_coord(nlp, x, y=y0, w=1.)`| Evaluate the sparse Hessian of the Lagrangian at (`x`,`y`) using the objective weight `w` in coordinate format
`hess(nlp, x, y=y0, w=1.)`      | Evaluate the sparse Hessian of the Lagrangian at (`x`,`y`) using the objective weight `w`
`write_sol(nlp, msg, x, y)`     | Write primal and dual solutions to file

## Missing Methods

* methods for LPs (sparse cost, sparse contraint matrix)
* methods to check optimality conditions.

This content is released under the [MIT](http://opensource.org/licenses/MIT) License.
<a rel="license" href="http://opensource.org/licenses/MIT">
<img alt="MIT license" height="40" src="http://upload.wikimedia.org/wikipedia/commons/c/c3/License_icon-mit.svg" /></a>
