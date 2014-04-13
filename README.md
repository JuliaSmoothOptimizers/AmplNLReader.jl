# ampl.jl: A [Julia](http://julialang.org) interface to [AMPL](http://www.ampl.com)

This is a rudimentary Julia interface to the AMPL Solver Library (ASL). Installing on OSX should be easy using [Homebrew](http://brew.sh):

Make sure you have the ASL:

    brew tap homebrew/science
    brew install asl

Clone this repository somewhere, and say `make`. That should do it.

You can test the installation with `julia test.jl`.

For an introduction to the AMPL modeling language, see

* R. Fourer, D. M. Gay, and B. W. Kernighan, [AMPL: A Mathematical Programming Language](http://ampl.com/REFS/amplmod.pdf), Management Science 36, pp. 519-554, 1990.
* R. Fourer, D. M. Gay, and B. W. Kernighan, [AMPL: A Modeling Language for Mathematical Programming](http://ampl.com/BOOK/download.html), Duxbury Press / Brooks/Cole Publishing Company, 2003.
* D. Orban, [The Lightning AMPL Tutorial. A Guide for Nonlinear Optimization Users](http://www.gerad.ca/fichiers/cahiers/G-2009-66.pdf), [GERAD](http://www.gerad.ca) Technical Report G-2009-66, 2009.

## Creating a Model

Suppose you have an AMPL model represented by the model and data files `mymodel.mod` and `mymodel.dat`. Decode this model as a so-called `nl` file using

    ampl -ogmymodel mymodel.mod mymodel.dat

In Julia, create an instance of `AmplModel` representing your model

    julia> include("ampl.jl")
    julia> mymodel = AmplModel("mymodel.nl")

At this point, you can examine the problem dimensions using `mymodel.nvar`, `mymodel.ncon`, etc.

## Optimization Problems

`ampl.jl` currently focuses on continuous problems written in the form

    minimize f(x)  subject to l ≤ x ≤ u,  L ≤ c(x) ≤ U,

where `f` is the objective function, `l` and `u` are vectors of lower and upper bounds on the variables, and `L` and `U` are vectors of lower and upper bounds on the general constraints.

## Attributes

`AmplModel` objects have the following attributes:

Attribute   | Type               | Notes
------------|--------------------|------------------------------------
`nvar`      | `Int             ` | number of variables
`x0  `      | `Array{Float64,1}` | initial guess
`lvar`      | `Array{Float64,1}` | vector of lower bounds
`uvar`      | `Array{Float64,1}` | vector of upper bounds
`ncon`      | `Int             ` | total number of general constraints
`nlc `      | `Int             ` | number of linear constraints
`nlnc`      | `Int             ` | number of nonlinear network constraints
`y0  `      | `Array{Float64,1}` | initial Lagrange multipliers
`lcon`      | `Array{Float64,1}` | vector of constraint lower bounds
`ucon`      | `Array{Float64,1}` | vector of constraint upper bounds
`lin `      | `Range1{Int64}   ` | indices of linear constraints
`nlin`      | `Range1{Int64}   ` | indices of nonlinear constraints (not network)
`nnzj`      | `Int             ` | number of nonzeros in the sparse Jacobian
`nnzh`      | `Int             ` | number of nonzeros in the sparse Hessian
`minimize`  | `Bool            ` | true if optimize == minimize
`islp`      | `Bool            ` | true if the problem is a linear program
`name`      | `ASCIIString     ` | problem name

## Methods

The following table lists the methods associated to an `AmplModel`. See [Hooking your Solver to AMPL](http://ampl.com/REFS/hooking2.pdf) for background.

Method                          | Notes
--------------------------------|--------------------------------
`varscale(nlp, s)`                | Scale the vector of variables by the vector `s`
`obj(nlp, x)`                     | Evaluate the objective function at `x`
`grad(nlp, x)`                    | Evaluate the objective function gradient at `x`
`lagscale(nlp, s)`                | Set the scaling factor in the Lagrangian
`conscale(nlp, s)`                | Scale the vector of constraints by the vector `s`
`cons(nlp, x)`                    | Evaluate the vector of constraints at `x`
`jth_con(nlp, x, j)`              | Evaluate the `j`-th constraint at `x`
`jth_congrad(nlp, x, j)`          | Evaluate the `j`-th constraint gradient at `x`
`jth_sparse_congrad(nlp, x, j)`   | Evaluate the `j`-th constraint sparse gradient at `x`
`jac(nlp, x)`                     | Evaluate the sparse Jacobian of the constraints at `x`
`hprod(nlp, x, v, y=y0, w=1)` | Evaluate the product of the Hessian of the Lagrangian at (`x`,`y`) with `v` using the objective weight `w`
`jth_hprod(nlp, x, v, j)`         | Compute the product of the Hessian of the `j`-th constraint at `x` with `v`
`ghjvprod(nlp, x, g, v)`          | Compute the vector of dot products (`g`, `Hj*v`)
`hess(nlp, x, y=y0, w=1.)`     | Evaluate the sparse Hessian of the Lagrangian at (`x`,`y`) using the objective weight `w`

## Missing Methods

* methods for LPs (sparse cost, sparse contraint matrix)
* methods to check optimality conditions.

## Todo

* Support for [multiple problems at once](http://ampl.com/REFS/HOOKING/index.html#Multipleproblemsandmultiplethreads)

[![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)](http://www.gnu.org/licenses/gpl.html "GPLv3")
