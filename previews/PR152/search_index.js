var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#NLPModels-API","page":"API","title":"NLPModels API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"AmplNLReader implements the NLPModels API.","category":"page"},{"location":"api/#Internal","page":"API","title":"Internal","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [AmplNLReader]\nPages   = [\"AmplNLReader.jl\"]\nOrder   = [:function]","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [AmplNLReader]","category":"page"},{"location":"reference/#AmplNLReader.AmplMPECModel","page":"Reference","title":"AmplNLReader.AmplMPECModel","text":"AmplMPECModel\n\nSmooth reformulation for an Ampl model with complementarity constraints.\n\nThe nonlinear program\n\nbeginaligned\n       min_x  quad   f(x)\nmathrmst quad   c_L  c(x)  c_U\n                      g_L  g(x)  x  x_L\nendaligned\n\nis reformulated as\n\nbeginaligned\n       min_x  quad   f(x)\nmathrmst quad   c_L  c(x)  c_U\n                      g_L  g(x)  \n                      x_L  x \n                      Diag(g(x)) x  0\nendaligned\n\nReturn the original model if it does not have complementarity constraints.\n\n\n\n\n\n","category":"type"},{"location":"reference/#AmplNLReader.AmplNLPMeta","page":"Reference","title":"AmplNLReader.AmplNLPMeta","text":"AmplNLPMeta <: AbstractNLPModelMeta\n\nA composite type that represents the main features of the optimization problem\n\noptimize    obj(x)\nsubject to  lvar ≤    x    ≤ uvar\n            lcon ≤ cons(x) ≤ ucon\n\nwhere x        is an nvar-dimensional vector,       obj      is the real-valued objective function,       cons     is the vector-valued constraint function,       optimize is either \"minimize\" or \"maximize\".\n\nHere, lvar, uvar, lcon and ucon are vectors. Some of their components may be infinite to indicate that the corresponding bound or general constraint is not present.\n\n\n\nAmplNLPMeta(nvar; kwargs...)\n\nCreate an AmplNLPMeta with nvar variables. The following keyword arguments are accepted:\n\nx0: initial guess\nlvar: vector of lower bounds\nuvar: vector of upper bounds\nnbv: number of linear binary variables\nniv: number of linear non-binary integer variables\nnlvb: number of nonlinear variables in both objectives and constraints\nnlvo: number of nonlinear variables in objectives (includes nlvb)\nnlvc: number of nonlinear variables in constraints (includes nlvb)\nnlvbi: number of integer nonlinear variables in both objectives and constraints\nnlvci: number of integer nonlinear variables in constraints only\nnlvoi: number of integer nonlinear variables in objectives only\nnwv: number of linear network (arc) variables\nncon: number of general constraints\ny0: initial Lagrange multipliers\nlcon: vector of constraint lower bounds\nucon: vector of constraint upper bounds\nnnzo: number of nonzeros in all objectives gradients\nnnzj: number of elements needed to store the nonzeros in the sparse Jacobian\nlin_nnzj: number of elements needed to store the nonzeros in the sparse Jacobian of linear constraints\nnln_nnzj: number of elements needed to store the nonzeros in the sparse Jacobian of nonlinear constraints\nnnzh: number of elements needed to store the nonzeros in the sparse Hessian\nnlin: number of linear constraints\nnnln: number of nonlinear general constraints\nnnnet: number of nonlinear network constraints\nnlnet: number of linear network constraints\nlin: indices of linear constraints\nnln: indices of nonlinear constraints\nnnet: indices of nonlinear network constraints\nlnet: indices of linear network constraints\nminimize: true if optimize == minimize\nnlo: number of nonlinear objectives\nislp: true if the problem is a linear program\nn_cc: number of complementarity constraints\ncvar: indices of variables appearing in complementarity constraints (0 if constraint is regular)\nname: problem name\n\n\n\n\n\n","category":"type"},{"location":"#AmplNLReader.jl-documentation","page":"Home","title":"AmplNLReader.jl documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides an interface to optimization problems modeled in the AMPL modeling language.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Currently, only smooth problems are supported.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package implements the NLPModels.jl API.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Interaction with AMPL models is realized by way of the AMPL Solver Library (ASL) as implemented by David Gay.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Most, but not all, AMPL models come in the form of a model file (problem.mod) and a data file (problem.dat). For this reason, AmplNLReader assumes that an nl file is available. An nl file is the result of decoding a model and, possibly, a data file to instantiate a concrete problem. One can be generated from the command line using","category":"page"},{"location":"","page":"Home","title":"Home","text":"ampl -ogproblem problem.mod problem.dat","category":"page"},{"location":"","page":"Home","title":"Home","text":"where problem should be replaced with your problem's name.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For an introduction to the AMPL modeling language and background on the ASL, see","category":"page"},{"location":"","page":"Home","title":"Home","text":"R. Fourer, D. M. Gay, and B. W. Kernighan, AMPL: A Mathematical Programming Language, Management Science 36, pp. 519-554, 1990.\nR. Fourer, D. M. Gay, and B. W. Kernighan, AMPL: A Modeling Language for Mathematical Programming, Duxbury Press / Brooks/Cole Publishing Company, 2003\nD. M. Gay, Hooking your Solver to AMPL, Technical Report 97-4-06, Computing Sciences Research Center, Bell Laboratories, Murray Hill, NJ, 1997.\nD. Orban, The Lightning AMPL Tutorial. A Guide for Nonlinear Optimization Users, GERAD Technical Report G-2009-66, 2009.","category":"page"},{"location":"#Installing","page":"Home","title":"Installing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The following command should automatically download and install AmplNLReader.jl and its dependencies (Use ] to enter pkg> mode):","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add AmplNLReader\npkg> build AmplNLReader\npkg> test AmplNLReader","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using AmplNLReader\n\njulia> hs33 = AmplModel(\"hs033.nl\")\nMinimization problem hs033.nl\nnvar = 3, ncon = 2 (0 linear)\n\njulia> print(hs33)\nMinimization problem hs033.nl\nnvar = 3, ncon = 2 (0 linear)\nlvar = 1x3 Array{Float64,2}:\n 0.0  0.0  0.0\nuvar = 1x3 Array{Float64,2}:\n Inf  Inf  5.0\nlcon = 1x2 Array{Float64,2}:\n -Inf  4.0\nucon = 1x2 Array{Float64,2}:\n 0.0  Inf\nx0 = 1x3 Array{Float64,2}:\n 0.0  0.0  3.0\ny0 = 1x2 Array{Float64,2}:\n -0.0  -0.0","category":"page"},{"location":"","page":"Home","title":"Home","text":"There is support for holding multiple models in memory simultaneously. This should be transparent to the user.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Check the NLPModels API for details on the complete API.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
