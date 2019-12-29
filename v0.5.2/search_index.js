var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#AmplNLReader.jl-documentation-1",
    "page": "Home",
    "title": "AmplNLReader.jl documentation",
    "category": "section",
    "text": "This package provides an interface to optimization problems modeled in the AMPL modeling language.Currently, only smooth problems are supported.This package implements the NLPModels.jl API.Interaction with AMPL models is realized by way of the AMPL Solver Library (ASL) as implemented by David Gay.Most, but not all, AMPL models come in the form of a model file (problem.mod) and a data file (problem.dat). For this reason, AmplNLReader assumes that an nl file is available. An nl file is the result of decoding a model and, possibly, a data file to instantiate a concrete problem. One can be generated from the command line usingampl -ogproblem problem.mod problem.datwhere problem should be replaced with your problem\'s name.For an introduction to the AMPL modeling language and background on the ASL, seeR. Fourer, D. M. Gay, and B. W. Kernighan, AMPL: A Mathematical Programming Language, Management Science 36, pp. 519-554, 1990.\nR. Fourer, D. M. Gay, and B. W. Kernighan, AMPL: A Modeling Language for Mathematical Programming, Duxbury Press / Brooks/Cole Publishing Company, 2003\nD. M. Gay, Hooking your Solver to AMPL, Technical Report 97-4-06, Computing Sciences Research Center, Bell Laboratories, Murray Hill, NJ, 1997.\nD. Orban, The Lightning AMPL Tutorial. A Guide for Nonlinear Optimization Users, GERAD Technical Report G-2009-66, 2009."
},

{
    "location": "#Installing-1",
    "page": "Home",
    "title": "Installing",
    "category": "section",
    "text": "The following command should automatically download and install AmplNLReader.jl and its dependencies (Use ] to enter pkg> mode):pkg> add AmplNLReader\npkg> build AmplNLReader\npkg> test AmplNLReader"
},

{
    "location": "#Usage-1",
    "page": "Home",
    "title": "Usage",
    "category": "section",
    "text": "julia> using AmplNLReader\n\njulia> hs33 = AmplModel(\"hs033.nl\")\nMinimization problem hs033.nl\nnvar = 3, ncon = 2 (0 linear)\n\njulia> print(hs33)\nMinimization problem hs033.nl\nnvar = 3, ncon = 2 (0 linear)\nlvar = 1x3 Array{Float64,2}:\n 0.0  0.0  0.0\nuvar = 1x3 Array{Float64,2}:\n Inf  Inf  5.0\nlcon = 1x2 Array{Float64,2}:\n -Inf  4.0\nucon = 1x2 Array{Float64,2}:\n 0.0  Inf\nx0 = 1x3 Array{Float64,2}:\n 0.0  0.0  3.0\ny0 = 1x2 Array{Float64,2}:\n -0.0  -0.0There is support for holding multiple models in memory simultaneously. This should be transparent to the user.Check the NLPModels API for details on the complete API."
},

{
    "location": "#Contents-1",
    "page": "Home",
    "title": "Contents",
    "category": "section",
    "text": ""
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": ""
},

{
    "location": "api/#NLPModels-API-1",
    "page": "API",
    "title": "NLPModels API",
    "category": "section",
    "text": "AmplNLReader implements the NLPModels API."
},

{
    "location": "api/#Internal-1",
    "page": "API",
    "title": "Internal",
    "category": "section",
    "text": "Modules = [AmplNLReader]\nPages   = [\"AmplNLReader.jl\"]\nOrder   = [:function]"
},

{
    "location": "reference/#",
    "page": "Reference",
    "title": "Reference",
    "category": "page",
    "text": ""
},

{
    "location": "reference/#Reference-1",
    "page": "Reference",
    "title": "Reference",
    "category": "section",
    "text": ""
},

]}
