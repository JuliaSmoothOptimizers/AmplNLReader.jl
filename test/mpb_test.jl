using Ipopt
using MathProgBase
using NLPModels

# pass an AmplModel to IPOPT
path = dirname(@__FILE__)
nlp = AmplModel(joinpath(path, "hs006.nl"))
model = NLPtoMPB(nlp, IpoptSolver())
@assert isa(model, Ipopt.IpoptMathProgModel)
MathProgBase.optimize!(model)
@assert MathProgBase.getobjval(model) ≈ 0.0
