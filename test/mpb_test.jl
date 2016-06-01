using Ipopt
using MathProgBase
# Explicitly use AMPLMathProgInterface to err if it's not installed
using AMPLMathProgInterface

# pass an AmplModel to IPOPT
path = dirname(@__FILE__)
nlp = AmplModel(joinpath(path, "hs006.nl"))
model = NLPtoMPB(nlp, IpoptSolver())
@assert isa(model, Ipopt.IpoptMathProgModel)
MathProgBase.optimize!(model)
@assert MathProgBase.getobjval(model) ≈ 0.0
