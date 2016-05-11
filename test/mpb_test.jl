using Ipopt
using MathProgBase

# pass an AmplModel to IPOPT
nlp = AmplModel("hs006.nl")
model = NLPtoMPB(nlp, IpoptSolver())
@assert isa(model, Ipopt.IpoptMathProgModel)
MathProgBase.optimize!(model)
@assert MathProgBase.getobjval(model) ≈ 0.0
