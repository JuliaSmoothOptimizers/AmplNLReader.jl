# An interface to pass a AmplModel to a MathProgBase solver

using MathProgBase
using AMPLMathProgInterface

export NLPtoMPB


"""Return a `MathProgBase` model corresponding to an `AmplModel`.

The second argument should be a solver instance, e.g., `IpoptSolver()`.
Currently, all models are treated as nonlinear models.
"""
function NLPtoMPB(nlp :: AmplModel, solver :: MathProgBase.AbstractMathProgSolver)
  model = MathProgBase.NonlinearModel(solver)
  AMPLMathProgInterface.loadamplproblem!(model, nlp)
  return model
end
