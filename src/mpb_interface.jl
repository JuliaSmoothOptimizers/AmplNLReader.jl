# An interface to pass a AmplModel to a MathProgBase solver

using MathProgBase

MathProgBase.features_available(::NLPModelEvaluator{AmplModel}) =
  [:Grad, :Jac, :HessVec, :Hess]
