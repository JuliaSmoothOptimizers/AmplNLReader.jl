require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))

import MathProgSolverInterface

type AmplNLPEvaluator <: MathProgSolverInterface.AbstractNLPEvaluator
  nlp::AmplModel
end


# ASL has methods to speed up initialization if certain functionalities aren't
# needed. Currently we always request everything.
MathProgSolverInterface.initialize(::AmplNLPEvaluator, requested_features) = nothing
MathProgSolverInterface.features_available(::AmplNLPEvaluator) = [:Grad, :Jac, :HessVec, :Hess]

MathProgSolverInterface.eval_f(d::AmplNLPEvaluator, x) = obj(d.nlp, x)

MathProgSolverInterface.eval_g(d::AmplNLPEvaluator, g, x) = copy!(g, cons(d.nlp, x))

MathProgSolverInterface.eval_grad_f(d::AmplNLPEvaluator, g, x) = copy!(g, grad(d.nlp, x))

function MathProgSolverInterface.jac_structure(d::AmplNLPEvaluator)
  rows, cols, vals = jac(d.nlp, d.nlp.x0)
  return rows, cols
end

function MathProgSolverInterface.hesslag_structure(d::AmplNLPEvaluator)
  rows, cols, vals = hess(d.nlp, d.nlp.x0, y=ones(d.nlp.ncon))
  return rows, cols
end


function MathProgSolverInterface.eval_jac_g(d::AmplNLPEvaluator, J, x)
  rows, cols, vals = jac(d.nlp, x)
  copy!(J, vals)
end

# Are there specialized methods for Jac-vec products?
# MathProgSolverInterface.eval_jac_prod(d::AmplNLPEvaluator, J, x)
# MathProgSolverInterface.eval_jac_prod_t(d::AmplNLPEvaluator, J, x)


function MathProgSolverInterface.eval_hesslag_prod(d::AmplNLPEvaluator, h, x, v, σ, μ)
  obj(d.nlp, x) # force hessian evaluation at this point
  result = hprod(d.nlp, x, v, y = -μ, obj_weight = σ)
  copy!(h, result)
end

function MathProgSolverInterface.eval_hesslag(d::AmplNLPEvaluator, H, x, σ, μ)
  obj(d.nlp, x) # force hessian evaluation at this point
  rows, cols, vals = hess(d.nlp, x, y = -μ, obj_weight = σ)
  copy!(H, vals)
end

# How do we extract this?
#MathProgSolverInterface.isobjlinear(d::AmplNLPEvaluator)
#MathProgSolverInterface.isobjquadratic(d::AmplNLPEvaluator)

MathProgSolverInterface.isconstrlinear(d::AmplNLPEvaluator,i) = (i in d.nlp.lin)

function loadamplproblem!(m::MathProgSolverInterface.AbstractMathProgModel, nlp::AmplModel)
  sense = nlp.minimize ? :Min : :Max
  MathProgSolverInterface.loadnonlinearproblem!(m, nlp.nvar, nlp.ncon, nlp.lvar, 
    nlp.uvar, nlp.lcon, nlp.ucon, sense, AmplNLPEvaluator(nlp))
end
