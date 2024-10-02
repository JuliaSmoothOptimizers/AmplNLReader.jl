
using Test
using NLPModels
using AmplNLReader
using LinearAlgebra
using SparseArrays

function test_mpec_wrapper(nlp)
  mp = AmplNLReader.AmplMPECModel(nlp)

  n = NLPModels.get_nvar(nlp)
  m = NLPModels.get_ncon(nlp)
  nnzj = NLPModels.get_nnzj(nlp)
  lvar = NLPModels.get_lvar(nlp)
  lcon = NLPModels.get_lcon(nlp)
  # Complementarity structure
  n_cc = nlp.meta.n_cc
  cvar = nlp.meta.cvar
  icc = mp.ind_cc_cons

  x0 = ones(n)

  # Test problems dimensions
  @test n == NLPModels.get_nvar(mp)
  @test NLPModels.get_ncon(nlp) + n_cc == NLPModels.get_ncon(mp)

  # Test objective and gradient match
  @test NLPModels.obj(nlp, x0) == NLPModels.obj(mp, x0)
  @test NLPModels.grad(nlp, x0) == NLPModels.grad(mp, x0)

  # Test evaluation of constraints
  ixc = cvar[icc]
  c = NLPModels.cons(nlp, x0)
  c_cc = NLPModels.cons(mp, x0)
  @test c_cc == [c; (x0[ixc] .- lvar[ixc]) .* (c[icc] .- lcon[icc])]

  # Test evaluation of Jacobian
  (I, J) = NLPModels.jac_structure(nlp)
  V = zeros(nnzj)
  # Evaluate Jacobian without complementarities
  NLPModels.jac_coord!(nlp, x0, V)
  jac_nlp = sparse(I, J, V)  # original Jacobian
  # Evaluate Jacobian with complementarities
  nnzj_cc = NLPModels.get_nnzj(mp)
  (I_cc, J_cc) = NLPModels.jac_structure(mp)
  V_cc = zeros(nnzj_cc)
  NLPModels.jac_coord!(mp, x0, V_cc)
  jac_cc = sparse(I_cc, J_cc, V_cc)  # Jacobian with complementarities
  @test size(jac_cc) == (m + n_cc, n)
  # Test non-complementarity part matches original Jacobian
  @test jac_cc[1:m, :] == jac_nlp
  # Test Jacobian of complementarity part
  # ∇(Diag(g(x)) x) = Diag(g(x)) + Diag(x) * ∇g(x)
  G = Diagonal(c[icc] .- lcon[icc])
  H = Diagonal(x0[ixc] .- lvar[ixc])
  ∇G = jac_nlp[icc, :]
  ∇H = sparse(1:n_cc, ixc, ones(n_cc), n_cc, n)
  @test jac_cc[(m + 1):(m + n_cc), :] == G * ∇H + H * ∇G

  # jtprod
  v = ones(m + n_cc)
  Jtv = zeros(n)
  NLPModels.jtprod!(mp, x0, v, Jtv)
  @test Jtv ≈ jac_cc' * v

  # jprod
  v = ones(n)
  Jv = zeros(m + n_cc)
  NLPModels.jprod!(mp, x0, v, Jv)
  @test Jv ≈ jac_cc * v

  # Hessian
  y = ones(m)
  (I, J) = NLPModels.hess_structure(nlp)
  # Evaluate original Hessian
  nnzh = NLPModels.get_nnzh(nlp)
  hvals = zeros(nnzh)
  NLPModels.hess_coord!(nlp, x0, y, hvals)
  hess_nlp = sparse(I, J, hvals, n, n)
  # Evaluate Hessian with complementarity constraints
  y_cc = ones(m + n_cc)
  (I_cc, J_cc) = NLPModels.hess_structure(mp)
  nnzh_cc = NLPModels.get_nnzh(mp)
  hvals_cc = zeros(nnzh_cc)
  NLPModels.hess_coord!(mp, x0, y_cc, hvals_cc)
  hess_cc = sparse(I_cc, J_cc, hvals_cc, n, n)
  # Test expression is correct
  y .= 0.0
  y[icc] .= (x0[ixc] .- lvar[ixc])
  hvals_tmp = zeros(nnzh)
  NLPModels.hess_coord!(nlp, x0, y, hvals_tmp; obj_weight = 0.0)
  ∇2G = sparse(I, J, hvals_tmp, n, n)
  D = Diagonal(y_cc[(m + 1):(m + n_cc)])
  JtJ = LowerTriangular(∇G' * D * ∇H + ∇H' * D * ∇G)
  @test hess_cc == hess_nlp + ∇2G + JtJ
  return
end

path = joinpath(dirname(@__FILE__), "problems")
instance = "bard1.nl"
nlp = AmplModel(joinpath(path, instance))
test_mpec_wrapper(nlp)
