export AmplMPECModel

@doc raw"""
    AmplMPECModel

Smooth reformulation for an Ampl model with complementarity constraints.

The nonlinear program
```math
\begin{aligned}
       min_x  \quad &  f(x)\\
\mathrm{s.t.} \quad &  c_L ≤ c(x) ≤ c_U,\\
                    &  g_L ≤ g(x) ⟂ x ≥ x_L,
\end{aligned}
```
is reformulated as
```math
\begin{aligned}
       min_x  \quad &  f(x)\\
\mathrm{s.t.} \quad &  c_L ≤ c(x) ≤ c_U,\\
                    &  g_L ≤ g(x)  \\
                    &  x_L ≤ x \\
                    &  Diag(g(x)) x ≤ 0
\end{aligned}
```

Return the original model if it does not have complementarity constraints.

"""
struct AmplMPECModel <: AbstractNLPModel{Float64, Vector{Float64}}
  mp::AmplModel
  meta::NLPModelMeta
  counters::Counters
  ind_cc_cons::Vector{Int}         # Dimension [m_cc]
  buffer_c1::Vector{Float64}       # Dimension [m]
  buffer_c2::Vector{Float64}       # Dimension [m]
  buffer_vars::Vector{Float64}     # Dimension [n]
  buffer_jac::Vector{Float64}      # Dimension [nnzj]
  # Sparsity pattern of complementarity constraints in CSR format.
  Jccp::Vector{Int}                # Dimension [m_cc]
  Jccj::Vector{Int}                # Dimension [nnzj_cc]
  Jccv::Vector{Int}                # Dimension [nnzj_cc]
end

function AmplMPECModel(mp::AmplModel)
  if mp.meta.n_cc == 0
    @warn("Model does not have any MPEC constraints, returning original model.")
    return mp
  end
  n_cc = mp.meta.n_cc
  cvar = mp.meta.cvar
  ind_cc_cons = findall(cvar .> 0)
  @lencheck mp.meta.n_cc ind_cc_cons

  n = NLPModels.get_nvar(mp)
  m = NLPModels.get_ncon(mp)
  lvar = NLPModels.get_lvar(mp)
  uvar = NLPModels.get_uvar(mp)
  lcon = NLPModels.get_lcon(mp)
  ucon = NLPModels.get_ucon(mp)
  nnzj = NLPModels.get_nnzj(mp)

  # Check MPEC model has been formatted by Ampl.
  for i = 1:n_cc
    ic = ind_cc_cons[i]
    @assert isinf(uvar[cvar[ic]])
    @assert isinf(ucon[ic])
  end

  buffer_c1 = zeros(m)
  buffer_c2 = zeros(m)
  buffer_vars = zeros(n)
  buffer_jac = zeros(nnzj)

  # Reverse mapping
  cnt = 1
  ind_cc = zeros(Int, m)
  for i = 1:m
    if cvar[i] > 0
      ind_cc[i] = cnt
      cnt += 1
    end
  end

  # Get sparsity pattern of original problem.
  rows = zeros(Int, nnzj)
  cols = zeros(Int, nnzj)
  NLPModels.jac_structure!(mp, rows, cols)
  # Analyze sparsity pattern of complementarity constraints
  cnt = 1
  # Count number of nonzeroes in each row
  Jccp = zeros(Int, n_cc + 1)
  # Terms Diag(x) * ∇g(x)
  for (i, j) in zip(rows, cols)
    if ind_cc[i] > 0
      Jccp[ind_cc[i]] += 1
    end
  end
  # Terms Diag(g(x))
  Jccp[1:n_cc] .+= 1
  nnzcc = sum(Jccp)
  # Cumsum nnz per row
  cnt = 1
  for i = 1:n_cc
    tmp = Jccp[i]
    Jccp[i] = cnt
    cnt += tmp
  end
  Jccp[n_cc + 1] = nnzcc + 1
  # Column indexes
  Jccj = zeros(Int, nnzcc)
  # Stores position in original Jacobian (-1 if terms associated to Diag(g(x)))
  Jccv = zeros(Int, nnzcc)
  # Terms Diag(x) * ∇g(x)
  for k = 1:nnzj
    i, j = rows[k], cols[k]
    ic = ind_cc[i]
    if cvar[i] > 0
      dest = Jccp[ic]
      Jccj[dest] = j
      Jccv[dest] = k
      Jccp[ic] += 1
    end
  end
  # Terms Diag(g(x))
  for i = 1:n_cc
    ic = ind_cc_cons[i]
    dest = Jccp[i]
    Jccj[dest] = cvar[ic]
    Jccv[dest] = -1   # store -1 as term is associated to Diag(g(x))
    Jccp[i] += 1
  end
  # Reorder Jccp
  last = 1
  for i = 1:(n_cc + 1)
    tmp = Jccp[i]
    Jccp[i] = last
    last = tmp
  end

  lcon_cc = fill(-Inf, n_cc)
  ucon_cc = fill(0.0, n_cc)
  # TODO: find a better initialization for complementarities' multipliers
  y0_cc = fill(0.0, n_cc)

  meta = NLPModels.NLPModelMeta(
    n;
    lvar = lvar,
    uvar = uvar,
    x0 = NLPModels.get_x0(mp),
    y0 = [NLPModels.get_y0(mp); y0_cc],
    nnzj = NLPModels.get_nnzj(mp) + nnzcc,
    nnzh = NLPModels.get_nnzh(mp) + nnzcc - n_cc,
    ncon = m + n_cc,
    lcon = [lcon; lcon_cc],
    ucon = [ucon; ucon_cc],
    minimize = mp.meta.minimize,
    name = "MPEC-" * mp.meta.name,
  )

  return AmplMPECModel(
    mp,
    meta,
    Counters(),
    ind_cc_cons,
    buffer_c1,
    buffer_c2,
    buffer_vars,
    buffer_jac,
    Jccp,
    Jccj,
    Jccv,
  )
end

function NLPModels.obj(nlp::AmplMPECModel, x::AbstractVector)
  @lencheck nlp.meta.nvar x
  return NLPModels.obj(nlp.mp, x)
end

function NLPModels.grad!(nlp::AmplMPECModel, x::AbstractVector, g::AbstractVector)
  @lencheck nlp.meta.nvar x g
  return NLPModels.grad!(nlp.mp, x, g)
end

function NLPModels.cons!(nlp::AmplMPECModel, x::AbstractVector, c::AbstractVector)
  @lencheck nlp.meta.nvar x
  @lencheck nlp.meta.ncon c

  m = NLPModels.get_ncon(nlp.mp)
  n_cc = nlp.mp.meta.n_cc
  cvar = nlp.mp.meta.cvar
  lcon = NLPModels.get_lcon(nlp.mp)
  lvar = NLPModels.get_lvar(nlp.mp)
  # Evaluate regular constraints
  c_reg = nlp.buffer_c1
  NLPModels.cons!(nlp.mp, x, c_reg)
  c[1:m] .= c_reg
  # Evaluate complementarity constraints
  for i = 1:n_cc
    ic = nlp.ind_cc_cons[i]
    c[m + i] = (c[ic] - lcon[ic]) * (x[cvar[ic]] - lvar[cvar[ic]])
  end
  return c
end

function NLPModels.jac_structure!(
  nlp::AmplMPECModel,
  rows::AbstractVector{<:Integer},
  cols::AbstractVector{<:Integer},
)
  n = NLPModels.get_nvar(nlp.mp)
  m = NLPModels.get_ncon(nlp.mp)
  nnzj = NLPModels.get_nnzj(nlp.mp)
  rows_reg, cols_reg = view(rows, 1:nnzj), view(cols, 1:nnzj)
  NLPModels.jac_structure!(nlp.mp, rows_reg, cols_reg)

  n_cc = nlp.mp.meta.n_cc
  cvar = nlp.mp.meta.cvar

  # Sparsity pattern of complementarity constraints is encoded
  # in CSR format in (Jccp, Jccj)
  cnt = nnzj
  for i = 1:n_cc
    for k = nlp.Jccp[i]:(nlp.Jccp[i + 1] - 1)
      cnt += 1
      rows[cnt] = i + m
      cols[cnt] = nlp.Jccj[k]
    end
  end

  return rows, cols
end

function NLPModels.jac_coord!(nlp::AmplMPECModel, x::AbstractVector, vals::AbstractVector)
  @lencheck nlp.meta.nvar x
  @lencheck nlp.meta.nnzj vals

  n = NLPModels.get_nvar(nlp.mp)
  m = NLPModels.get_ncon(nlp.mp)
  nnzj = NLPModels.get_nnzj(nlp.mp)
  lcon = NLPModels.get_lcon(nlp.mp)
  lvar = NLPModels.get_lvar(nlp.mp)

  vals_reg = view(vals, 1:nnzj)
  NLPModels.jac_coord!(nlp.mp, x, vals_reg)

  cnt = nnzj

  c = nlp.buffer_c2
  NLPModels.cons!(nlp.mp, x, c)

  # Add contribution of complementarity constraints
  n_cc = nlp.mp.meta.n_cc
  cvar = nlp.mp.meta.cvar
  v = nlp.buffer_c1
  Jv = nlp.buffer_vars
  for i = 1:n_cc
    fill!(v, 0.0)
    fill!(Jv, 0.0)
    ic = nlp.ind_cc_cons[i]
    v[ic] = 1.0
    NLPModels.jtprod!(nlp.mp, x, v, Jv)
    # Unpack Jv
    for k = nlp.Jccp[i]:(nlp.Jccp[i + 1] - 1)
      cnt += 1
      j = nlp.Jccj[k]
      if nlp.Jccv[k] >= 1
        # Terms Diag(x) ∇g(x)
        vals[cnt] = (x[cvar[ic]] - lvar[cvar[ic]]) * Jv[j]
      else
        # Terms Diag(g(x))
        vals[cnt] = (c[ic] - lcon[ic])
      end
    end
  end
  return vals
end

function NLPModels.jprod!(
  nlp::AmplMPECModel,
  x::AbstractVector,
  v::AbstractVector,
  Jv::AbstractVector,
)
  @lencheck nlp.meta.nvar x v
  @lencheck nlp.meta.ncon Jv

  n_cc = nlp.mp.meta.n_cc
  cvar = nlp.mp.meta.cvar

  m = NLPModels.get_ncon(nlp.mp)
  lcon = NLPModels.get_lcon(nlp.mp)
  lvar = NLPModels.get_lvar(nlp.mp)

  Jv_buf = nlp.buffer_c1

  NLPModels.jprod!(nlp.mp, x, v, Jv_buf)
  Jv[1:m] .= Jv_buf

  # Add contribution of complementarity constraints
  c = nlp.buffer_c2
  NLPModels.cons!(nlp.mp, x, c)
  for i = 1:n_cc
    ic = nlp.ind_cc_cons[i]
    xj = x[cvar[ic]] - lvar[cvar[ic]]
    Jv[m + i] = Jv[ic] * xj + (c[ic] - lcon[ic]) * v[cvar[ic]]
  end

  return Jv
end

function NLPModels.jtprod!(
  nlp::AmplMPECModel,
  x::AbstractVector,
  v::AbstractVector,
  Jtv::AbstractVector,
)
  @lencheck nlp.meta.nvar x Jtv
  @lencheck nlp.meta.ncon v

  m = NLPModels.get_ncon(nlp.mp)
  lcon = NLPModels.get_lcon(nlp.mp)
  lvar = NLPModels.get_lvar(nlp.mp)

  n_cc = nlp.mp.meta.n_cc
  cvar = nlp.mp.meta.cvar

  v_buf = nlp.buffer_c1
  v_buf .= @view v[1:m]
  for i = 1:n_cc
    ic = nlp.ind_cc_cons[i]
    v_buf[ic] += (x[cvar[ic]] - lvar[cvar[ic]])
  end
  NLPModels.jtprod!(nlp.mp, x, v_buf, Jtv)

  # Add contribution of complementarity constraints
  c = nlp.buffer_c2
  NLPModels.cons!(nlp.mp, x, c)
  for i = 1:n_cc
    ic = nlp.ind_cc_cons[i]
    Jtv[cvar[ic]] += (c[ic] - lcon[ic])
  end
  return Jtv
end

function NLPModels.hess_structure!(nlp::AmplMPECModel, hrows::AbstractVector, hcols::AbstractVector)
  n = NLPModels.get_nvar(nlp.mp)
  n_cc = nlp.mp.meta.n_cc
  cvar = nlp.mp.meta.cvar
  # Evaluate Hessian structure of original model
  nnzh = NLPModels.get_nnzh(nlp.mp)
  rows_reg, cols_reg = view(hrows, 1:nnzh), view(hcols, 1:nnzh)
  NLPModels.hess_structure!(nlp.mp, rows_reg, cols_reg)

  cnt = nnzh
  # Add complementarity contributions ∇g(x)' ∇h(x) + ∇h(x)' ∇g(x)
  # with h(x) := x[cvar[icc]]
  for i = 1:n_cc
    ic = nlp.ind_cc_cons[i]
    for k = nlp.Jccp[i]:(nlp.Jccp[i + 1] - 1)
      j = nlp.Jccj[k]
      # Only store lower-triangular entries
      if nlp.Jccv[k] >= 1
        if cvar[ic] >= j
          cnt += 1
          hrows[cnt] = cvar[ic]
          hcols[cnt] = j
        else
          cnt += 1
          hrows[cnt] = j
          hcols[cnt] = cvar[ic]
        end
      end
    end
  end
  return hrows, hcols
end

function NLPModels.hess_coord!(
  nlp::AmplMPECModel,
  x::AbstractVector,
  vals::AbstractVector;
  obj_weight::Real = one(eltype(x)),
)
  @lencheck nlp.meta.nvar x
  @lencheck nlp.meta.nnzh vals
  NLPModels.hess_coord!(nlp.mp, x, vals; obj_weight = obj_weight)
  return vals
end

function NLPModels.hess_coord!(
  nlp::AmplMPECModel,
  x::AbstractVector,
  y::AbstractVector,
  vals::AbstractVector;
  obj_weight::Real = one(eltype(x)),
)
  @lencheck nlp.meta.nvar x
  @lencheck nlp.meta.ncon y
  @lencheck nlp.meta.nnzh vals

  m = NLPModels.get_ncon(nlp.mp)
  lvar = NLPModels.get_lvar(nlp.mp)
  nnzh = NLPModels.get_nnzh(nlp.mp)
  n_cc = nlp.mp.meta.n_cc
  cvar = nlp.mp.meta.cvar

  y_buf = nlp.buffer_c1
  y_buf .= @view y[1:m]
  for i = 1:n_cc
    ic = nlp.ind_cc_cons[i]
    y_buf[ic] += y[m + i] * (x[cvar[ic]] - lvar[cvar[ic]])
  end
  NLPModels.hess_coord!(nlp.mp, x, y_buf, vals; obj_weight = obj_weight)
  # Add complementarity contributions
  jac_buf = nlp.buffer_jac
  NLPModels.jac_coord!(nlp.mp, x, jac_buf)
  cnt = nnzh
  for i = 1:n_cc
    ic = nlp.ind_cc_cons[i]
    mul_i = y[m + i]
    for k = nlp.Jccp[i]:(nlp.Jccp[i + 1] - 1)
      j = nlp.Jccj[k]
      if nlp.Jccv[k] >= 1
        if cvar[ic] >= j
          cnt += 1
          coef = (cvar[ic] == j) ? 2.0 : 1.0
          vals[cnt] = coef * mul_i * jac_buf[nlp.Jccv[k]]
        else
        end
      end
    end
  end
  return vals
end
