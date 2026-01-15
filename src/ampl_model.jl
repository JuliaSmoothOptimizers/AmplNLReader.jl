export AmplModel, AmplException, write_sol, amplmodel_finalize

# Import methods we override.
import Base.show, Base.print

struct AmplException
  msg::String
end

mutable struct AmplModel <: AbstractNLPModel{Float64, Vector{Float64}}
  meta::AmplNLPMeta   # Problem metadata.
  __asl::Ptr{Cvoid}   # Pointer to internal ASL structure. Do not touch.
  counters::Counters  # Evaluation counters
  safe::Bool          # Always evaluate the objective before the Hessian.

  function AmplModel(stub::AbstractString; safe::Bool = false)

    # check that stub or stub.nl exists
    fname = basename(stub)
    ext = occursin(".", fname) ? split(fname, '.')[end] : ""
    if ext == "nl"
      isfile(stub) || throw(AmplException("cannot find $(stub)"))
    else
      isfile("$(stub).nl") || throw(AmplException("cannot find $(stub).nl"))
    end

    asl = asl_init(stub)
    (asl == C_NULL) && error("Error allocating ASL structure")

    minimize = asl_objtype(asl) == 0
    islp = asl_islp(asl) != 0

    nlo = asl_nlo(asl) |> Int

    nvar = asl_nvar(asl) |> Int
    ncon = asl_ncon(asl) |> Int

    x0 = unsafe_wrap(Array, asl_x0(asl), (nvar,), own = false)
    y0 = unsafe_wrap(Array, asl_y0(asl), (ncon,), own = false)

    lvar = unsafe_wrap(Array, asl_lvar(asl), (nvar,), own = false)
    uvar = unsafe_wrap(Array, asl_uvar(asl), (nvar,), own = false)

    nzo = asl_nzo(asl) |> Int
    nbv = asl_nbv(asl) |> Int
    niv = asl_niv(asl) |> Int
    nlvb = asl_nlvb(asl) |> Int
    nlvo = asl_nlvo(asl) |> Int
    nlvc = asl_nlvc(asl) |> Int
    nlvbi = asl_nlvbi(asl) |> Int
    nlvci = asl_nlvci(asl) |> Int
    nlvoi = asl_nlvoi(asl) |> Int
    nwv = asl_nwv(asl) |> Int
    n_cc = asl_n_cc(asl) |> Int

    lcon = unsafe_wrap(Array, asl_lcon(asl), (ncon,), own = false)
    ucon = unsafe_wrap(Array, asl_ucon(asl), (ncon,), own = false)

    if n_cc > 0
      cvar = unsafe_wrap(Array, asl_cvar(asl), (ncon,), own = false)

      # Check complementarity constraints are well specified:
      cc_cons = cvar .> 0
      @assert all(isinf, ucon[cc_cons])
      @assert all(isfinite, lcon[cc_cons])
    else
      cvar = Int[]
    end

    nlnet = asl_lnc(asl) |> Int
    nnnet = asl_nlnc(asl) |> Int
    nnln = (asl_nlc(asl) |> Int) - nnnet

    nln = 1:nnln
    nnet = (nnln + 1):(nnln + nnnet)
    lnet = (nnln + nnnet + 1):(nnln + nnnet + nlnet)
    lin = (nnln + nnnet + nlnet + 1):ncon

    nnzj = asl_nnzj(asl) |> Int
    lin_nnzj = 0
    for j in lin
      lin_nnzj += asl_sparse_congrad_nnz(asl, j-1) |> Cint
    end
    nln_nnzj = 0
    for j in nln
      nln_nnzj += asl_sparse_congrad_nnz(asl, j-1) |> Cint
    end
    nnzh = (asl_nnzh(asl)) |> Int

    meta = AmplNLPMeta(
      nvar,
      x0 = x0,
      lvar = lvar,
      uvar = uvar,
      nlo = nlo,
      nnzo = nzo,
      ncon = ncon,
      y0 = y0,
      lcon = lcon,
      ucon = ucon,
      nnzj = nnzj,
      lin_nnzj = lin_nnzj,
      nln_nnzj = nln_nnzj,
      nnzh = nnzh,
      nbv = nbv,
      niv = niv,
      nlvb = nlvb,
      nlvo = nlvo,
      nlvc = nlvc,
      nlvbi = nlvbi,
      nlvci = nlvci,
      nlvoi = nlvoi,
      nwv = nwv,
      lin = lin,
      nln = nln,
      nnet = nnet,
      lnet = lnet,
      nlnet = nlnet,
      minimize = minimize,
      islp = islp,
      n_cc = n_cc,
      cvar = cvar,
      name = split(fname, ".")[1],  # do not include path or extension in model name
    )

    nlp = new(meta, asl, Counters(), safe)

    finalizer(amplmodel_finalize, nlp)
    return nlp
  end
end

function check_ampl_model(nlp::AmplModel)
  (nlp.__asl == C_NULL) && throw(AmplException("Uninitialized AMPL model"))
end

# Methods associated to AmplModel instances.

function NLPModels.reset!(nlp::AmplModel)
  reset!(nlp.counters)
  return nlp
end

function write_sol(nlp::AmplModel, msg::String, x::AbstractVector, y::AbstractVector)
  check_ampl_model(nlp)
  length(x) == nlp.meta.nvar || error("x must have length $(nlp.meta.nvar)")
  length(y) == nlp.meta.ncon || error("y must have length $(nlp.meta.ncon)")

  asl_write_sol(nlp.__asl, msg, x, y)
end

function amplmodel_finalize(nlp::AmplModel)
  if nlp.__asl != C_NULL
    asl_finalize(nlp.__asl)::Cvoid
    nlp.__asl = C_NULL
  end
  return
end

# Displaying AmplModel instances.

function show(io::IO, nlp::AmplModel)
  check_ampl_model(nlp)
  show(io, nlp.meta)
end

function print(io::IO, nlp::AmplModel)
  check_ampl_model(nlp)
  print(io, nlp.meta)
end

# Scaling AmplModel instances.

function NLPModels.varscale(nlp::AmplModel, s::Vector{Cdouble})
  check_ampl_model(nlp)
  length(s) >= nlp.meta.nvar || error("s must have length at least $(nlp.meta.nvar)")

  err = Ref{Cint}(0)
  asl_varscale(nlp.__asl, s, err)
  err[] == 0 || throw(AmplException("Error while scaling variables"))
end

NLPModels.varscale(nlp::AmplModel, s::AbstractVector) = varscale(nlp, Vector{Cdouble}(s))

function NLPModels.lagscale(nlp::AmplModel, σ::Float64)
  check_ampl_model(nlp)
  err = Ref{Cint}(0)
  asl_lagscale(nlp.__asl, σ, err)
  err[] == 0 || throw(AmplException("Error while scaling Lagrangian"))
end

function NLPModels.conscale(nlp::AmplModel, s::Vector{Cdouble})
  check_ampl_model(nlp)
  length(s) >= nlp.meta.ncon || error("s must have length at least $(nlp.meta.ncon)")

  err = Ref{Cint}(0)
  asl_conscale(nlp.__asl, s, err)
  err[] == 0 || throw(AmplException("Error while scaling constraints"))
end

NLPModels.conscale(nlp::AmplModel, s::AbstractVector) = conscale(nlp, Vector{Cdouble}(s))

# Evaluating objective, constraints and derivatives.

function NLPModels.obj(nlp::AmplModel, x::Vector{Cdouble})
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Ref{Cint}(0)
  f = asl_obj(nlp.__asl, x, err)
  nlp.counters.neval_obj += 1
  err[] == 0 || throw(AmplException("Error while evaluating objective"))
  return f
end

NLPModels.obj(nlp::AmplModel, x::AbstractVector) = obj(nlp, Vector{Cdouble}(x))

function NLPModels.grad!(nlp::AmplModel, x::Vector{Cdouble}, g::Vector{Cdouble})
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Ref{Cint}(0)
  asl_grad(nlp.__asl, x, g, err)
  nlp.counters.neval_grad += 1
  err[] == 0 || throw(AmplException("Error while evaluating objective gradient"))
  return g
end

function NLPModels.grad!(nlp::AmplModel, x::AbstractVector, g::AbstractVector)
  g_ = Vector{Cdouble}(undef, nlp.meta.nvar)
  grad!(nlp, Vector{Cdouble}(x), g_)
  g[1:(nlp.meta.nvar)] .= g_
  return g
end

function NLPModels.cons!(nlp::AmplModel, x::Vector{Cdouble}, c::Vector{Cdouble})
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Ref{Cint}(0)
  asl_cons(nlp.__asl, x, c, err)
  nlp.counters.neval_cons += 1
  err[] == 0 || throw(AmplException("Error while evaluating constraints"))
  return c
end

function NLPModels.cons!(nlp::AmplModel, x::AbstractVector, c::AbstractVector)
  c_ = Vector{Cdouble}(undef, nlp.meta.ncon)
  cons!(nlp, Vector{Cdouble}(x), c_)
  c[1:(nlp.meta.ncon)] .= c_
  return c
end

function NLPModels.jth_con(nlp::AmplModel, x::Vector{Cdouble}, j::Int)
  check_ampl_model(nlp)
  (1 <= j <= nlp.meta.ncon) || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Ref{Cint}(0)
  cj = asl_jcon(nlp.__asl, x, j-1, err)
  nlp.counters.neval_jcon += 1
  err[] == 0 || throw(AmplException("Error while evaluating $j-th constraint"))
  return cj
end

NLPModels.jth_con(nlp::AmplModel, x::AbstractVector, j::Int) = jth_con(nlp, Vector{Cdouble}(x), j)

function NLPModels.jth_congrad!(nlp::AmplModel, x::Vector{Cdouble}, j::Int, g::Vector{Cdouble})
  check_ampl_model(nlp)
  (1 <= j <= nlp.meta.ncon) || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Ref{Cint}(0)
  asl_jcongrad(nlp.__asl, x, g, j-1, err)
  nlp.counters.neval_jgrad += 1
  err[] == 0 || throw(AmplException("Error while evaluating $j-th constraint gradient"))
  return g
end

function NLPModels.jth_congrad!(nlp::AmplModel, x::AbstractVector, j::Int, g::AbstractVector)
  g_ = Vector{Cdouble}(undef, nlp.meta.nvar)
  jth_congrad!(nlp, Vector{Cdouble}(x), j, g_)
  g[1:(nlp.meta.nvar)] .= g_
  return g
end

function NLPModels.jth_sparse_congrad(nlp::AmplModel, x::Vector{Cdouble}, j::Int)
  check_ampl_model(nlp)
  (1 <= j <= nlp.meta.ncon) || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  nnz = asl_sparse_congrad_nnz(nlp.__asl, j-1)

  err = Ref{Cint}(0)
  inds = Vector{Cint}(undef, nnz)
  vals = Vector{Cdouble}(undef, nnz)
  asl_sparse_congrad(nlp.__asl, x, j-1, inds, vals, err)
  nlp.counters.neval_jgrad += 1
  err[] == 0 || throw(AmplException("Error while evaluating $j-th sparse constraint gradient"))

  # Use 1-based indexing.
  inds .+= Cint(1)
  return sparsevec(inds, vals, nlp.meta.nvar)
end

NLPModels.jth_sparse_congrad(nlp::AmplModel, x::AbstractVector, j::Int) =
  jth_sparse_congrad(nlp, Vector{Cdouble}(x), j)

function NLPModels.jac_structure!(nlp::AmplModel, rows::Vector{Cint}, cols::Vector{Cint})
  asl_jac_structure(nlp.__asl, rows, cols)

  # Use 1-based indexing.
  rows[1:(nlp.meta.nnzj)] .+= Cint(1)
  cols[1:(nlp.meta.nnzj)] .+= Cint(1)
  return rows, cols
end

function NLPModels.jac_structure!(
  nlp::AmplModel,
  rows::AbstractVector{<:Integer},
  cols::AbstractVector{<:Integer},
)
  rows_ = Vector{Cint}(undef, nlp.meta.nnzj)
  cols_ = Vector{Cint}(undef, nlp.meta.nnzj)
  jac_structure!(nlp, rows_, cols_)
  rows[1:(nlp.meta.nnzj)] .= rows_
  cols[1:(nlp.meta.nnzj)] .= cols_
  return rows, cols
end

function NLPModels.jac_coord!(nlp::AmplModel, x::Vector{Cdouble}, vals::Vector{Cdouble})
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  _ = cons(nlp, x)
  nlp.counters.neval_cons -= 1

  err = Ref{Cint}(0)
  asl_jacval(nlp.__asl, x, vals, err)
  nlp.counters.neval_jac += 1
  err[] == 0 || throw(AmplException("Error while evaluating constraints Jacobian"))
  return vals
end

function NLPModels.jac_coord!(
  nlp::AmplModel,
  x::AbstractVector,
  vals::AbstractVector{<:AbstractFloat},
)
  vals_ = Vector{Cdouble}(undef, nlp.meta.nnzj)
  jac_coord!(nlp, Vector{Cdouble}(x), vals_)
  vals[1:(nlp.meta.nnzj)] .= vals_
  return vals
end

function NLPModels.jprod!(nlp::AmplModel, x::AbstractVector, v::AbstractVector, Jv::AbstractVector)
  nlp.counters.neval_jac -= 1
  nlp.counters.neval_jprod += 1
  Jv[1:(nlp.meta.ncon)] = jac(nlp, Vector{Cdouble}(x)) * v
  return Jv
end

function NLPModels.jtprod!(
  nlp::AmplModel,
  x::AbstractVector,
  v::AbstractVector,
  Jtv::AbstractVector,
)
  nlp.counters.neval_jac -= 1
  nlp.counters.neval_jtprod += 1
  Jtv[1:(nlp.meta.nvar)] = jac(nlp, Vector{Cdouble}(x))' * v
  return Jtv
end

function NLPModels.hprod!(
  nlp::AmplModel,
  x::AbstractVector,
  y::Vector{Cdouble},
  v::Vector{Cdouble},
  hv::Vector{Cdouble};
  obj_weight::Float64 = 1.0,
)
  # Note: x is in fact not used in hprod.
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(y) >= nlp.meta.ncon || error("y must have length at least $(nlp.meta.ncon)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")

  if nlp.safe
    _ = obj(nlp, x)
    nlp.counters.neval_obj -= 1
    _ = cons(nlp, x)
    nlp.counters.neval_cons -= 1
  end
  asl_hprod(nlp.__asl, y, v, hv, obj_weight)
  nlp.counters.neval_hprod += 1
  return hv
end

function NLPModels.hprod!(
  nlp::AmplModel,
  x::AbstractVector,
  y::AbstractVector,
  v::AbstractVector,
  hv::AbstractVector;
  obj_weight::Float64 = 1.0,
)
  hv_ = Vector{Cdouble}(undef, nlp.meta.nvar)
  hprod!(nlp, x, Vector{Cdouble}(y), Vector{Cdouble}(v), hv_; obj_weight = obj_weight)
  hv[1:(nlp.meta.nvar)] .= hv_
  return hv
end

function NLPModels.jth_hprod!(
  nlp::AmplModel,
  x::AbstractVector,
  v::Vector{Cdouble},
  j::Int,
  hv::Vector{Cdouble},
)
  # Note: x is in fact not used in hprod.
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")
  (0 <= j <= nlp.meta.ncon) || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")

  # if nlp.safe
  if j == 0
    _ = obj(nlp, Vector{Cdouble}(x))
    nlp.counters.neval_obj -= 1
  else
    _ = cons(nlp, Vector{Cdouble}(x))
    nlp.counters.neval_cons -= 1
  end
  # end
  asl_hvcompd(nlp.__asl, v, hv, j-1)
  nlp.counters.neval_jhprod += 1
  return hv
end

function NLPModels.jth_hprod!(
  nlp::AmplModel,
  x::AbstractVector,
  v::AbstractVector,
  j::Int,
  hv::AbstractVector,
)
  hv_ = Vector{Cdouble}(undef, nlp.meta.nvar)
  jth_hprod!(nlp, x, Vector{Cdouble}(v), j, hv_)
  hv[1:(nlp.meta.nvar)] .= hv_
  return hv
end

function NLPModels.ghjvprod!(
  nlp::AmplModel,
  x::AbstractVector,
  g::Vector{Cdouble},
  v::Vector{Cdouble},
  gHv::Vector{Cdouble},
)
  # Note: x is in fact not used.
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(g) >= nlp.meta.nvar || error("g must have length at least $(nlp.meta.nvar)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")

  if nlp.safe
    _ = cons(nlp, x)
    nlp.counters.neval_cons -= 1
  end
  asl_ghjvprod(nlp.__asl, g, v, gHv)
  nlp.counters.neval_hprod += nlp.meta.ncon
  return gHv
end

function NLPModels.ghjvprod!(
  nlp::AmplModel,
  x::AbstractVector,
  g::AbstractVector,
  v::AbstractVector,
  gHv::AbstractVector,
)
  gHv_ = Vector{Cdouble}(undef, nlp.meta.ncon)
  ghjvprod!(nlp, x, Vector{Cdouble}(g), Vector{Cdouble}(v), gHv_)
  gHv[1:(nlp.meta.ncon)] .= gHv_
  return gHv
end

function NLPModels.hess_structure!(nlp::AmplModel, rows::Vector{Cint}, cols::Vector{Cint})
  # Swap rows and cols to obtain the lower triangle.
  asl_hess_structure(nlp.__asl, cols, rows)

  # Use 1-based indexing.
  cols .+= Cint(1)
  rows .+= Cint(1)
  return (rows, cols)
end

function NLPModels.hess_structure!(
  nlp::AmplModel,
  rows::AbstractVector{<:Integer},
  cols::AbstractVector{<:Integer},
)
  rows_ = Vector{Cint}(undef, nlp.meta.nnzh)
  cols_ = Vector{Cint}(undef, nlp.meta.nnzh)
  hess_structure!(nlp, rows_, cols_)
  rows[1:(nlp.meta.nnzh)] .= rows_
  cols[1:(nlp.meta.nnzh)] .= cols_
  return (rows, cols)
end

function NLPModels.hess_coord!(
  nlp::AmplModel,
  x::AbstractVector,
  y::Vector{Cdouble},
  vals::Vector{Cdouble};
  obj_weight::Float64 = 1.0,
)
  # Note: x is in fact not used.
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(y) >= nlp.meta.ncon || error("y must have length at least $(nlp.meta.ncon)")

  # if nlp.safe
  _ = obj(nlp, x)
  nlp.counters.neval_obj -= 1
  _ = cons(nlp, x)
  nlp.counters.neval_cons -= 1
  # end

  asl_hessval(nlp.__asl, y, obj_weight, vals)
  nlp.counters.neval_hess += 1
  return vals
end

function NLPModels.hess_coord!(
  nlp::AmplModel,
  x::AbstractVector,
  y::AbstractVector{<:AbstractFloat},
  vals::AbstractVector{<:AbstractFloat};
  kwargs...,
)
  vals_ = Vector{Cdouble}(undef, nlp.meta.nnzh)
  hess_coord!(nlp, x, Vector{Cdouble}(y), vals_; kwargs...)
  vals[1:(nlp.meta.nnzh)] .= vals_
  return vals
end

# evaluate the objective Hessian
function NLPModels.hess_coord!(
  nlp::AmplModel,
  x::AbstractVector,
  vals::Vector{Cdouble};
  obj_weight::Float64 = 1.0,
)
  # Note: x is in fact not used.
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  # if nlp.safe
  _ = obj(nlp, x)
  nlp.counters.neval_obj -= 1
  _ = cons(nlp, x)
  nlp.counters.neval_cons -= 1
  # end

  asl_hessval(nlp.__asl, C_NULL, obj_weight, vals)
  nlp.counters.neval_hess += 1
  return vals
end

function NLPModels.hess_coord!(
  nlp::AmplModel,
  x::AbstractVector,
  vals::AbstractVector{<:AbstractFloat};
  kwargs...,
)
  vals_ = Vector{Cdouble}(undef, nlp.meta.nnzh)
  hess_coord!(nlp, x, vals_; kwargs...)
  vals[1:(nlp.meta.nnzh)] .= vals_
  return vals
end
