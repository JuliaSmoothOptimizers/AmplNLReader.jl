export AmplModel, AmplException, write_sol, amplmodel_finalize

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

    asl = @ccall libasl.asl_init(stub::Ptr{UInt8})::Ptr{Cvoid}
    (asl == C_NULL) && error("Error allocating ASL structure")

    minimize = (@ccall libasl.asl_objtype(asl::Ptr{Cvoid})::Cint) == 0
    islp = (@ccall libasl.asl_islp(asl::Ptr{Cvoid})::Cint) != 0

    nlo = (@ccall libasl.asl_nlo(asl::Ptr{Cvoid})::Cint) |> Int

    nvar = (@ccall libasl.asl_nvar(asl::Ptr{Cvoid})::Cint) |> Int
    ncon = (@ccall libasl.asl_ncon(asl::Ptr{Cvoid})::Cint) |> Int

    x0 = unsafe_wrap(
      Array,
      (@ccall libasl.asl_x0(asl::Ptr{Cvoid})::Ptr{Cdouble}),
      (nvar,),
      own = false,
    )
    y0 = unsafe_wrap(
      Array,
      (@ccall libasl.asl_y0(asl::Ptr{Cvoid})::Ptr{Cdouble}),
      (ncon,),
      own = false,
    )

    lvar = unsafe_wrap(
      Array,
      (@ccall libasl.asl_lvar(asl::Ptr{Cvoid})::Ptr{Cdouble}),
      (nvar,),
      own = false,
    )
    uvar = unsafe_wrap(
      Array,
      (@ccall libasl.asl_uvar(asl::Ptr{Cvoid})::Ptr{Cdouble}),
      (nvar,),
      own = false,
    )

    nzo = (@ccall libasl.asl_nzo(asl::Ptr{Cvoid})::Cint) |> Int
    nbv = (@ccall libasl.asl_nbv(asl::Ptr{Cvoid})::Cint) |> Int
    niv = (@ccall libasl.asl_niv(asl::Ptr{Cvoid})::Cint) |> Int
    nlvb = (@ccall libasl.asl_nlvb(asl::Ptr{Cvoid})::Cint) |> Int
    nlvo = (@ccall libasl.asl_nlvo(asl::Ptr{Cvoid})::Cint) |> Int
    nlvc = (@ccall libasl.asl_nlvc(asl::Ptr{Cvoid})::Cint) |> Int
    nlvbi = (@ccall libasl.asl_nlvbi(asl::Ptr{Cvoid})::Cint) |> Int
    nlvci = (@ccall libasl.asl_nlvci(asl::Ptr{Cvoid})::Cint) |> Int
    nlvoi = (@ccall libasl.asl_nlvoi(asl::Ptr{Cvoid})::Cint) |> Int
    nwv = (@ccall libasl.asl_nwv(asl::Ptr{Cvoid})::Cint) |> Int
    n_cc = (@ccall libasl.asl_n_cc(asl::Ptr{Cvoid})::Cint) |> Int

    lcon = unsafe_wrap(
      Array,
      (@ccall libasl.asl_lcon(asl::Ptr{Cvoid})::Ptr{Cdouble}),
      (ncon,),
      own = false,
    )
    ucon = unsafe_wrap(
      Array,
      (@ccall libasl.asl_ucon(asl::Ptr{Cvoid})::Ptr{Cdouble}),
      (ncon,),
      own = false,
    )

    if n_cc > 0
      cvar = unsafe_wrap(
        Array,
        (@ccall libasl.asl_cvar(asl::Ptr{Cvoid})::Ptr{Cint}),
        (ncon,),
        own = false,
      )

      # Check complementarity constraints are well specified:
      cc_cons = cvar .> 0
      @assert all(isinf, ucon[cc_cons])
      @assert all(isfinite, lcon[cc_cons])
    else
      cvar = Int[]
    end

    nlnet = (@ccall libasl.asl_lnc(asl::Ptr{Cvoid})::Cint) |> Int
    nnnet = (@ccall libasl.asl_nlnc(asl::Ptr{Cvoid})::Cint) |> Int
    nnln = ((@ccall libasl.asl_nlc(asl::Ptr{Cvoid})::Cint) |> Int) - nnnet

    nln = 1:nnln
    nnet = (nnln + 1):(nnln + nnnet)
    lnet = (nnln + nnnet + 1):(nnln + nnnet + nlnet)
    lin = (nnln + nnnet + nlnet + 1):ncon

    nnzj = (@ccall libasl.asl_nnzj(asl::Ptr{Cvoid})::Cint) |> Int
    lin_nnzj = 0
    for j in lin
      lin_nnzj +=
        (@ccall libasl.asl_sparse_congrad_nnz(asl::Ptr{Cvoid}, (j - 1)::Cint)::Csize_t) |> Cint
    end
    nln_nnzj = 0
    for j in nln
      nln_nnzj +=
        (@ccall libasl.asl_sparse_congrad_nnz(asl::Ptr{Cvoid}, (j - 1)::Cint)::Csize_t) |> Cint
    end
    nnzh = (@ccall libasl.asl_nnzh(asl::Ptr{Cvoid})::Cint) |> Int

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

# Import methods we override.
import Base.show, Base.print

# Methods associated to AmplModel instances.

function NLPModels.reset!(nlp::AmplModel)
  reset!(nlp.counters)
  return nlp
end

function write_sol(nlp::AmplModel, msg::String, x::AbstractVector, y::AbstractVector)
  check_ampl_model(nlp)
  length(x) == nlp.meta.nvar || error("x must have length $(nlp.meta.nvar)")
  length(y) == nlp.meta.ncon || error("y must have length $(nlp.meta.ncon)")

  @ccall libasl.asl_write_sol(
    nlp.__asl::Ptr{Cvoid},
    msg::Ptr{UInt8},
    x::Ptr{Cdouble},
    y::Ptr{Cdouble},
  )::Cvoid
end

function amplmodel_finalize(nlp::AmplModel)
  if nlp.__asl != C_NULL
    @ccall libasl.asl_finalize(nlp.__asl::Ptr{Cvoid})::Cvoid
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
  @ccall libasl.asl_varscale(nlp.__asl::Ptr{Cvoid}, s::Ptr{Cdouble}, err::Ref{Cint})::Cvoid
  err[] == 0 || throw(AmplException("Error while scaling variables"))
end

NLPModels.varscale(nlp::AmplModel, s::AbstractVector) = varscale(nlp, Vector{Cdouble}(s))

function NLPModels.lagscale(nlp::AmplModel, σ::Float64)
  check_ampl_model(nlp)
  err = Ref{Cint}(0)
  @ccall libasl.asl_lagscale(nlp.__asl::Ptr{Cvoid}, σ::Cdouble, err::Ref{Cint})::Cvoid
  err[] == 0 || throw(AmplException("Error while scaling Lagrangian"))
end

function NLPModels.conscale(nlp::AmplModel, s::Vector{Cdouble})
  check_ampl_model(nlp)
  length(s) >= nlp.meta.ncon || error("s must have length at least $(nlp.meta.ncon)")

  err = Ref{Cint}(0)
  @ccall libasl.asl_conscale(nlp.__asl::Ptr{Cvoid}, s::Ptr{Cdouble}, err::Ref{Cint})::Cvoid
  err[] == 0 || throw(AmplException("Error while scaling constraints"))
end

NLPModels.conscale(nlp::AmplModel, s::AbstractVector) = conscale(nlp, Vector{Cdouble}(s))

# Evaluating objective, constraints and derivatives.

function NLPModels.obj(nlp::AmplModel, x::Vector{Cdouble})
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Ref{Cint}(0)
  f = @ccall libasl.asl_obj(nlp.__asl::Ptr{Cvoid}, x::Ptr{Cdouble}, err::Ref{Cint})::Float64
  nlp.counters.neval_obj += 1
  err[] == 0 || throw(AmplException("Error while evaluating objective"))
  return f
end

NLPModels.obj(nlp::AmplModel, x::AbstractVector) = obj(nlp, Vector{Cdouble}(x))

function NLPModels.grad!(nlp::AmplModel, x::Vector{Cdouble}, g::Vector{Cdouble})
  check_ampl_model(nlp)
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Ref{Cint}(0)
  @ccall libasl.asl_grad(
    nlp.__asl::Ptr{Cvoid},
    x::Ptr{Cdouble},
    g::Ptr{Cdouble},
    err::Ref{Cint},
  )::Cvoid
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
  @ccall libasl.asl_cons(
    nlp.__asl::Ptr{Cvoid},
    x::Ptr{Cdouble},
    c::Ptr{Cdouble},
    err::Ref{Cint},
  )::Cvoid
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
  cj = @ccall libasl.asl_jcon(
    nlp.__asl::Ptr{Cvoid},
    x::Ptr{Cdouble},
    (j - 1)::Cint,
    err::Ref{Cint},
  )::Float64
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
  @ccall libasl.asl_jcongrad(
    nlp.__asl::Ptr{Cvoid},
    x::Ptr{Cdouble},
    g::Ptr{Cdouble},
    (j - 1)::Cint,
    err::Ref{Cint},
  )::Cvoid
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

  nnz = @ccall libasl.asl_sparse_congrad_nnz(nlp.__asl::Ptr{Cvoid}, (j - 1)::Cint)::Csize_t

  err = Ref{Cint}(0)
  inds = Vector{Cint}(undef, nnz)
  vals = Vector{Cdouble}(undef, nnz)
  @ccall libasl.asl_sparse_congrad(
    nlp.__asl::Ptr{Cvoid},
    x::Ptr{Cdouble},
    (j - 1)::Cint,
    inds::Ptr{Cint},
    vals::Ptr{Cdouble},
    err::Ref{Cint},
  )::Cvoid
  nlp.counters.neval_jgrad += 1
  err[] == 0 || throw(AmplException("Error while evaluating $j-th sparse constraint gradient"))

  # Use 1-based indexing.
  inds .+= Cint(1)
  return sparsevec(inds, vals, nlp.meta.nvar)
end

NLPModels.jth_sparse_congrad(nlp::AmplModel, x::AbstractVector, j::Int) =
  jth_sparse_congrad(nlp, Vector{Cdouble}(x), j)

function NLPModels.jac_structure!(nlp::AmplModel, rows::Vector{Cint}, cols::Vector{Cint})
  @ccall libasl.asl_jac_structure(nlp.__asl::Ptr{Cvoid}, rows::Ptr{Cint}, cols::Ptr{Cint})::Cvoid

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
  @ccall libasl.asl_jacval(
    nlp.__asl::Ptr{Cvoid},
    x::Ptr{Cdouble},
    vals::Ptr{Cdouble},
    err::Ref{Cint},
  )::Cvoid
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
  v::Vector{Cdouble},
  hv::Vector{Cdouble};
  obj_weight::Float64 = 1.0,
)
  # Note: x is in fact not used in hprod.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")

  if nlp.safe
    _ = obj(nlp, x)
    nlp.counters.neval_obj -= 1
    _ = cons(nlp, x)
    nlp.counters.neval_cons -= 1
  end
  @asl_call(
    :asl_hprod,
    Nothing,
    (Ptr{Nothing}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble),
    nlp.__asl,
    v,
    hv,
    obj_weight
  )
  nlp.counters.neval_hprod += 1
  return hv
end

function NLPModels.hprod!(
  nlp::AmplModel,
  x::AbstractVector,
  v::AbstractVector,
  hv::AbstractVector;
  obj_weight::Float64 = 1.0,
)
  hv_ = Vector{Cdouble}(undef, nlp.meta.nvar)
  hprod!(nlp, x,Vector{Cdouble}(v), hv_; obj_weight = obj_weight)
  hv[1:(nlp.meta.nvar)] .= hv_
  return hv
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
  @ccall libasl.asl_hprod(
    nlp.__asl::Ptr{Cvoid},
    y::Ptr{Cdouble},
    v::Ptr{Cdouble},
    hv::Ptr{Cdouble},
    obj_weight::Cdouble,
  )::Cvoid
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
  @ccall libasl.asl_hvcompd(
    nlp.__asl::Ptr{Cvoid},
    v::Ptr{Cdouble},
    hv::Ptr{Cdouble},
    (j - 1)::Cint,
  )::Cvoid
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
  @ccall libasl.asl_ghjvprod(
    nlp.__asl::Ptr{Cvoid},
    g::Ptr{Cdouble},
    v::Ptr{Cdouble},
    gHv::Ptr{Cdouble},
  )::Cvoid
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
  @ccall libasl.asl_hess_structure(nlp.__asl::Ptr{Cvoid}, cols::Ptr{Cint}, rows::Ptr{Cint})::Cvoid

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

  @ccall libasl.asl_hessval(
    nlp.__asl::Ptr{Cvoid},
    y::Ptr{Cdouble},
    obj_weight::Cdouble,
    vals::Ptr{Cdouble},
  )::Cvoid
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

  @ccall libasl.asl_hessval(
    nlp.__asl::Ptr{Cvoid},
    C_NULL::Ptr{Cdouble},
    obj_weight::Cdouble,
    vals::Ptr{Cdouble},
  )::Cvoid
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
