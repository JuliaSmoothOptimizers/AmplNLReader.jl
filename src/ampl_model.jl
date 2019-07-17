export AmplModel, AmplException, write_sol, amplmodel_finalize

# Convenience macro.
macro asl_call(func, args...)
  args = map(esc,args)
  quote
    ccall(($(esc(func)), libasl), $(args...))
  end
end

struct AmplException
  msg :: String
end

macro check_ampl_model()
  esc(:(nlp.__asl == C_NULL && throw(AmplException("Uninitialized AMPL model"))))
end

mutable struct AmplModel <: AbstractNLPModel
  meta  :: NLPModelMeta;     # Problem metadata.
  __asl :: Ptr{Nothing};        # Pointer to internal ASL structure. Do not touch.

  counters :: Counters       # Evaluation counters
  safe :: Bool               # Always evaluate the objective before the Hessian.

  function AmplModel(stub :: AbstractString; safe :: Bool=false)

    # check that stub or stub.nl exists
    fname = basename(stub)
    ext = occursin(".", fname) ? split(fname, '.')[2] : ""
    if ext == "nl"
      isfile(stub) || throw(AmplException("cannot find $(stub)"))
    else
      isfile("$(stub).nl") || throw(AmplException("cannot find $(stub).nl"))
    end

    asl = @asl_call(:asl_init, Ptr{Nothing}, (Ptr{UInt8},), stub);
    asl == C_NULL && error("Error allocating ASL structure")

    minimize = @asl_call(:asl_objtype, Int32, (Ptr{Nothing},), asl) == 0;
    islp = @asl_call(:asl_islp, Int32, (Ptr{Nothing},), asl) != 0;

    nlo = Int(@asl_call(:asl_nlo, Int32, (Ptr{Nothing},), asl));

    nvar = Int(@asl_call(:asl_nvar, Int32, (Ptr{Nothing},), asl));
    ncon = Int(@asl_call(:asl_ncon, Int32, (Ptr{Nothing},), asl));

    x0   = unsafe_wrap(Array, @asl_call(:asl_x0,   Ptr{Float64}, (Ptr{Nothing},), asl),
                            (nvar,), own=false)
    y0   = unsafe_wrap(Array, @asl_call(:asl_y0,   Ptr{Float64}, (Ptr{Nothing},), asl),
                            (ncon,), own=false)

    lvar = unsafe_wrap(Array, @asl_call(:asl_lvar, Ptr{Float64}, (Ptr{Nothing},), asl),
                            (nvar,), own=false)
    uvar = unsafe_wrap(Array, @asl_call(:asl_uvar, Ptr{Float64}, (Ptr{Nothing},), asl),
                            (nvar,), own=false)

    nzo = Int(@asl_call(:asl_nzo, Int32, (Ptr{Nothing},), asl))
    nbv = Int(@asl_call(:asl_nbv, Int32, (Ptr{Nothing},), asl))
    niv = Int(@asl_call(:asl_niv, Int32, (Ptr{Nothing},), asl))
    nlvb = Int(@asl_call(:asl_nlvb, Int32, (Ptr{Nothing},), asl))
    nlvo = Int(@asl_call(:asl_nlvo, Int32, (Ptr{Nothing},), asl))
    nlvc = Int(@asl_call(:asl_nlvc, Int32, (Ptr{Nothing},), asl))
    nlvbi = Int(@asl_call(:asl_nlvbi, Int32, (Ptr{Nothing},), asl))
    nlvci = Int(@asl_call(:asl_nlvci, Int32, (Ptr{Nothing},), asl))
    nlvoi = Int(@asl_call(:asl_nlvoi, Int32, (Ptr{Nothing},), asl))
    nwv = Int(@asl_call(:asl_nwv, Int32, (Ptr{Nothing},), asl))

    lcon = unsafe_wrap(Array, @asl_call(:asl_lcon, Ptr{Float64}, (Ptr{Nothing},), asl),
                            (ncon,), own=false)
    ucon = unsafe_wrap(Array, @asl_call(:asl_ucon, Ptr{Float64}, (Ptr{Nothing},), asl),
                            (ncon,), own=false)

    nlnet = Int(@asl_call(:asl_lnc, Int32, (Ptr{Nothing},), asl))
    nnnet = Int(@asl_call(:asl_nlnc, Int32, (Ptr{Nothing},), asl))
    nnln = Int(@asl_call(:asl_nlc,  Int32, (Ptr{Nothing},), asl)) - nnnet
    nlin = ncon - nnln - nnnet

    nln  = 1 : nnln
    nnet = nnln+1 : nnln+nnnet
    lnet = nnln+nnnet+1 : nnln+nnnet+nlnet
    lin  = nnln+nnnet+nlnet+1 : ncon

    nnzj = Int(@asl_call(:asl_nnzj, Int32, (Ptr{Nothing},), asl))
    nnzh = Int(@asl_call(:asl_nnzh, Int32, (Ptr{Nothing},), asl))

    meta = NLPModelMeta(nvar, x0=x0, lvar=lvar, uvar=uvar,
                        nlo=nlo, nnzo=nzo,
                        ncon=ncon, y0=y0, lcon=lcon, ucon=ucon,
                        nnzj=nnzj, nnzh=nnzh,
                        nbv=nbv, niv=niv,
                        nlvb=nlvb, nlvo=nlvo, nlvc=nlvc,
                        nlvbi=nlvbi, nlvci=nlvci, nlvoi=nlvoi, nwv=nwv,
                        lin=lin, nln=nln, nnet=nnet, lnet=lnet,
                        nlin=nlin, nnln=nnln, nlnet=nlnet,
                        minimize=minimize, islp=islp, name=stub)

    nlp = new(meta, asl, Counters(), safe)

    finalizer(amplmodel_finalize, nlp)
    return nlp
  end

end


# Import methods we override.
import Base.show, Base.print

# Methods associated to AmplModel instances.

function NLPModels.reset!(nlp :: AmplModel)
  reset!(nlp.counters)
  return nlp
end

function write_sol(nlp :: AmplModel, msg :: String, x :: AbstractVector, y :: AbstractVector)
  @check_ampl_model
  length(x) == nlp.meta.nvar || error("x must have length $(nlp.meta.nvar)")
  length(y) == nlp.meta.ncon || error("y must have length $(nlp.meta.ncon)")

  @asl_call(:asl_write_sol, Nothing,
                    (Ptr{Nothing}, Ptr{UInt8}, Ptr{Float64}, Ptr{Float64}),
                     nlp.__asl, msg,        x,            y)
end

function amplmodel_finalize(nlp :: AmplModel)
  if nlp.__asl == C_NULL
    return
  end
  @asl_call(:asl_finalize, Nothing, (Ptr{Nothing},), nlp.__asl)
  nlp.__asl = C_NULL
end

# Displaying AmplModel instances.

function show(io :: IO, nlp :: AmplModel)
  @check_ampl_model
  show(io, nlp.meta)
end

function print(io :: IO, nlp :: AmplModel)
  @check_ampl_model
  print(io, nlp.meta)
end


# Scaling AmplModel instances.

function NLPModels.varscale(nlp :: AmplModel, s :: Vector{Float64})
  @check_ampl_model
  length(s) >= nlp.meta.nvar || error("s must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  @asl_call(:asl_varscale, Nothing, (Ptr{Nothing}, Ptr{Float64}, Ref{Cint}), nlp.__asl, s, err)
  err == 0 || throw(AmplException("Error while scaling variables"))
end

NLPModels.varscale(nlp :: AmplModel, s :: AbstractVector) = varscale(nlp, Vector{Float64}(s))

function NLPModels.lagscale(nlp :: AmplModel, σ :: Float64)
  @check_ampl_model
  err = Cint(0)
  @asl_call(:asl_lagscale, Nothing, (Ptr{Nothing}, Float64, Ref{Cint}), nlp.__asl, σ, err)
  err == 0 || throw(AmplException("Error while scaling Lagrangian"))
end

function NLPModels.conscale(nlp :: AmplModel, s :: Vector{Float64})
  @check_ampl_model
  length(s) >= nlp.meta.ncon || error("s must have length at least $(nlp.meta.ncon)")

  err = Cint(0)
  @asl_call(:asl_conscale, Nothing, (Ptr{Nothing}, Ptr{Float64}, Ref{Cint}), nlp.__asl, s, err)
  err == 0 || throw(AmplException("Error while scaling constraints"))
end

NLPModels.conscale(nlp :: AmplModel, s :: AbstractVector) = conscale(nlp, Vector{Float64}(s))

# Evaluating objective, constraints and derivatives.

function NLPModels.obj(nlp :: AmplModel, x :: Vector{Float64})
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  f = @asl_call(:asl_obj, Float64, (Ptr{Nothing}, Ptr{Float64}, Ref{Cint}), nlp.__asl, x, err)
  nlp.counters.neval_obj += 1
  err == 0 || throw(AmplException("Error while evaluating objective"))
  return f
end

NLPModels.obj(nlp :: AmplModel, x :: AbstractVector) = obj(nlp, Vector{Float64}(x))

function NLPModels.grad!(nlp :: AmplModel, x :: Vector{Float64}, g :: Vector{Float64})
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  @asl_call(:asl_grad, Ptr{Float64},
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Ref{Cint}),
             nlp.__asl,    x,            g,            err)
  nlp.counters.neval_grad += 1
  err == 0 || throw(AmplException("Error while evaluating objective gradient"))
  return g
end

function NLPModels.grad!(nlp :: AmplModel, x :: AbstractVector, g :: AbstractVector)
  g_ = Vector{Float64}(undef, nlp.meta.nvar)
  grad!(nlp, Vector{Float64}(x), g_)
  g[1 : nlp.meta.nvar] .= g_
  return g
end

function NLPModels.cons!(nlp :: AmplModel, x :: Vector{Float64}, c :: Vector{Float64})
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  @asl_call(:asl_cons, Nothing,
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Ref{Cint}),
             nlp.__asl,    x,            c,            err)
  nlp.counters.neval_cons += 1
  err == 0 || throw(AmplException("Error while evaluating constraints"))
  return c
end

function NLPModels.cons!(nlp :: AmplModel, x :: AbstractVector, c :: AbstractVector)
  c_ = Vector{Float64}(undef, nlp.meta.ncon)
  cons!(nlp, Vector{Float64}(x), c_)
  c[1 : nlp.meta.ncon] .= c_
  return c
end

function NLPModels.jth_con(nlp :: AmplModel, x :: Vector{Float64}, j :: Int)
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  cj = @asl_call(:asl_jcon, Float64,
                 (Ptr{Nothing}, Ptr{Float64}, Int32, Ref{Cint}),
                  nlp.__asl,    x,            j-1,   err)
  nlp.counters.neval_jcon += 1
  err == 0 || throw(AmplException("Error while evaluating $j-th constraint"))
  return cj
end

NLPModels.jth_con(nlp :: AmplModel, x :: AbstractVector, j :: Int) = jth_con(nlp, Vector{Float64}(x), j)

function NLPModels.jth_congrad!(nlp :: AmplModel, x :: Vector{Float64}, j :: Int, g :: Vector{Float64})
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  @asl_call(:asl_jcongrad, Ptr{Float64},
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Int32, Ref{Cint}),
             nlp.__asl,    x,            g,            j-1,   err)
  nlp.counters.neval_jgrad += 1
  err == 0 || throw(AmplException("Error while evaluating $j-th constraint gradient"))
  return g
end

function NLPModels.jth_congrad!(nlp :: AmplModel, x :: AbstractVector, j :: Int, g :: AbstractVector)
  g_ = Vector{Float64}(undef, nlp.meta.nvar)
  jth_congrad!(nlp, Vector{Float64}(x), j, g_)
  g[1 : nlp.meta.nvar] .= g_
  return g
end

function NLPModels.jth_sparse_congrad(nlp :: AmplModel, x :: Vector{Float64}, j :: Int)
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  nnz = @asl_call(:asl_sparse_congrad_nnz, Csize_t,
                  (Ptr{Nothing}, Cint), nlp.__asl, j-1)

  err = Cint(0)
  inds = Vector{Cint}(undef, nnz)
  vals = Vector{Float64}(undef, nnz)
  @asl_call(:asl_sparse_congrad, Nothing,
            (Ptr{Nothing}, Ptr{Float64}, Int32, Ptr{Cint}, Ptr{Float64}, Ref{Cint}),
             nlp.__asl, x,            j-1,   inds,      vals,         err)
  nlp.counters.neval_jgrad += 1
  err == 0 || throw(AmplException("Error while evaluating $j-th sparse constraint gradient"))
  # Use 1-based indexing.
  @. inds += Cint(1)
  return sparsevec(inds, vals, nlp.meta.nvar)
end

NLPModels.jth_sparse_congrad(nlp :: AmplModel, x :: AbstractVector, j :: Int) = jth_sparse_congrad(nlp, Vector{Float64}(x), j)

function NLPModels.jac_structure!(nlp :: AmplModel, rows :: Vector{Cint}, cols :: Vector{Cint})
  @asl_call(:asl_jac_structure, Nothing,
            (Ptr{Nothing}, Ptr{Int32}, Ptr{Int32}),
             nlp.__asl,    rows,      cols)
  # Use 1-based indexing.
  @. rows[1 : nlp.meta.nnzj] += Cint(1)
  @. cols[1 : nlp.meta.nnzj] += Cint(1)
  return rows, cols
end

function NLPModels.jac_structure!(nlp :: AmplModel, rows :: AbstractVector{<: Integer}, cols :: AbstractVector{<: Integer})
  rows_ = Vector{Cint}(undef, nlp.meta.nnzj)
  cols_ = Vector{Cint}(undef, nlp.meta.nnzj)
  jac_structure!(nlp, rows_, cols_)
  rows[1 : nlp.meta.nnzj] .= rows_
  cols[1 : nlp.meta.nnzj] .= cols_
  return rows, cols
end

function NLPModels.jac_coord!(nlp :: AmplModel, x :: Vector{Float64}, rows :: AbstractVector{<: Integer}, cols :: AbstractVector{<: Integer}, vals :: Vector{Cdouble})
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  _ = cons(nlp, x) ; nlp.counters.neval_cons -= 1

  err = Cint(0)
  @asl_call(:asl_jacval, Nothing,
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Cdouble}, Ref{Cint}),
             nlp.__asl,    x,            vals,         err)
  nlp.counters.neval_jac += 1
  err == 0 || throw(AmplException("Error while evaluating constraints Jacobian"))
  return (rows, cols, vals)
end

function NLPModels.jac_coord!(nlp :: AmplModel,
                              x :: AbstractVector,
                              rows :: AbstractVector{<: Integer},
                              cols :: AbstractVector{<: Integer},
                              vals :: AbstractVector{<: AbstractFloat})
  vals_ = Vector{Cdouble}(undef, nlp.meta.nnzj)
  jac_coord!(nlp, Vector{Float64}(x), rows, cols, vals_)
  vals[1 : nlp.meta.nnzj] .= vals_
  return (rows, cols, vals)
end

function NLPModels.jprod!(nlp :: AmplModel,
                          x :: AbstractVector,
                          v :: AbstractVector,
                          Jv :: AbstractVector)
  nlp.counters.neval_jac -= 1
  nlp.counters.neval_jprod += 1
  Jv[1:nlp.meta.ncon] = jac(nlp, Vector{Float64}(x)) * v
  return Jv
end

function NLPModels.jtprod!(nlp :: AmplModel,
                           x :: AbstractVector,
                           v :: AbstractVector,
                           Jtv :: AbstractVector)
  nlp.counters.neval_jac -= 1
  nlp.counters.neval_jtprod += 1
  Jtv[1:nlp.meta.nvar] = jac(nlp, Vector{Float64}(x))' * v
  return Jtv
end

function NLPModels.hprod!(nlp :: AmplModel,
                          x :: AbstractVector,
                          v :: Vector{Float64},
                          hv :: Vector{Float64};
                          y :: Vector{Float64} = nlp.meta.y0,
                          obj_weight :: Float64 = 1.0)
  # Note: x is in fact not used in hprod.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(y) >= nlp.meta.ncon || error("y must have length at least $(nlp.meta.ncon)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")

  if nlp.safe
    _ = obj(nlp, x) ; nlp.counters.neval_obj -= 1
    _ = cons(nlp, x) ; nlp.counters.neval_cons -= 1
  end
  @asl_call(:asl_hprod, Ptr{Float64},
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
             nlp.__asl,    y,            v,            hv,           obj_weight);
  nlp.counters.neval_hprod += 1
  return hv
end

function NLPModels.hprod!(nlp :: AmplModel,
                          x :: AbstractVector,
                          v :: AbstractVector,
                          hv :: AbstractVector;
                          y :: AbstractVector = nlp.meta.y0,
                          obj_weight :: Float64=1.0)
  hv_ = Vector{Float64}(undef, nlp.meta.nvar)
  hprod!(nlp, x, Vector{Float64}(v), hv_; y=Vector{Float64}(y), obj_weight=obj_weight)
  hv[1 : nlp.meta.nvar] .= hv_
  return hv
end

function NLPModels.jth_hprod!(nlp :: AmplModel,
                              x :: AbstractVector,
                              v :: Vector{Float64},
                              j :: Int,
                              hv :: Vector{Float64})
  # Note: x is in fact not used in hprod.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")

  if nlp.safe
    if j == 0
      _ = obj(nlp, Vector{Float64}(x)) ; nlp.counters.neval_obj -= 1
    else
      _ = cons(nlp, Vector{Float64}(x)) ; nlp.counters.neval_cons -= 1
    end
  end
  @asl_call(:asl_hvcompd, Ptr{Float64},
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Int),
             nlp.__asl,    v,            hv,           j-1);
  nlp.counters.neval_jhprod += 1
  return hv
end

function NLPModels.jth_hprod!(nlp :: AmplModel,
                              x :: AbstractVector,
                              v :: AbstractVector,
                              j :: Int,
                              hv :: AbstractVector)
  hv_ = Vector{Float64}(undef, nlp.meta.nvar)
  jth_hprod!(nlp, x, Vector{Float64}(v), j, hv_)
  hv[1 : nlp.meta.nvar] .= hv_
  return hv
end

function NLPModels.ghjvprod!(nlp :: AmplModel,
                             x :: AbstractVector,
                             g :: Vector{Float64},
                             v :: Vector{Float64},
                             gHv :: Vector{Float64})
  # Note: x is in fact not used.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(g) >= nlp.meta.nvar || error("g must have length at least $(nlp.meta.nvar)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")

  if nlp.safe
    _ = cons(nlp, x) ; nlp.counters.neval_cons -= 1
  end
  @asl_call(:asl_ghjvprod, Ptr{Float64},
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
             nlp.__asl,    g,            v,            gHv);
  nlp.counters.neval_hprod += nlp.meta.ncon
end

function NLPModels.ghjvprod!(nlp :: AmplModel,
                             x :: AbstractVector,
                             g :: AbstractVector,
                             v :: AbstractVector,
                             gHv :: AbstractVector)
  gHv_ = Vector{Float64}(undef, nlp.meta.nvar)
  ghjvprod!(nlp, x, Vector{Float64}(x), Vector{Float64}(g), Vector{Float64}(v), gHv_)
  gHv[1 : nlp.meta.var] .= gHv_
  return gHv
end

function NLPModels.hess_structure!(nlp :: AmplModel, rows :: Vector{Cint}, cols :: Vector{Cint})
  # Swap rows and cols to obtain the lower triangle.
  @asl_call(:asl_hess_structure, Nothing,
            (Ptr{Nothing}, Ptr{Int32}, Ptr{Int32}),
             nlp.__asl,    cols,      rows)
  # Use 1-based indexing.
  @. cols += Cint(1)
  @. rows += Cint(1)
  return (rows, cols)
end

function NLPModels.hess_structure!(nlp :: AmplModel, rows :: AbstractVector{<: Integer}, cols :: AbstractVector{<: Integer})
  rows_ = Vector{Cint}(undef, nlp.meta.nnzh)
  cols_ = Vector{Cint}(undef, nlp.meta.nnzh)
  hess_structure!(nlp, rows_, cols_)
  rows[1 : nlp.meta.nnzh] .= rows_
  cols[1 : nlp.meta.nnzh] .= cols_
  return (rows, cols)
end

function NLPModels.hess_coord!(nlp :: AmplModel,
                               x :: AbstractVector,
                               rows :: AbstractVector{<: Integer},
                               cols :: AbstractVector{<: Integer},
                               vals :: Vector{Cdouble};
                               y :: Vector{Float64} = nlp.meta.y0,
                               obj_weight :: Float64 = 1.0)
  # Note: x is in fact not used.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(y) >= nlp.meta.ncon || error("y must have length at least $(nlp.meta.ncon)")

  # if nlp.safe
    _ = obj(nlp, x) ; nlp.counters.neval_obj -= 1
    _ = cons(nlp, x) ; nlp.counters.neval_cons -= 1
  # end

  @asl_call(:asl_hessval, Nothing,
            (Ptr{Nothing}, Ptr{Float64}, Float64,    Ptr{Cdouble}),
             nlp.__asl,    y,            obj_weight, vals)
  nlp.counters.neval_hess += 1
  return (rows, cols, vals)
end

function NLPModels.hess_coord!(nlp :: AmplModel,
                               x :: AbstractVector,
                               rows :: AbstractVector{<: Integer},
                               cols :: AbstractVector{<: Integer},
                               vals :: AbstractVector{<: AbstractFloat}; kwargs...)
  vals_ = Vector{Cdouble}(undef, nlp.meta.nnzh)
  hess_coord!(nlp, x, rows, cols, vals_; kwargs...)
  vals[1 : nlp.meta.nnzh] .= vals_
  return (rows, cols, vals)
end
