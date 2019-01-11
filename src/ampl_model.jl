export AmplModel, AmplException,
       reset!,
       write_sol, amplmodel_finalize, varscale, lagscale, conscale,
       obj, grad, grad!,
       cons, cons!, jth_con, jth_congrad, jth_congrad!, jth_sparse_congrad,
       jac_coord, jac, jprod, jprod!, jtprod, jtprod!,
       jth_hprod, jth_hprod!, ghjvprod, ghjvprod!,
       hess_coord, hess, hprod, hprod!

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
import NLPModels.reset!
import NLPModels.varscale, NLPModels.lagscale, NLPModels.conscale
import NLPModels.obj, NLPModels.grad, NLPModels.grad!
import NLPModels.cons, NLPModels.cons!, NLPModels.jth_con
import NLPModels.jth_congrad, NLPModels.jth_congrad!, NLPModels.jth_sparse_congrad
import NLPModels.jac_coord, NLPModels.jac
import NLPModels.jprod, NLPModels.jtprod, NLPModels.jprod!, NLPModels.jtprod!
import NLPModels.jth_hprod, NLPModels.jth_hprod!
import NLPModels.ghjvprod, NLPModels.ghjvprod!
import NLPModels.hess_coord, NLPModels.hess, NLPModels.hprod, NLPModels.hprod!

# Methods associated to AmplModel instances.

"Reset evaluation counters in `nlp`."
function reset!(nlp :: AmplModel)
  reset!(nlp.counters)
  return nlp
end

"Write message `msg` along with primal and dual variables `x` and `y` to file."
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

"Scale the vector of variables by the vector `s`."
function varscale(nlp :: AmplModel, s :: AbstractVector)
  @check_ampl_model
  length(s) >= nlp.meta.nvar || error("s must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  @asl_call(:asl_varscale, Nothing, (Ptr{Nothing}, Ptr{Float64}, Ref{Cint}), nlp.__asl, s, err)
  err == 0 || throw(AmplException("Error while scaling variables"))
end

"""Set the scaling factor σ in the Lagrangian:
    L(x,y) = f(x) + σ ∑ yi ci(x).
"""
function lagscale(nlp :: AmplModel, σ :: Float64)
  @check_ampl_model
  err = Cint(0)
  @asl_call(:asl_lagscale, Nothing, (Ptr{Nothing}, Float64, Ref{Cint}), nlp.__asl, σ, err)
  err == 0 || throw(AmplException("Error while scaling Lagrangian"))
end

"Scale the vector of constraints by the vector `s`."
function conscale(nlp :: AmplModel, s :: AbstractVector)
  @check_ampl_model
  length(s) >= nlp.meta.ncon || error("s must have length at least $(nlp.meta.ncon)")

  err = Cint(0)
  @asl_call(:asl_conscale, Nothing, (Ptr{Nothing}, Ptr{Float64}, Ref{Cint}), nlp.__asl, s, err)
  err == 0 || throw(AmplException("Error while scaling constraints"))
end

# Evaluating objective, constraints and derivatives.

"Evaluate the objective function of `nlp` at `x`."
function obj(nlp :: AmplModel, x :: AbstractVector)
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  f = @asl_call(:asl_obj, Float64, (Ptr{Nothing}, Ptr{Float64}, Ref{Cint}), nlp.__asl, x, err)
  nlp.counters.neval_obj += 1
  err == 0 || throw(AmplException("Error while evaluating objective"))
  return f
end

"Evaluate the gradient of the objective function at `x`."
function grad(nlp :: AmplModel, x :: AbstractVector)
  g = Vector{Float64}(undef, nlp.meta.nvar)
  return grad!(nlp, x, g)
end

"Evaluate the gradient of the objective function at `x` in place."
function grad!(nlp :: AmplModel, x :: AbstractVector, g :: AbstractVector)
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  @asl_call(:asl_grad, Ptr{Float64},
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Ref{Cint}),
             nlp.__asl, x,            g,            err)
  nlp.counters.neval_grad += 1
  err == 0 || throw(AmplException("Error while evaluating objective gradient"))
  return g
end

"Evaluate the constraints at `x`."
function cons(nlp :: AmplModel, x :: AbstractVector)
  c = Vector{Float64}(undef, nlp.meta.ncon)
  return cons!(nlp, x, c)
end

"Evaluate the constraints at `x` in place."
function cons!(nlp :: AmplModel, x :: AbstractVector, c :: AbstractVector)
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  @asl_call(:asl_cons, Nothing,
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Ref{Cint}),
             nlp.__asl, x,            c,            err)
  nlp.counters.neval_cons += 1
  err == 0 || throw(AmplException("Error while evaluating constraints"))
  return c
end

"Evaluate the `j`-th constraint at `x`."
function jth_con(nlp :: AmplModel, x :: AbstractVector, j :: Int)
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  cj = @asl_call(:asl_jcon, Float64,
                 (Ptr{Nothing}, Ptr{Float64}, Int32, Ref{Cint}),
                  nlp.__asl, x,            j-1,   err)
  nlp.counters.neval_jcon += 1
  err == 0 || throw(AmplException("Error while evaluating $j-th constraint"))
  return cj
end

"Evaluate the `j`-th constraint gradient at `x`."
function jth_congrad(nlp :: AmplModel, x :: AbstractVector, j :: Int)
  g = Vector{Float64}(undef, nlp.meta.nvar)
  return jth_congrad!(nlp, x, j, g)
end

"Evaluate the `j`-th constraint gradient at `x` in place."
function jth_congrad!(nlp :: AmplModel, x :: AbstractVector, j :: Int, g :: AbstractVector)
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  @asl_call(:asl_jcongrad, Ptr{Float64},
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Int32, Ref{Cint}),
             nlp.__asl, x,            g,            j-1,   err)
  nlp.counters.neval_jgrad += 1
  err == 0 || throw(AmplException("Error while evaluating $j-th constraint gradient"))
  return g
end

"Evaluate the `j`-th constraint sparse gradient at `x`."
function jth_sparse_congrad(nlp :: AmplModel, x :: AbstractVector, j :: Int)
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
  @. inds += 1
  return sparsevec(inds, vals, nlp.meta.nvar)
end

"Evaluate the constraints Jacobian at `x` in sparse coordinate format."
function jac_coord(nlp :: AmplModel, x :: AbstractVector)
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint(0)
  rows = Vector{Cint}(undef, nlp.meta.nnzj)
  cols = Vector{Cint}(undef, nlp.meta.nnzj)
  vals = Vector{Float64}(undef, nlp.meta.nnzj)
  @asl_call(:asl_jac, Nothing,
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ref{Cint}),
             nlp.__asl, x,            rows,      cols,      vals,         err)
  nlp.counters.neval_jac += 1
  err == 0 || throw(AmplException("Error while evaluating constraints Jacobian"))
  # Use 1-based indexing.
  @. rows += 1
  @. cols += 1
  return (rows, cols, vals)
end

"Evaluate the constraints Jacobian at `x` as a sparse matrix."
function jac(nlp :: AmplModel, x :: AbstractVector)
  @check_ampl_model
  (rows, cols, vals) = jac_coord(nlp, x)
  return sparse(rows, cols, vals, nlp.meta.ncon, nlp.meta.nvar)
end

"""
Evaluate the Jacobian-vector product at `x`.
Warning: Currently building the Jacobian for this.
"""
function jprod(nlp :: AmplModel, x :: AbstractVector, v :: AbstractVector)
  Jv = zeros(nlp.meta.ncon)
  return jprod!(nlp, x, v, Jv)
end

"""
Evaluate the Jacobian-vector product at `x` in place.
Warning: Currently building the Jacobian for this.
"""
function jprod!(nlp :: AmplModel,
                x :: AbstractVector,
                v :: AbstractVector,
                Jv :: AbstractVector)
  nlp.counters.neval_jac -= 1
  nlp.counters.neval_jprod += 1
  Jv[1:nlp.meta.ncon] = jac(nlp, x) * v
  return Jv
end

"""
Evaluate the transposed-Jacobian-vector product at `x`.
Warning: Currently building the Jacobian for this.
"""
function jtprod(nlp :: AmplModel, x :: AbstractVector, v :: AbstractVector)
  Jtv = zeros(nlp.meta.nvar)
  return jtprod!(nlp, x, v, Jtv)
end

"""
Evaluate the transposed-Jacobian-vector product at `x` in place.
Warning: Currently building the Jacobian for this.
"""
function jtprod!(nlp :: AmplModel,
                 x :: AbstractVector,
                 v :: AbstractVector,
                 Jtv :: AbstractVector)
  nlp.counters.neval_jac -= 1
  nlp.counters.neval_jtprod += 1
  Jtv[1:nlp.meta.nvar] = jac(nlp, x)' * v
  return Jtv
end

"Evaluate the product of the Lagrangian Hessian at `(x,y)` with the vector `v`."
function hprod(nlp :: AmplModel,
               x :: AbstractVector,
               v :: AbstractVector;
               y :: AbstractVector = nlp.meta.y0,
               obj_weight :: Float64 = 1.0)
  hv = Vector{Float64}(undef, nlp.meta.nvar);
  return hprod!(nlp, x, v, hv, y=y, obj_weight=obj_weight)
end

"Evaluate the product of the Lagrangian Hessian at `(x,y)` with the vector `v` in place."
function hprod!(nlp :: AmplModel,
                x :: AbstractVector,
                v :: AbstractVector,
                hv :: AbstractVector;
                y :: AbstractVector = nlp.meta.y0,
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
             nlp.__asl, y,            v,            hv,           obj_weight);
  nlp.counters.neval_hprod += 1
  return hv
end

"""Evaluate the product of the `j`-th constraint Hessian at `x` with the vector `v`.
The objective Hessian is used if `j=0`.
"""
function jth_hprod(nlp :: AmplModel,
                   x :: AbstractVector, v :: AbstractVector, j :: Int)
  hv = Vector{Float64}(undef, nlp.meta.nvar)
  return jth_hprod!(nlp, x, v, j, hv)
end

"""Evaluate the product of the `j`-th constraint Hessian at `x` with the vector `v` in place.
The objective Hessian is used if `j=0`.
"""
function jth_hprod!(nlp :: AmplModel,
                    x :: AbstractVector, v :: AbstractVector,
                    j :: Int, hv :: AbstractVector)
# Note: x is in fact not used in hprod.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")

  if nlp.safe
    if j == 0
      _ = obj(nlp, x) ; nlp.counters.neval_obj -= 1
    else
      _ = cons(nlp, x) ; nlp.counters.neval_cons -= 1
    end
  end
  @asl_call(:asl_hvcompd, Ptr{Float64},
            (Ptr{Nothing}, Ptr{Float64}, Ptr{Float64}, Int),
             nlp.__asl, v,            hv,           j-1);
  nlp.counters.neval_jhprod += 1
  return hv
end

"""Compute the vector of dot products `(g, Hj*v)`
where `Hj` is the Hessian of the `j`-th constraint at `x`.
"""
function ghjvprod(nlp :: AmplModel,
                  x :: AbstractVector, g :: AbstractVector, v :: AbstractVector)
  gHv = Vector{Float64}(undef, nlp.meta.ncon);
  return ghjvprod!(nlp, x, g, v, gHv)
end

"""Compute the vector of dot products `(g, Hj*v)` in place
where `Hj` is the Hessian of the `j`-th constraint at `x`.
"""
function ghjvprod!(nlp :: AmplModel,
                   x :: AbstractVector, g :: AbstractVector,
                   v :: AbstractVector, gHv :: AbstractVector)
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
             nlp.__asl, g,            v,            gHv);
  nlp.counters.neval_hprod += nlp.meta.ncon
end

"""Evaluate the Lagrangian Hessian at `(x,y)` in sparse coordinate format.
Only the lower triangle is returned.
"""
function hess_coord(nlp :: AmplModel,
                    x :: AbstractVector;
                    y :: AbstractVector = nlp.meta.y0,
                    obj_weight :: Float64 = 1.0)
  # Note: x is in fact not used.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(y) >= nlp.meta.ncon || error("y must have length at least $(nlp.meta.ncon)")

  if nlp.safe
    _ = obj(nlp, x) ; nlp.counters.neval_obj -= 1
    _ = cons(nlp, x) ; nlp.counters.neval_cons -= 1
  end

  rows = Vector{Cint}(undef, nlp.meta.nnzh)
  cols = Vector{Cint}(undef, nlp.meta.nnzh)
  vals = Vector{Float64}(undef, nlp.meta.nnzh)
  @asl_call(:asl_hess, Nothing,
            (Ptr{Nothing}, Ptr{Float64}, Float64,    Ptr{Cint}, Ptr{Cint}, Ptr{Float64}),
             nlp.__asl, y,            obj_weight, rows,      cols,      vals)
  nlp.counters.neval_hess += 1
  # Use 1-based indexing.
  # Swap rows and cols to obtain the lower triangle.
  @. cols += 1
  @. rows += 1
  return (cols, rows, vals)
end

"""Evaluate the Lagrangian Hessian at `(x,y)` as a sparse matrix.
Only the lower triangle is returned.
"""
function hess(nlp :: AmplModel,
              x :: AbstractVector;
              y :: AbstractVector = nlp.meta.y0,
              obj_weight :: Float64 = 1.0)
  @check_ampl_model
  (rows, cols, vals) = hess_coord(nlp, x, y=y, obj_weight=obj_weight);
  return sparse(rows, cols, vals, nlp.meta.nvar, nlp.meta.nvar)
end
