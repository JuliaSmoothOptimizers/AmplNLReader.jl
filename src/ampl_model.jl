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

type AmplException
  msg :: String
end

macro check_ampl_model()
  esc(:(nlp.__asl == C_NULL && throw(AmplException("Uninitialized AMPL model"))))
end

type AmplModel <: AbstractNLPModel
  meta  :: NLPModelMeta;     # Problem metadata.
  __asl :: Ptr{Void};        # Pointer to internal ASL structure. Do not touch.

  counters :: Counters       # Evaluation counters
  safe :: Bool               # Always evaluate the objective before the Hessian.

  function AmplModel(stub :: AbstractString; safe :: Bool=false)
    asl = @compat @asl_call(:asl_init, Ptr{Void}, (Ptr{UInt8},), stub);
    asl == C_NULL && error("Error allocating ASL structure")

    minimize = @asl_call(:asl_objtype, Int32, (Ptr{Void},), asl) == 0;
    islp = @asl_call(:asl_islp, Int32, (Ptr{Void},), asl) != 0;

    nlo = @compat Int(@asl_call(:asl_nlo, Int32, (Ptr{Void},), asl));

    nvar = @compat Int(@asl_call(:asl_nvar, Int32, (Ptr{Void},), asl));
    ncon = @compat Int(@asl_call(:asl_ncon, Int32, (Ptr{Void},), asl));

    x0   = @compat unsafe_wrap(Array, @asl_call(:asl_x0,   Ptr{Float64}, (Ptr{Void},), asl),
                            (nvar,), false)
    y0   = @compat unsafe_wrap(Array, @asl_call(:asl_y0,   Ptr{Float64}, (Ptr{Void},), asl),
                            (ncon,), false)

    lvar = @compat unsafe_wrap(Array, @asl_call(:asl_lvar, Ptr{Float64}, (Ptr{Void},), asl),
                            (nvar,), false)
    uvar = @compat unsafe_wrap(Array, @asl_call(:asl_uvar, Ptr{Float64}, (Ptr{Void},), asl),
                            (nvar,), false)

    nzo = @compat Int(@asl_call(:asl_nzo, Int32, (Ptr{Void},), asl))
    nbv = @compat Int(@asl_call(:asl_nbv, Int32, (Ptr{Void},), asl))
    niv = @compat Int(@asl_call(:asl_niv, Int32, (Ptr{Void},), asl))
    nlvb = @compat Int(@asl_call(:asl_nlvb, Int32, (Ptr{Void},), asl))
    nlvo = @compat Int(@asl_call(:asl_nlvo, Int32, (Ptr{Void},), asl))
    nlvc = @compat Int(@asl_call(:asl_nlvc, Int32, (Ptr{Void},), asl))
    nlvbi = @compat Int(@asl_call(:asl_nlvbi, Int32, (Ptr{Void},), asl))
    nlvci = @compat Int(@asl_call(:asl_nlvci, Int32, (Ptr{Void},), asl))
    nlvoi = @compat Int(@asl_call(:asl_nlvoi, Int32, (Ptr{Void},), asl))
    nwv = @compat Int(@asl_call(:asl_nwv, Int32, (Ptr{Void},), asl))

    lcon = @compat unsafe_wrap(Array, @asl_call(:asl_lcon, Ptr{Float64}, (Ptr{Void},), asl),
                            (ncon,), false)
    ucon = @compat unsafe_wrap(Array, @asl_call(:asl_ucon, Ptr{Float64}, (Ptr{Void},), asl),
                            (ncon,), false)

    nlnet = @compat Int(@asl_call(:asl_lnc, Int32, (Ptr{Void},), asl))
    nnnet = @compat Int(@asl_call(:asl_nlnc, Int32, (Ptr{Void},), asl))
    nnln = @compat(Int(@asl_call(:asl_nlc,  Int32, (Ptr{Void},), asl))) - nnnet
    nlin = ncon - nnln - nnnet

    nln  = 1 : nnln
    nnet = nnln+1 : nnln+nnnet
    lnet = nnln+nnnet+1 : nnln+nnnet+nlnet
    lin  = nnln+nnnet+nlnet+1 : ncon

    nnzj = @compat Int(@asl_call(:asl_nnzj, Int32, (Ptr{Void},), asl))
    nnzh = @compat Int(@asl_call(:asl_nnzh, Int32, (Ptr{Void},), asl))

    meta = NLPModelMeta(nvar, x0=x0, lvar=lvar, uvar=uvar,
                        nlo=nlo, nnzo=nzo,
                        ncon=ncon, y0=y0, lcon=lcon, ucon=ucon,
                        nnzj=nnzj, nnzh=nnzh,
                        nbv=nbv, niv=niv,
                        nlvb=nlvb, nlvo=nlvo, nlvc=nlvc,
                        nlvbi=nlvbi, nlvci=nlvci, nlvoi=nlvoi, nwv=nwv,
                        lin=lin, nln=nln, nnet=nnet, lnet=lnet,
                        nlin=nlin, nnln=nnln, nnet=nnet, nlnet=nlnet,
                        minimize=minimize, islp=islp, name=stub)

    nlp = new(meta, asl, Counters(), safe)

    finalizer(nlp, amplmodel_finalize)
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
import NLPModels.NLPtoMPB

# Methods associated to AmplModel instances.

"Reset evaluation counters in `nlp`."
function reset!(nlp :: AmplModel)
  reset!(nlp.counters)
  return nlp
end

"Write message `msg` along with primal and dual variables `x` and `y` to file."
function write_sol(nlp :: AmplModel, msg :: String, x :: Vector{Float64}, y :: Vector{Float64})
  @check_ampl_model
  length(x) == nlp.meta.nvar || error("x must have length $(nlp.meta.nvar)")
  length(y) == nlp.meta.ncon || error("y must have length $(nlp.meta.ncon)")

  @compat @asl_call(:asl_write_sol, Void,
                    (Ptr{Void}, Ptr{UInt8}, Ptr{Float64}, Ptr{Float64}),
                     nlp.__asl, msg,        x,            y)
end

function amplmodel_finalize(nlp :: AmplModel)
  if nlp.__asl == C_NULL
    return
  end
  @asl_call(:asl_finalize, Void, (Ptr{Void},), nlp.__asl)
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
function varscale(nlp :: AmplModel, s :: Vector{Float64})
  @check_ampl_model
  length(s) >= nlp.meta.nvar || error("s must have length at least $(nlp.meta.nvar)")

  err = Cint[0]
  @asl_call(:asl_varscale, Void, (Ptr{Void}, Ptr{Float64}, Ptr{Cint}), nlp.__asl, s, err)
  err[1] == 0 || throw(AmplException("Error while scaling variables"))
end

"""Set the scaling factor σ in the Lagrangian:
    L(x,y) = f(x) + σ ∑ yi ci(x).
"""
function lagscale(nlp :: AmplModel, σ :: Float64)
  @check_ampl_model
  err = Cint[0]
  @asl_call(:asl_lagscale, Void, (Ptr{Void}, Float64, Ptr{Cint}), nlp.__asl, σ, err)
  err[1] == 0 || throw(AmplException("Error while scaling Lagrangian"))
end

"Scale the vector of constraints by the vector `s`."
function conscale(nlp :: AmplModel, s :: Vector{Float64})
  @check_ampl_model
  length(s) >= nlp.meta.ncon || error("s must have length at least $(nlp.meta.ncon)")

  err = Cint[0]
  @asl_call(:asl_conscale, Void, (Ptr{Void}, Ptr{Float64}, Ptr{Cint}), nlp.__asl, s, err)
  err[1] == 0 || throw(AmplException("Error while scaling constraints"))
end

# Evaluating objective, constraints and derivatives.

"Evaluate the objective function of `nlp` at `x`."
function obj(nlp :: AmplModel, x :: Vector{Float64})
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint[0]
  f = @asl_call(:asl_obj, Float64, (Ptr{Void}, Ptr{Float64}, Ptr{Cint}), nlp.__asl, x, err)
  nlp.counters.neval_obj += 1
  err[1] == 0 || throw(AmplException("Error while evaluating objective"))
  return f
end

"Evaluate the gradient of the objective function at `x`."
function grad(nlp :: AmplModel, x :: Vector{Float64})
  g = Array{Float64}(nlp.meta.nvar)
  return grad!(nlp, x, g)
end

"Evaluate the gradient of the objective function at `x` in place."
function grad!(nlp :: AmplModel, x :: Vector{Float64}, g :: Vector{Float64})
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint[0]
  @asl_call(:asl_grad, Ptr{Float64},
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}),
             nlp.__asl, x,            g,            err)
  nlp.counters.neval_grad += 1
  err[1] == 0 || throw(AmplException("Error while evaluating objective gradient"))
  return g
end

"Evaluate the constraints at `x`."
function cons(nlp :: AmplModel, x :: Vector{Float64})
  c = Array{Float64}(nlp.meta.ncon)
  return cons!(nlp, x, c)
end

"Evaluate the constraints at `x` in place."
function cons!(nlp :: AmplModel, x :: Vector{Float64}, c :: Vector{Float64})
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint[0]
  @asl_call(:asl_cons, Void,
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}),
             nlp.__asl, x,            c,            err)
  nlp.counters.neval_cons += 1
  err[1] == 0 || throw(AmplException("Error while evaluating constraints"))
  return c
end

"Evaluate the `j`-th constraint at `x`."
function jth_con(nlp :: AmplModel, x :: Vector{Float64}, j :: Int)
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint[0]
  cj = @asl_call(:asl_jcon, Float64,
                 (Ptr{Void}, Ptr{Float64}, Int32, Ptr{Cint}),
                  nlp.__asl, x,            j-1,   err)
  nlp.counters.neval_jcon += 1
  err[1] == 0 || throw(AmplException("Error while evaluating $j-th constraint"))
  return cj
end

"Evaluate the `j`-th constraint gradient at `x`."
function jth_congrad(nlp :: AmplModel, x :: Vector{Float64}, j :: Int)
  g = Array{Float64}(nlp.meta.nvar)
  return jth_congrad!(nlp, x, j, g)
end

"Evaluate the `j`-th constraint gradient at `x` in place."
function jth_congrad!(nlp :: AmplModel, x :: Vector{Float64}, j :: Int, g :: Vector{Float64})
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint[0]
  @asl_call(:asl_jcongrad, Ptr{Float64},
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Int32, Ptr{Cint}),
             nlp.__asl, x,            g,            j-1,   err)
  nlp.counters.neval_jgrad += 1
  err[1] == 0 || throw(AmplException("Error while evaluating $j-th constraint gradient"))
  return g
end

"Evaluate the `j`-th constraint sparse gradient at `x`."
function jth_sparse_congrad(nlp :: AmplModel, x :: Vector{Float64}, j :: Int)
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  nnz = @asl_call(:asl_sparse_congrad_nnz, Csize_t,
                  (Ptr{Void}, Cint), nlp.__asl, j-1)

  err = Cint[0]
  inds = Array{Int64}(nnz)
  vals = Array{Float64}(nnz)
  @asl_call(:asl_sparse_congrad, Void,
            (Ptr{Void}, Ptr{Float64}, Int32, Ptr{Int64}, Ptr{Float64}, Ptr{Cint}),
             nlp.__asl, x,            j-1,   inds,       vals,         err)
  nlp.counters.neval_jgrad += 1
  err[1] == 0 || throw(AmplException("Error while evaluating $j-th sparse constraint gradient"))
  # Use 1-based indexing.
  return sparsevec(inds+1, vals, nlp.meta.nvar)
end

"Evaluate the constraints Jacobian at `x` in sparse coordinate format."
function jac_coord(nlp :: AmplModel, x :: Vector{Float64})
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  err = Cint[0]
  rows = Array{Int64}(nlp.meta.nnzj)
  cols = Array{Int64}(nlp.meta.nnzj)
  vals = Array{Float64}(nlp.meta.nnzj)
  @asl_call(:asl_jac, Void,
            (Ptr{Void}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Cint}),
             nlp.__asl, x,            rows,       cols,       vals,         err)
  nlp.counters.neval_jac += 1
  err[1] == 0 || throw(AmplException("Error while evaluating constraints Jacobian"))
  # Use 1-based indexing.
  return (rows+1, cols+1, vals)
end

"Evaluate the constraints Jacobian at `x` as a sparse matrix."
function jac(nlp :: AmplModel, x :: Vector{Float64})
  @check_ampl_model
  (rows, cols, vals) = jac_coord(nlp, x)
  return sparse(rows, cols, vals, nlp.meta.ncon, nlp.meta.nvar)
end

"""
Evaluate the Jacobian-vector product at `x`.
Warning: Currently building the Jacobian for this.
"""
function jprod(nlp :: AmplModel, x :: Vector{Float64}, v :: Vector{Float64})
  Jv = zeros(nlp.meta.ncon)
  return jprod!(nlp, x, v, Jv)
end

"""
Evaluate the Jacobian-vector product at `x` in place.
Warning: Currently building the Jacobian for this.
"""
function jprod!(nlp :: AmplModel,
                x :: Vector{Float64},
                v :: Vector{Float64},
                Jv :: Vector{Float64})
  nlp.counters.neval_jac -= 1
  nlp.counters.neval_jprod += 1
  Jv[1:nlp.meta.ncon] = jac(nlp, x) * v
  return Jv
end

"""
Evaluate the transposed-Jacobian-vector product at `x`.
Warning: Currently building the Jacobian for this.
"""
function jtprod(nlp :: AmplModel, x :: Vector{Float64}, v :: Vector{Float64})
  Jtv = zeros(nlp.meta.nvar)
  return jtprod!(nlp, x, v, Jtv)
end

"""
Evaluate the transposed-Jacobian-vector product at `x` in place.
Warning: Currently building the Jacobian for this.
"""
function jtprod!(nlp :: AmplModel,
                 x :: Vector{Float64},
                 v :: Vector{Float64},
                 Jtv :: Vector{Float64})
  nlp.counters.neval_jac -= 1
  nlp.counters.neval_jtprod += 1
  Jtv[1:nlp.meta.nvar] = jac(nlp, x)' * v
  return Jtv
end

"Evaluate the product of the Lagrangian Hessian at `(x,y)` with the vector `v`."
function hprod(nlp :: AmplModel,
               x :: Vector{Float64},
               v :: Vector{Float64};
               y :: Vector{Float64} = nlp.meta.y0,
               obj_weight :: Float64 = 1.0)
  hv = Array{Float64}(nlp.meta.nvar);
  return hprod!(nlp, x, v, hv, y=y, obj_weight=obj_weight)
end

"Evaluate the product of the Lagrangian Hessian at `(x,y)` with the vector `v` in place."
function hprod!(nlp :: AmplModel,
                x :: Vector{Float64},
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
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
             nlp.__asl, y,            v,            hv,           obj_weight);
  nlp.counters.neval_hprod += 1
  return hv
end

"""Evaluate the product of the `j`-th constraint Hessian at `x` with the vector `v`.
The objective Hessian is used if `j=0`.
"""
function jth_hprod(nlp :: AmplModel,
                   x :: Vector{Float64}, v :: Vector{Float64}, j :: Int)
  hv = Array{Float64}(nlp.meta.nvar)
  return jth_hprod!(nlp, x, v, j, hv)
end

"""Evaluate the product of the `j`-th constraint Hessian at `x` with the vector `v` in place.
The objective Hessian is used if `j=0`.
"""
function jth_hprod!(nlp :: AmplModel,
                    x :: Vector{Float64}, v :: Vector{Float64},
                    j :: Int, hv :: Vector{Float64})
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
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Int),
             nlp.__asl, v,            hv,           j-1);
  nlp.counters.neval_jhprod += 1
  return hv
end

"""Compute the vector of dot products `(g, Hj*v)`
where `Hj` is the Hessian of the `j`-th constraint at `x`.
"""
function ghjvprod(nlp :: AmplModel,
                  x :: Vector{Float64}, g :: Vector{Float64}, v :: Vector{Float64})
  gHv = Array{Float64}(nlp.meta.ncon);
  return ghjvprod!(nlp, x, g, v, gHv)
end

"""Compute the vector of dot products `(g, Hj*v)` in place
where `Hj` is the Hessian of the `j`-th constraint at `x`.
"""
function ghjvprod!(nlp :: AmplModel,
                   x :: Vector{Float64}, g :: Vector{Float64},
                   v :: Vector{Float64}, gHv :: Vector{Float64})
  # Note: x is in fact not used.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(g) >= nlp.meta.nvar || error("g must have length at least $(nlp.meta.nvar)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")

  if nlp.safe
    _ = cons(nlp, x) ; nlp.counters.neval_cons -= 1
  end
  @asl_call(:asl_ghjvprod, Ptr{Float64},
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
             nlp.__asl, g,            v,            gHv);
  nlp.counters.neval_hprod += nlp.meta.ncon
end

"""Evaluate the Lagrangian Hessian at `(x,y)` in sparse coordinate format.
Only the lower triangle is returned.
"""
function hess_coord(nlp :: AmplModel,
                    x :: Vector{Float64};
                    y :: Vector{Float64} = nlp.meta.y0,
                    obj_weight :: Float64 = 1.0)
  # Note: x is in fact not used.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(y) >= nlp.meta.ncon || error("y must have length at least $(nlp.meta.ncon)")

  if nlp.safe
    _ = obj(nlp, x) ; nlp.counters.neval_obj -= 1
    _ = cons(nlp, x) ; nlp.counters.neval_cons -= 1
  end

  rows = Array{Int64}(nlp.meta.nnzh)
  cols = Array{Int64}(nlp.meta.nnzh)
  vals = Array{Float64}(nlp.meta.nnzh)
  @asl_call(:asl_hess, Void,
            (Ptr{Void}, Ptr{Float64}, Float64,    Ptr{Int64}, Ptr{Int64}, Ptr{Float64}),
             nlp.__asl, y,            obj_weight, rows,       cols,       vals)
  nlp.counters.neval_hess += 1
  # Use 1-based indexing.
  # Swap rows and cols to obtain the lower triangle.
  return (cols+1, rows+1, vals)
end

"""Evaluate the Lagrangian Hessian at `(x,y)` as a sparse matrix.
Only the lower triangle is returned.
"""
function hess(nlp :: AmplModel,
              x :: Vector{Float64};
              y :: Vector{Float64} = nlp.meta.y0,
              obj_weight :: Float64 = 1.0)
  @check_ampl_model
  (rows, cols, vals) = hess_coord(nlp, x, y=y, obj_weight=obj_weight);
  return sparse(rows, cols, vals, nlp.meta.nvar, nlp.meta.nvar)
end
