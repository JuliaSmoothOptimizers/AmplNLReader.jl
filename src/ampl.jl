# A rudimentary AMPL interface in Julia.
# D. Orban, Vancouver, April 2014.
module ampl

using NLP  # Defines NLPModelMeta.

export AmplModel, AmplException,
       write_sol, amplmodel_finalize, varscale, lagscale, conscale,
       obj, grad, cons, jth_con, jth_congrad, jth_sparse_congrad,
       jac_coord, jac, hprod, jth_hprod, ghjvprod, hess_coord, hess

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
  include("../deps/deps.jl")
else
  error("ASL library not properly installed. Please run Pkg.build(\"ampl\")")
end

# Convenience macro.
const asl = "libasl";
macro asl_call(func, args...)
  quote
    ccall(($func, $asl), $(args...))
  end
end

type AmplException
  msg :: ASCIIString
end

macro check_ampl_model()
  :(nlp.__asl == C_NULL && throw(AmplException("Uninitialized AMPL model")))
end

type AmplModel
  meta  :: NLPModelMeta;     # Problem metadata.
  __asl :: Ptr{Void};        # Pointer to internal ASL structure. Do not touch.

  function AmplModel(stub :: ASCIIString)
    asl = @asl_call(:asl_init, Ptr{Void}, (Ptr{Uint8},), stub);
    asl == C_NULL && error("Error allocating ASL structure")

    minimize = !bool(@asl_call(:asl_objtype, Int32, (Ptr{Void},), asl));
    islp = bool(@asl_call(:asl_islp, Int32, (Ptr{Void},), asl));

    nvar = int(@asl_call(:asl_nvar, Int32, (Ptr{Void},), asl));
    ncon = int(@asl_call(:asl_ncon, Int32, (Ptr{Void},), asl));

    x0   = pointer_to_array(@asl_call(:asl_x0,   Ptr{Float64}, (Ptr{Void},), asl),
                            (nvar,), false);
    y0   = pointer_to_array(@asl_call(:asl_y0,   Ptr{Float64}, (Ptr{Void},), asl),
                            (ncon,), false);

    lvar = pointer_to_array(@asl_call(:asl_lvar, Ptr{Float64}, (Ptr{Void},), asl),
                            (nvar,), false);
    uvar = pointer_to_array(@asl_call(:asl_uvar, Ptr{Float64}, (Ptr{Void},), asl),
                            (nvar,), false);

    lcon = pointer_to_array(@asl_call(:asl_lcon, Ptr{Float64}, (Ptr{Void},), asl),
                            (ncon,), false);
    ucon = pointer_to_array(@asl_call(:asl_ucon, Ptr{Float64}, (Ptr{Void},), asl),
                            (ncon,), false);

    nnet = int(@asl_call(:asl_nlnc, Int32, (Ptr{Void},), asl));
    nnln = int(@asl_call(:asl_nlc,  Int32, (Ptr{Void},), asl)) - nnet;
    nlin = ncon - nnln - nnet

    nln  = 1 : nnln
    net  = nnln+1 : nnln+nnet
    lin  = nnln+nnet+1 : ncon

    nnzj = int(@asl_call(:asl_nnzj, Int32, (Ptr{Void},), asl));
    nnzh = int(@asl_call(:asl_nnzh, Int32, (Ptr{Void},), asl));

    meta = NLPModelMeta(nvar, x0=x0, lvar=lvar, uvar=uvar,
                        ncon=ncon, y0=y0, lcon=lcon, ucon=ucon,
                        nnzj=nnzj, nnzh=nnzh,
                        lin=lin, nln=nln, net=net,
                        nlin=nlin, nnln=nnln, nnet=nnet,
                        minimize=minimize, islp=islp, name=stub);

    nlp = new(meta, asl);

    lagscale(nlp, -1.0)  # Lagrangian L(x,y) = f(x) - ∑ yi ci(x)

    finalizer(nlp, amplmodel_finalize)
    return nlp
  end

end

# Methods associated to AmplModel instances.

function write_sol(nlp :: AmplModel, msg :: ASCIIString, x :: Array{Float64,1}, y :: Array{Float64,1})
  @check_ampl_model
  length(x) == nlp.meta.nvar || error("x must have length $(nlp.meta.nvar)")
  length(y) == nlp.meta.ncon || error("y must have length $(nlp.meta.ncon)")

  @asl_call(:asl_write_sol, Void,
            (Ptr{Void}, Ptr{Uint8}, Ptr{Float64}, Ptr{Float64}),
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

import Base.show, Base.print
function show(io :: IO, nlp :: AmplModel)
  @check_ampl_model
  show(io, nlp.meta);
end

function print(io :: IO, nlp :: AmplModel)
  @check_ampl_model
  print(io, nlp.meta);
end

# Scaling AmplModel instances.

function varscale(nlp :: AmplModel, s :: Array{Float64,1})
  # Scale the vector of variables by the vector s.
  @check_ampl_model
  length(s) >= nlp.meta.nvar || error("s must have length at least $(nlp.meta.nvar)")

  @asl_call(:asl_varscale, Void, (Ptr{Void}, Ptr{Float64}), nlp.__asl, s)
end

function lagscale(nlp :: AmplModel, s :: Float64)
  # Set the scaling factor σ in the Lagrangian:
  # L(x,y) = f(x) + σ ∑ yi ci(x).
  @check_ampl_model
  @asl_call(:asl_lagscale, Void, (Ptr{Void}, Float64), nlp.__asl, s)
end

function conscale(nlp :: AmplModel, s :: Array{Float64,1})
  # Scale the vector of constraints by the vector s.
  @check_ampl_model
  length(s) >= nlp.meta.ncon || error("s must have length at least $(nlp.meta.ncon)")

  @asl_call(:asl_conscale, Void, (Ptr{Void}, Ptr{Float64}), nlp.__asl, s)
end

# Evaluating objective, constraints and derivatives.

function obj(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the objective function at x.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  @asl_call(:asl_obj, Float64, (Ptr{Void}, Ptr{Float64}), nlp.__asl, x)
end

function grad(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the objective function gradient at x.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  g = Array(Float64, nlp.meta.nvar)
  @asl_call(:asl_grad, Ptr{Float64}, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), nlp.__asl, x, g)
  return g
end

function cons(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the vector of constraints at x.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  c = Array(Float64, nlp.meta.ncon)
  @asl_call(:asl_cons, Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), nlp.__asl, x, c)
  return c
end

function jth_con(nlp :: AmplModel, x :: Array{Float64,1}, j :: Int)
  # Evaluate the j-th constraint at x.
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  @asl_call(:asl_jcon, Float64, (Ptr{Void}, Ptr{Float64}, Int32), nlp.__asl, x, j-1)
end

function jth_congrad(nlp :: AmplModel, x :: Array{Float64,1}, j :: Int)
  # Evaluate the j-th constraint gradient at x.
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  g = Array(Float64, nlp.meta.nvar)
  @asl_call(:asl_jcongrad, Ptr{Float64},
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Int32),
             nlp.__asl, x,            g,            j-1)
  return g
end

function jth_sparse_congrad(nlp :: AmplModel, x :: Array{Float64,1}, j :: Int)
  # Evaluate the j-th constraint sparse gradient at x.
  @check_ampl_model
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  nnz = @asl_call(:asl_sparse_congrad_nnz, Csize_t,
                  (Ptr{Void}, Cint), nlp.__asl, j-1)
  inds = Array(Int64, nnz)
  vals = Array(Float64, nnz)
  @asl_call(:asl_sparse_congrad, Void,
            (Ptr{Void}, Ptr{Float64}, Int32, Ptr{Int64}, Ptr{Float64}),
             nlp.__asl, x,            j-1,   inds,       vals)
  # Use 1-based indexing.
  return sparsevec(inds+1, vals, nlp.meta.nvar)
end

function jac_coord(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the sparse Jacobian of the constraints at x in coordinate format.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")

  rows = Array(Int64, nlp.meta.nnzj)
  cols = Array(Int64, nlp.meta.nnzj)
  vals = Array(Float64, nlp.meta.nnzj)
  @asl_call(:asl_jac, Void,
            (Ptr{Void}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}),
             nlp.__asl, x,            rows,       cols,       vals)
  # Use 1-based indexing.
  return (rows+1, cols+1, vals)
end

function jac(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the sparse Jacobian of the constraints.
  @check_ampl_model
  (rows, cols, vals) = jac_coord(nlp, x);
  return sparse(rows, cols, vals, nlp.meta.ncon, nlp.meta.nvar)
end

function hprod(nlp :: AmplModel,
               x :: Array{Float64,1},
               v :: Array{Float64,1};
               y :: Array{Float64,1} = nlp.meta.y0,
               obj_weight :: Float64 = 1.0)
  # Evaluate the product of the Hessian of the Lagrangian at (x,y) with v.
  # Note: x is in fact not used.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(y) >= nlp.meta.ncon || error("y must have length at least $(nlp.meta.ncon)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")

  hv = Array(Float64, nlp.meta.nvar);
  @asl_call(:asl_hprod, Ptr{Float64},
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
             nlp.__asl, y,            v,            hv,           obj_weight);
  return hv
end

function jth_hprod(nlp :: AmplModel,
                   x :: Array{Float64,1}, v :: Array{Float64,1}, j :: Int)
  # Compute the product of the Hessian of the j-th constraint at x with v.
  # If j=0, compute the product of the Hessian of the objective at x with v.
  # Note: x is in fact not used.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")
  (1 <= j <= nlp.meta.ncon)  || error("expected 0 ≤ j ≤ $(nlp.meta.ncon)")

  hv = Array(Float64, nlp.meta.nvar);
  @asl_call(:asl_hvcompd, Ptr{Float64},
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Int),
             nlp.__asl, v,            hv,           j-1);
  return (j > 0) ? -hv : hv  # lagscale() flipped the sign of each constraint.
end

function ghjvprod(nlp :: AmplModel,
                  x :: Array{Float64,1}, g :: Array{Float64,1}, v :: Array{Float64,1})
  # Compute the vector of dot products (g, Hj*v)
  # where Hj is the Hessian of the j-th constraint at x.
  # Note: x is in fact not used.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(g) >= nlp.meta.nvar || error("g must have length at least $(nlp.meta.nvar)")
  length(v) >= nlp.meta.nvar || error("v must have length at least $(nlp.meta.nvar)")

  gHv = Array(Float64, nlp.meta.ncon);
  @asl_call(:asl_ghjvprod, Ptr{Float64},
            (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
             nlp.__asl, g,            v,            gHv);
  return -gHv  # lagscale() flipped the sign of each constraint.
end

function hess_coord(nlp :: AmplModel,
                    x :: Array{Float64,1};
                    y :: Array{Float64,1} = nlp.y0,
                    obj_weight :: Float64 = 1.0)
  # Evaluate the sparse Hessian of the Lagrangian at (x,y) in coordinate format.
  # Note: x is in fact not used.
  @check_ampl_model
  length(x) >= nlp.meta.nvar || error("x must have length at least $(nlp.meta.nvar)")
  length(y) >= nlp.meta.ncon || error("y must have length at least $(nlp.meta.ncon)")

  rows = Array(Int64, nlp.meta.nnzh)
  cols = Array(Int64, nlp.meta.nnzh)
  vals = Array(Float64, nlp.meta.nnzh)
  @asl_call(:asl_hess, Void,
            (Ptr{Void}, Ptr{Float64}, Float64,    Ptr{Int64}, Ptr{Int64}, Ptr{Float64}),
             nlp.__asl, y,            obj_weight, rows,       cols,       vals)
  # Use 1-based indexing.
  return (rows+1, cols+1, vals)
end

function hess(nlp :: AmplModel,
              x :: Array{Float64,1};
              y :: Array{Float64,1} = nlp.y0,
              obj_weight :: Float64 = 1.0)
  @check_ampl_model
  (rows, cols, vals) = hess_coord(nlp, x, y=y, obj_weight=obj_weight);
  return sparse(rows, cols, vals, nlp.meta.nvar, nlp.meta.nvar)
end

end  # Module ampl
