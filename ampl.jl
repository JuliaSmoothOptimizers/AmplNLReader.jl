# A rudimentary AMPL interface in Julia.
# D. Orban, Vancouver, April 2014.

# Convenience macros and libasl specifics.
jampl = "jampl";
macro jampl_call(func, args...)
  quote
    ccall(($func, $jampl), $(args...))
  end
end

# TODO: When / If Julia allows abstract types with fields,
# AmplModel should inherit from NLPModel.
# include("nlp.jl")  # import NLPModel.

type AmplModel

  # A composite type that represents the optimization problem
  #
  #  optimize   obj(x)
  #  subject to lvar ≤    x    ≤ uvar
  #             lcon ≤ cons(x) ≤ ucon
  #
  # where x        is an nvar-dimensional vector,
  #       obj      is the real-valued objective function,
  #       cons     is the vector-valued constraint function,
  #       optimize is either "minimize" or "maximize".
  #
  # Here, lvar, uvar, lcon and ucon are vectors. Some of their
  # components may be infinite to indicate that the corresponding
  # bound or general constraint is not present.

  nvar :: Int               # number of variables
  x0   :: Array{Float64,1}  # initial guess
  lvar :: Array{Float64,1}  # vector of lower bounds
  uvar :: Array{Float64,1}  # vector of upper bounds

  ncon :: Int               # total number of general constraints
  nlc  :: Int               # number of linear constraints
  nlnc :: Int               # number of nonlinear network constraints
  y0   :: Array{Float64,1}  # initial Lagrange multipliers
  lcon :: Array{Float64,1}  # vector of constraint lower bounds
  ucon :: Array{Float64,1}  # vector of constraint upper bounds

  lin  :: Range1{Int64}     # indices of linear constraints
  nlin :: Range1{Int64}     # indices of nonlinear constraints (not network)

  nnzj :: Int               # number of nonzeros in the sparse Jacobian
  nnzh :: Int               # number of nonzeros in the sparse Hessian

  minimize :: Bool          # true if optimize == minimize
  islp :: Bool              # true if the problem is a linear program
  name :: ASCIIString       # problem name

  function AmplModel(stub :: ASCIIString)
    init_status = @jampl_call(:jampl_init, Int32, (Ptr{Uint8},), stub);
    if init_status != 0
      msg = "Error allocating ASL structure"
      error(msg)
    end

    minimize = !bool(@jampl_call(:jampl_objtype, Int32, ()));
    islp = bool(@jampl_call(:jampl_islp, Int32, ()));

    nvar = int(@jampl_call(:jampl_nvar, Int32, ()));
    x0   = pointer_to_array(@jampl_call(:jampl_x0,   Ptr{Float64}, ()), (nvar,), false);
    lvar = pointer_to_array(@jampl_call(:jampl_lvar, Ptr{Float64}, ()), (nvar,), false);
    uvar = pointer_to_array(@jampl_call(:jampl_uvar, Ptr{Float64}, ()), (nvar,), false);

    ncon = int(@jampl_call(:jampl_ncon, Int32, ()));
    nlc  = int(@jampl_call(:jampl_nlc,  Int32, ()));
    nlnc = int(@jampl_call(:jampl_nlnc, Int32, ()));
    y0   = pointer_to_array(@jampl_call(:jampl_y0,   Ptr{Float64}, ()), (ncon,), false);
    lcon = pointer_to_array(@jampl_call(:jampl_lcon, Ptr{Float64}, ()), (ncon,), false);
    ucon = pointer_to_array(@jampl_call(:jampl_ucon, Ptr{Float64}, ()), (ncon,), false);

    lin  = nlc+nlnc+1 : ncon
    nlin = 1 : nlc

    nnzj = int(@jampl_call(:jampl_nnzj, Int32, ()));
    nnzh = int(@jampl_call(:jampl_nnzh, Int32, ()));

    nlp = new(nvar, x0, lvar, uvar, ncon, nlc, nlnc, y0, lcon, ucon,
              lin, nlin, nnzj, nnzh, minimize, islp, stub)

    lagscale(nlp, -1.0)  # Lagrangian L(x,y) = f(x) - ∑ yi ci(x)

    finalizer(nlp, amplmodel_finalize)
    return nlp
  end

end

# Methods associated to AmplModel objects.

function amplmodel_finalize(nlp :: AmplModel)
  @jampl_call(:jampl_finalize, Void, ())
end

function varscale(nlp :: AmplModel, s :: Array{Float64,1})
  # Scale the vector of variables by the vector s.
  if length(s) < nlp.nvar
    error("s must have length at least $(nlp.nvar)")
  end
  @jampl_call(:jampl_varscale, Void, (Ptr{Float64},), s)
end

function obj(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the objective function at x.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  @jampl_call(:jampl_obj, Float64, (Ptr{Float64},), x)
end

function grad(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the objective function gradient at x.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  # Require that Julia be responsible for freeing up this chunk of memory.
  pointer_to_array(@jampl_call(:jampl_grad, Ptr{Float64}, (Ptr{Float64},), x),
                   (nlp.nvar,), true)
end

function lagscale(nlp :: AmplModel, s :: Float64)
  # Set the scaling factor σ in the Lagrangian:
  # L(x,y) = f(x) + σ ∑ yi ci(x).
  @jampl_call(:jampl_lagscale, Void, (Float64,), s)
end

function conscale(nlp :: AmplModel, s :: Array{Float64,1})
  # Scale the vector of constraints by the vector s.
  if length(s) < nlp.ncon
    error("s must have length at least $(nlp.ncon)")
  end
  @jampl_call(:jampl_conscale, Void, (Ptr{Float64},), s)
end

function cons(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the vector of constraints at x.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  # Require that Julia be responsible for freeing up this chunk of memory.
  pointer_to_array(@jampl_call(:jampl_cons, Ptr{Float64}, (Ptr{Float64},), x),
                   (nlp.ncon,), true)
end

function jth_con(nlp :: AmplModel, x :: Array{Float64,1}, j :: Int)
  # Evaluate the j-th constraint at x.
  if (j < 1) || (j > nlp.ncon)
    msg = "Invalid constraint index $j"
    error(msg)
  end
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  @jampl_call(:jampl_jcon, Float64, (Ptr{Float64}, Int32), x, j-1)
end

function jth_congrad(nlp :: AmplModel, x :: Array{Float64,1}, j :: Int)
  # Evaluate the j-th constraint gradient at x.
  if (j < 1) || (j > nlp.ncon)
    msg = "Invalid constraint index $j"
    error(msg)
  end
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  # Require that Julia be responsible for freeing up this chunk of memory.
  pointer_to_array(@jampl_call(:jampl_jcongrad, Ptr{Float64}, (Ptr{Float64}, Int32), x, j-1),
                   (nlp.nvar,), true)
end

function jth_sparse_congrad(nlp :: AmplModel, x :: Array{Float64,1}, j :: Int)
  # Evaluate the j-th constraint sparse gradient at x.
  if (j < 1) || (j > nlp.ncon)
    msg = "Invalid constraint index $j"
    error(msg)
  end
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  inds, vals = @jampl_call(:jampl_sparse_congrad, Any, (Ptr{Float64},Int), x, j-1)
  sparsevec(inds, vals, nlp.nvar)
end

function jac(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the sparse Jacobian of the constraints at x.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  rows, cols, vals = @jampl_call(:jampl_jac, Any, (Ptr{Float64},), x)
  sparse(rows, cols, vals, nlp.ncon, nlp.nvar)
end

function hprod(nlp :: AmplModel,
               x :: Array{Float64,1},
               v :: Array{Float64,1};
               y :: Array{Float64,1} = nlp.y0,
               obj_weight :: Float64 = 1.0)
  # Evaluate the product of the Hessian of the Lagrangian at (x,y) with v.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  if length(y) < nlp.ncon
    error("y must have length at least $(nlp.ncon)")
  end
  if length(v) < nlp.nvar
    error("v must have length at least $(nlp.nvar)")
  end
  # Require that Julia be responsible for freeing up this chunk of memory.
  pointer_to_array(@jampl_call(:jampl_hprod, Ptr{Float64},
                               (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
                               x, y, v, obj_weight),
                   (nlp.nvar,), true)
end

function jth_hprod(nlp :: AmplModel,
                x :: Array{Float64,1}, v :: Array{Float64,1}, j :: Int)
  # Compute the product of the Hessian of the j-th constraint at x with v.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  if length(v) < nlp.nvar
    error("v must have length at least $(nlp.nvar)")
  end
  ej = zeros(nlp.ncon,);
  ej[j] = 1;
  hprod(nlp, x, v, y=ej);
end

function ghjvprod(nlp :: AmplModel,
                  x :: Array{Float64,1}, g :: Array{Float64,1}, v :: Array{Float64,1})
  # Compute the vector of dot products (g, Hj*v)
  # where Hj is the Hessian of the j-th constraint at x.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  if length(g) < nlp.nvar
    error("g must have length at least $(nlp.nvar)")
  end
  if length(v) < nlp.nvar
    error("v must have length at least $(nlp.nvar)")
  end
  # Require that Julia be responsible for freeing up this chunk of memory.
  pointer_to_array(@jampl_call(:jampl_ghjvprod, Ptr{Float64},
                               (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
                               x, g, v),
                   (nlp.ncon,), true)
end

function hess(nlp :: AmplModel,
              x :: Array{Float64,1};
              y :: Array{Float64,1} = nlp.y0,
              obj_weight :: Float64 = 1.0)
  # Evaluate the sparse Hessian of the Lagrangian at (x,y).
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  if length(y) < nlp.ncon
    error("y must have length at least $(nlp.ncon)")
  end
  rows, cols, vals = @jampl_call(:jampl_hess, Any,
                                 (Ptr{Float64}, Ptr{Float64}, Float64),
                                  x,            y,            obj_weight)
  sparse(rows, cols, vals, nlp.nvar, nlp.nvar)
end
