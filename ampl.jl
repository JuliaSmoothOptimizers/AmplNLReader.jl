# A rudimentary AMPL interface in Julia.
# D. Orban, Vancouver, April 2014.

include("ampl_utils.jl")

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

  __asl :: Ptr{Void}         # pointer to internal ASL structure. Do not touch.

  nvar  :: Int               # number of variables
  ncon  :: Int               # total number of general constraints
  x0    :: Array{Float64,1}  # initial guess
  y0    :: Array{Float64,1}  # initial Lagrange multipliers

  lvar  :: Array{Float64,1}  # vector of lower bounds
  uvar  :: Array{Float64,1}  # vector of upper bounds

  ifix  :: Array{Int64,1}    # indices of fixed variables
  ilow  :: Array{Int64,1}    # indices of variables with lower bound only
  iupp  :: Array{Int64,1}    # indices of variables with upper bound only
  irng  :: Array{Int64,1}    # indices of variables with lower and upper bound (range)
  ifree :: Array{Int64,1}    # indices of free variables
  iinf  :: Array{Int64,1}    # indices of infeasible bounds

  lcon  :: Array{Float64,1}  # vector of constraint lower bounds
  ucon  :: Array{Float64,1}  # vector of constraint upper bounds

  jfix  :: Array{Int64,1}    # indices of equality constraints
  jlow  :: Array{Int64,1}    # indices of constraints of the form c(x) ≥ cl
  jupp  :: Array{Int64,1}    # indices of constraints of the form c(x) ≤ cu
  jrng  :: Array{Int64,1}    # indices of constraints of the form cl ≤ c(x) ≤ cu
  jfree :: Array{Int64,1}    # indices of "free" constraints (there shouldn't be any)
  jinf  :: Array{Int64,1}    # indices of the visibly infeasible constraints

  nlin  :: Int               # number of linear constraints
  nnln  :: Int               # number of nonlinear general constraints
  nnet  :: Int               # number of nonlinear network constraints

  lin   :: Range1{Int64}     # indices of linear constraints
  nln   :: Range1{Int64}     # indices of nonlinear constraints
  net   :: Range1{Int64}     # indices of nonlinear network constraints

  nnzj  :: Int               # number of nonzeros in the sparse Jacobian
  nnzh  :: Int               # number of nonzeros in the sparse Hessian

  minimize :: Bool           # true if optimize == minimize
  islp  :: Bool              # true if the problem is a linear program
  name  :: ASCIIString       # problem name

  function AmplModel(stub :: ASCIIString)
    asl = @jampl_call(:jampl_init, Ptr{Void}, (Ptr{Uint8},), stub);
    if asl == C_NULL
      msg = "Error allocating ASL structure"
      error(msg)
    end

    minimize = !bool(@jampl_call(:jampl_objtype, Int32, (Ptr{Void},), asl));
    islp = bool(@jampl_call(:jampl_islp, Int32, (Ptr{Void},), asl));

    nvar = int(@jampl_call(:jampl_nvar, Int32, (Ptr{Void},), asl));
    ncon = int(@jampl_call(:jampl_ncon, Int32, (Ptr{Void},), asl));

    x0   = pointer_to_array(@jampl_call(:jampl_x0,   Ptr{Float64}, (Ptr{Void},), asl),
                            (nvar,), false);
    y0   = pointer_to_array(@jampl_call(:jampl_y0,   Ptr{Float64}, (Ptr{Void},), asl),
                            (ncon,), false);

    lvar = pointer_to_array(@jampl_call(:jampl_lvar, Ptr{Float64}, (Ptr{Void},), asl),
                            (nvar,), false);
    uvar = pointer_to_array(@jampl_call(:jampl_uvar, Ptr{Float64}, (Ptr{Void},), asl),
                            (nvar,), false);

    ifix = find(lvar .== uvar);
    ilow = find((lvar .> -Inf) & (uvar .== Inf));
    iupp = find((lvar .== -Inf) & (uvar .< Inf));
    irng = find((lvar .> -Inf) & (uvar .< Inf) & (lvar .< uvar));
    ifree = find((lvar .== -Inf) & (uvar .== Inf));
    iinf = find(lvar .> uvar);

    lcon = pointer_to_array(@jampl_call(:jampl_lcon, Ptr{Float64}, (Ptr{Void},), asl),
                            (ncon,), false);
    ucon = pointer_to_array(@jampl_call(:jampl_ucon, Ptr{Float64}, (Ptr{Void},), asl),
                            (ncon,), false);

    jfix = find(lcon .== ucon);
    jlow = find((lcon .> -Inf) & (ucon .== Inf));
    jupp = find((lcon .== -Inf) & (ucon .< Inf));
    jrng = find((lcon .> -Inf) & (ucon .< Inf) & (lcon .< ucon));
    jfree = find((lcon .== -Inf) & (ucon .== Inf));
    jinf = find(lcon .> ucon);

    nnet = int(@jampl_call(:jampl_nlnc, Int32, (Ptr{Void},), asl));
    nnln = int(@jampl_call(:jampl_nlc,  Int32, (Ptr{Void},), asl)) - nnet;
    nlin = ncon - nnln - nnet

    nln  = 1 : nnln
    net  = nnln+1 : nnln+nnet
    lin  = nnln+nnet+1 : ncon

    nnzj = int(@jampl_call(:jampl_nnzj, Int32, (Ptr{Void},), asl));
    nnzh = int(@jampl_call(:jampl_nnzh, Int32, (Ptr{Void},), asl));

    nlp = new(asl, nvar, ncon, x0, y0,
              lvar, uvar, ifix, ilow, iupp, irng, ifree, iinf,
              lcon, ucon, jfix, jlow, jupp, jrng, jfree, jinf,
              nlin, nnln, nnet, lin, nln, net,
              nnzj, nnzh, minimize, islp, stub)

    lagscale(nlp, -1.0)  # Lagrangian L(x,y) = f(x) - ∑ yi ci(x)

    finalizer(nlp, amplmodel_finalize)
    return nlp
  end

end

# Methods associated to AmplModel instances.

function write_sol(nlp :: AmplModel, msg :: ASCIIString, x :: Array{Float64,1}, y :: Array{Float64,1})
  if length(x) != nlp.nvar
    error("x must have length $(nlp.nvar)")
  end
  if length(y) != nlp.ncon
    error("y must have length $(nlp.ncon)")
  end
  @jampl_call(:jampl_write_sol, Void,
                                (Ptr{Void}, Ptr{Uint8}, Ptr{Float64}, Ptr{Float64}),
                                nlp.__asl,  msg,        x,            y)
end

function amplmodel_finalize(nlp :: AmplModel)
  @jampl_call(:jampl_finalize, Void, (Ptr{Void},), nlp.__asl)
end

# Displaying AmplModel instances.

import Base.show, Base.print
function show(io :: IO, nlp :: AmplModel)
  s  = sprint_formatted (nlp.minimize ? "Minimization " : "Maximization ")
  s *= @sprintf("problem %s\n", nlp.name)
  s *= @sprintf("nvar = %d, ncon = %d (%d linear)\n", nlp.nvar, nlp.ncon, nlp.nlin)
  print(io, s)
end

function print(io :: IO, nlp :: AmplModel)
  print_formatted (nlp.minimize ? "Minimization " : "Maximization ")
  @printf("problem %s\n", nlp.name)
  @printf("nvar = %d, ncon = %d (%d linear)\n", nlp.nvar, nlp.ncon, nlp.nlin)
  @printf("lvar = "); display(nlp.lvar'); @printf("\n")
  @printf("uvar = "); display(nlp.uvar'); @printf("\n")
  @printf("lcon = "); display(nlp.lcon'); @printf("\n")
  @printf("ucon = "); display(nlp.ucon'); @printf("\n")
  @printf("x0 = ");   display(nlp.x0'); @printf("\n")
  @printf("y0 = ");   display(nlp.y0'); @printf("\n")
  if nlp.nlin > 0
    @printf("linear constraints:    "); display(nlp.lin); @printf("\n");
  end
  if nlp.nnln > 0
    @printf("nonlinear constraints: "); display(nlp.nln); @printf("\n");
  end
  if nlp.nnet > 0
    @printf("network constraints:   "); display(nlp.net); @printf("\n");
  end
end

# Scaling AmplModel instances.

function varscale(nlp :: AmplModel, s :: Array{Float64,1})
  # Scale the vector of variables by the vector s.
  if length(s) < nlp.nvar
    error("s must have length at least $(nlp.nvar)")
  end
  @jampl_call(:jampl_varscale, Void,
              (Ptr{Void}, Ptr{Float64}),
               nlp.__asl, s)
end

function lagscale(nlp :: AmplModel, s :: Float64)
  # Set the scaling factor σ in the Lagrangian:
  # L(x,y) = f(x) + σ ∑ yi ci(x).
  @jampl_call(:jampl_lagscale, Void, (Ptr{Void}, Float64), nlp.__asl, s)
end

function conscale(nlp :: AmplModel, s :: Array{Float64,1})
  # Scale the vector of constraints by the vector s.
  if length(s) < nlp.ncon
    error("s must have length at least $(nlp.ncon)")
  end
  @jampl_call(:jampl_conscale, Void,
              (Ptr{Void}, Ptr{Float64}), nlp.__asl, s)
end

# Evaluating objective, constraints and derivatives.

function obj(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the objective function at x.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  @jampl_call(:jampl_obj, Float64, (Ptr{Void}, Ptr{Float64}), nlp.__asl, x)
end

function grad(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the objective function gradient at x.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  g = Array(Float64, nlp.nvar)
  @jampl_call(:jampl_grad, Ptr{Float64}, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), nlp.__asl, x, g)
  return g
end

function cons(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the vector of constraints at x.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  c = Array(Float64, nlp.ncon)
  @jampl_call(:jampl_cons, Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), nlp.__asl, x, c)
  return c
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
  @jampl_call(:jampl_jcon, Float64, (Ptr{Void}, Ptr{Float64}, Int32), nlp.__asl, x, j-1)
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
  pointer_to_array(@jampl_call(:jampl_jcongrad, Ptr{Float64},
                               (Ptr{Void}, Ptr{Float64}, Int32),
                                nlp.__asl, x,            j-1),
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
  nnz = @jampl_call(:jampl_sparse_congrad_nnz, Csize_t,
                     (Ptr{Void}, Cint), nlp.__asl, j-1)
  inds = Array(Int64, nnz)
  vals = Array(Float64, nnz)
  @jampl_call(:jampl_sparse_congrad, Void,
                           (Ptr{Void}, Ptr{Float64}, Int32, Ptr{Int64}, Ptr{Float64}),
                            nlp.__asl, x,            j-1,   inds,       vals)
  sparsevec(inds, vals, nlp.nvar)
end

function jac(nlp :: AmplModel, x :: Array{Float64,1})
  # Evaluate the sparse Jacobian of the constraints at x.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  
  rows = Array(Int64, nlp.nnzj)
  cols = Array(Int64, nlp.nnzj)
  vals = Array(Float64, nlp.nnzj)
  @jampl_call(:jampl_jac, Void, (Ptr{Void}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}), nlp.__asl, x, rows, cols, vals)
  sparse(rows, cols, vals, nlp.ncon, nlp.nvar)
end

function hprod(nlp :: AmplModel,
               x :: Array{Float64,1},
               v :: Array{Float64,1};
               y :: Array{Float64,1} = nlp.y0,
               obj_weight :: Float64 = 1.0)
  # Evaluate the product of the Hessian of the Lagrangian at (x,y) with v.
  # Note: x is in fact not used.
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
                               (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Float64),
                                nlp.__asl, y,            v,            obj_weight),
                   (nlp.nvar,), true)
end

function jth_hprod(nlp :: AmplModel,
                x :: Array{Float64,1}, v :: Array{Float64,1}, j :: Int)
  # Compute the product of the Hessian of the j-th constraint at x with v.
  # If j=0, compute the product of the Hessian of the objective at x with v.
  # Note: x is in fact not used.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  if length(v) < nlp.nvar
    error("v must have length at least $(nlp.nvar)")
  end
  if j < 0 | j > nlp.ncon
    error("expected 0 ≤ j ≤ $(nlp.ncon)")
  end
  # Require that Julia be responsible for freeing up this chunk of memory.
  Hv = pointer_to_array(@jampl_call(:jampl_hvcompd, Ptr{Float64},
                                     (Ptr{Void}, Ptr{Float64}, Int),
                                      nlp.__asl, v,            j-1),
                         (nlp.nvar,), true)
  return (j > 0) ? -Hv : Hv  # lagscale() flipped the sign of each constraint.
end

function ghjvprod(nlp :: AmplModel,
                  x :: Array{Float64,1}, g :: Array{Float64,1}, v :: Array{Float64,1})
  # Compute the vector of dot products (g, Hj*v)
  # where Hj is the Hessian of the j-th constraint at x.
  # Note: x is in fact not used.
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
  gHv = pointer_to_array(@jampl_call(:jampl_ghjvprod, Ptr{Float64},
                                     (Ptr{Void}, Ptr{Float64}, Ptr{Float64}),
                                      nlp.__asl, g,            v),
                         (nlp.ncon,), true)
  return -gHv  # lagscale() flipped the sign of each constraint.
end

function hess(nlp :: AmplModel,
              x :: Array{Float64,1};
              y :: Array{Float64,1} = nlp.y0,
              obj_weight :: Float64 = 1.0)
  # Evaluate the sparse Hessian of the Lagrangian at (x,y).
  # Note: x is in fact not used.
  if length(x) < nlp.nvar
    error("x must have length at least $(nlp.nvar)")
  end
  if length(y) < nlp.ncon
    error("y must have length at least $(nlp.ncon)")
  end
  rows = Array(Int64, nlp.nnzh)
  cols = Array(Int64, nlp.nnzh)
  vals = Array(Float64, nlp.nnzh)
  @jampl_call(:jampl_hess, Void,
                                 (Ptr{Void}, Ptr{Float64}, Float64, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}),
                                  nlp.__asl, y, obj_weight, rows, cols, vals)
  sparse(rows, cols, vals, nlp.nvar, nlp.nvar)
end
