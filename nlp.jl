type NLPModel  # immutable?

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

  ncon :: Int               # number of general constraints
  y0   :: Array{Float64,1}  # initial Lagrange multipliers
  lcon :: Array{Float64,1}  # vector of constraint lower bounds
  ucon :: Array{Float64,1}  # vector of constraint upper bounds

  nnzj :: Int               # number of nonzeros in the sparse Jacobian
  nnzh :: Int               # number of nonzeros in the sparse Hessian

  minimize :: Bool          # true if optimize == minimize
  islp :: Bool              # true if the problem is a linear program
  name :: ASCIIString       # nroblem name

  function NLPModel(nvar;
                    x0=zeros(nvar,),
                    lvar=-Inf * ones(nvar,),
                    uvar=Inf * ones(nvar,),
                    ncon=0,
                    y0=zeros(ncon,),
                    lcon=-Inf * ones(ncon,),
                    ucon=Inf * ones(ncon,),
                    nnzj=0,
                    nnzh=0,
                    minimize=true,
                    islp=false,
                    name="Generic")
    if (nvar < 1) || (ncon < 0)
      error("Nonsensical dimensions")
    end
    if size(x0) != (nvar,)
      error("x0 has size inconsistent with nvar")
    end
    if size(lvar) != (nvar,)
      error("lvar has size inconsistent with nvar")
    end
    if size(uvar) != (nvar,)
      error("uvar has size inconsistent with nvar")
    end
    if size(y0) != (ncon,)
      error("y0 has size inconsistent with ncon")
    end
    if size(lcon) != (ncon,)
      error("lcon has size inconsistent with ncon")
    end
    if size(ucon) != (ncon,)
      error("ucon has size inconsistent with ncon")
    end

    new(nvar, x0, lvar, uvar, ncon, y0, lcon, ucon,
        nnzj, nnzh, minimize, islp, name)
  end
end
