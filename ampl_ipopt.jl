# Example script to solve AMPL model with IPOPT.
# Needs some restructuring to be callable from AMPL

require("ampl.jl")
require("ampl_mathprogbase.jl")

using Ipopt
import MathProgBase

function solve_with_ipopt(nlfile::ASCIIString)
  nlp = AmplModel(nlfile)
  # options can be set here
  m = MathProgBase.model(IpoptSolver())
  loadamplproblem!(m, nlp)
  MathProgBase.optimize!(m)

  objval = MathProgBase.getobjval(m)
  x = MathProgBase.getsolution(m)

  println("Optimal value: $objval")
  println("Solution: $x")
  
end

if length(ARGS) != 1
  error("Usage: julia ampl_ipopt.jl nlfile")
end

solve_with_ipopt(ascii(ARGS[1]))
