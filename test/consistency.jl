using JuMP, NLPModels

nlppath = joinpath(Pkg.dir("NLPModels"), "test")
include(joinpath(nlppath, "consistency.jl"))
testpath = dirname(@__FILE__)

problems = [:brownden, :hs5, :hs6, :hs10, :hs11, :hs14]
for problem in problems
  problem_s = string(problem)
  include(joinpath(nlppath, "$problem_s.jl"))

  problem_f = eval(problem)
  nlp_ampl = AmplModel(joinpath(testpath, "$problem_s.nl"), safe=true)
  nlp_mpb = MathProgNLPModel(problem_f())
  nlps = [nlp_ampl, nlp_mpb]

  @printf("Checking problem %-15s%12s\t", problem_s, "")
  consistent_nlps(nlps)
end
