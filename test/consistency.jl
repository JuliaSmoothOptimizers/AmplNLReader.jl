using NLPModels

nlppath = joinpath(dirname(pathof(NLPModels)), "..", "test")
include(joinpath(nlppath, "consistency.jl"))
testpath = dirname(@__FILE__)

problems = [:brownden, :hs5, :hs6, :hs10, :hs11, :hs14]
for problem in problems
  problem_s = string(problem)
  include(joinpath(nlppath, "$problem_s.jl"))

  nlp_ampl = AmplModel(joinpath(testpath, "$problem_s.nl"), safe=true)
  nlp_autodiff = eval(Meta.parse("$(problem)_autodiff"))()
  nlps = [nlp_ampl, nlp_autodiff]

  @printf("Checking problem %-15s%12s\t", problem_s, "")
  consistent_nlps(nlps)
end
