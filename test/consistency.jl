@testset "Consistency" begin
  problems = ["BROWNDEN", "HS5", "HS6", "HS10", "HS11", "HS14"]
  for problem in problems
    @testset "Problem $problem" begin
      nlp_ampl = AmplModel(joinpath("problems", lowercase(problem) * ".nl"), safe=true)
      nlp_man = eval(Symbol(problem))()
      nlps = [nlp_ampl, nlp_man]

      consistent_nlps(nlps)
    end
  end
end