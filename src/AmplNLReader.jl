# A rudimentary AMPL interface in Julia.
# D. Orban, Vancouver, April 2014.
module AmplNLReader

using NLPModels
using Compat

include(Pkg.dir("MathProgBase","src","NLP","NLP.jl"))
using .NLP  # Defines NLPModelMeta.

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
  include("../deps/deps.jl")
else
  error("ASL library not properly installed. Please run Pkg.build(\"AmplNLReader\")")
end

include("ampl_model.jl")


end  # Module AmplNLReader
