# A rudimentary AMPL interface in Julia.
# D. Orban, Vancouver, April 2014.
module AmplNLReader

using NLPModels
using Compat
import Compat.String

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
  include("../deps/deps.jl")
else
  error("ASL library not properly installed. Please run Pkg.build(\"AmplNLReader\")")
end

include("ampl_model.jl")

if Pkg.installed("MathProgBase") != nothing
  include("mpb_interface.jl")
end

end  # Module AmplNLReader
