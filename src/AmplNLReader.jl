# A rudimentary AMPL interface in Julia.
# D. Orban, Vancouver, April 2014.
module AmplNLReader

using NLPModels
using SparseArrays

import ASL_jll

const libasl = ASL_jll.libasl

include("libasl.jl")
include("ampl_meta.jl")
include("ampl_model.jl")
include("ampl_mpec_model.jl")

end  # Module AmplNLReader
