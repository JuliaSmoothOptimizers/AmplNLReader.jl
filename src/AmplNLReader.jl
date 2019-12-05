# A rudimentary AMPL interface in Julia.
# D. Orban, Vancouver, April 2014.
module AmplNLReader

using LinearAlgebra
using NLPModels
using SparseArrays

using Libdl
using ASL_jll

const libasl = joinpath(dirname(ASL_jll.libasl_path), "libasl." * dlext)

include("ampl_model.jl")

end  # Module AmplNLReader
