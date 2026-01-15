# Script to parse AMPL headers and generate Julia wrappers.
using ASL_jll
using Clang
using Clang.Generators
using JuliaFormatter

function main()
  cd(@__DIR__)
  include_dir = joinpath(ASL_jll.artifact_dir, "include")
  headers = map(header -> joinpath(include_dir, header), ["aslinterface.h"])

  options = load_options(joinpath(@__DIR__, "asl.toml"))
  options["general"]["output_file_path"] = joinpath("..", "src", "libasl.jl")
  options["general"]["output_ignorelist"] = ["ASL"]

  args = get_default_args()
  push!(args, "-I$include_dir")

  ctx = create_context(headers, args, options)
  build!(ctx)

  path = options["general"]["output_file_path"]
  code = read(path, String)
  code = replace(code, "Ptr{ASL}" => "Ptr{Cvoid}")
  write(path, code)
  format_file(path, YASStyle())
  return nothing
end

# If we want to use the file as a script with `julia wrapper.jl`
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
