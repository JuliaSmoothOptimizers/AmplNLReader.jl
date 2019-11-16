using Test
using NLPModels
using AmplNLReader
using LinearAlgebra
using Printf
using SparseArrays

function exercise_ampl_model(nlp :: AmplModel)
  show(stdout, nlp)
  print(stdout, nlp)

  # Perform dummy scaling.
  dummy_varscaling = ones(2 * nlp.meta.nvar)
  dummy_conscaling = ones(2 * nlp.meta.ncon)
  varscale(nlp, @view dummy_varscaling[1:2:end])  # to exercise AbstractArray input
  conscale(nlp, @view dummy_conscaling[1:2:end])
  lagscale(nlp, 1.0)

  f = obj( nlp, nlp.meta.x0)
  g = grad(nlp, nlp.meta.x0)
  c = cons(nlp, nlp.meta.x0)
  J = jac( nlp, nlp.meta.x0)

  x = Vector{Cdouble}(undef, 2 * nlp.meta.nvar)
  x[1:2:end] .= nlp.meta.x0
  xview = @view x[1:2:end]

  for j = 1 : nlp.meta.ncon
    cj = jth_con(nlp, xview, j)
    @test cj == c[j]
    gj = jth_congrad(nlp, xview, j)
    @test all(gj .== Vector{Float64}(J[j, :]))
    sgj = jth_sparse_congrad(nlp, xview, j)
    @test all(gj .== Vector{Float64}(sgj))
  end

  jrows, jcols = jac_structure(nlp)
  jvals = Vector{Float64}(undef, nlp.meta.nnzj)
  jac_coord!(nlp, xview, jvals)
  @test norm(sparse(jrows, jcols, jvals, nlp.meta.ncon, nlp.meta.nvar) - J) ≤ sqrt(eps()) * norm(J)

  jvals2 = Vector{Float64}(undef, 2 * nlp.meta.nnzj)
  jac_coord!(nlp, xview, @view jvals2[1:2:end])
  @test all(jvals .== jvals2[1:2:end])

  e2 = ones(nlp.meta.ncon)
  H = hess(nlp, xview, e2)
  hrows, hcols = hess_structure(nlp)
  hvals = Vector{Float64}(undef, nlp.meta.nnzh)
  hess_coord!(nlp, xview, e2, hvals)
  @test norm(sparse(hrows, hcols, hvals, nlp.meta.nvar, nlp.meta.nvar) - H) ≤ sqrt(eps()) * norm(H)

  hvals2 = Vector{Float64}(undef, 2 * nlp.meta.nnzh)
  hess_coord!(nlp, xview, e2, @view hvals2[1:2:end])
  @test all(hvals .== hvals2[1:2:end])

  H = hess(nlp, xview)
  hvals = Vector{Float64}(undef, nlp.meta.nnzh)
  hess_coord!(nlp, xview, hvals)
  @test norm(sparse(hrows, hcols, hvals, nlp.meta.nvar, nlp.meta.nvar) - H) ≤ sqrt(eps()) * norm(H)

  hvals2 = Vector{Float64}(undef, 2 * nlp.meta.nnzh)
  hess_coord!(nlp, xview, @view hvals2[1:2:end])
  @test all(hvals .== hvals2[1:2:end])

  e = ones(nlp.meta.nvar)
  je = jprod(nlp, xview, e)
  @test norm(je - J * e) ≤ sqrt(eps()) * norm(je)

  jte = jtprod(nlp, xview, e2)
  @test norm(jte - J' * e2) ≤ sqrt(eps()) * norm(jte)

  je2 = Vector{Float64}(undef, 2 * nlp.meta.ncon)
  jprod!(nlp, xview, e, @view je2[1:2:end])
  @test all(je .== je2[1:2:end])

  jte2 = Vector{Float64}(undef, 2 * nlp.meta.nvar)
  jtprod!(nlp, xview, e2, @view jte2[1:2:end])
  @test all(jte .== jte2[1:2:end])

  he = hprod(nlp, xview, e)
  @test norm(he - Symmetric(hess(nlp, xview), :L) * e) ≤ sqrt(eps()) * norm(he)

  he2 = Vector{Float64}(undef, 2 * nlp.meta.nvar)
  hprod!(nlp, xview, e, @view he2[1:2:end])
  @test all(he .== he2[1:2:end])

  he2 = Vector{Float64}(undef, 2 * nlp.meta.nvar)
  jth_hprod!(nlp, xview, e, 0, @view he2[1:2:end])  # same as objective hprod
  @test all(he .== he2[1:2:end])

  ghje = Vector{Float64}(undef, nlp.meta.ncon)
  y = zeros(nlp.meta.ncon)
  for j = 1 : nlp.meta.ncon
    Hje = jth_hprod(nlp, xview, e, j)
    y[j] = 1
    H = hess(nlp, xview, y, obj_weight=0.0)
    @test norm(Hje - H * e) ≤ sqrt(eps()) * norm(Hje)
    y[j] = 0
    ghje[j] = dot(g, Hje)
  end

  ghje2 = ghjvprod(nlp, xview, g, e)
  @test norm(ghje - ghje2) ≤ sqrt(eps()) * norm(ghje)

  ghje = Vector{Float64}(undef, 2 * nlp.meta.ncon)
  ghjvprod!(nlp, xview, g, e, @view ghje[1:2:end])
  @test all(ghje2 .== ghje[1:2:end])

  write_sol(nlp, "And the winner is...", rand(nlp.meta.nvar), rand(nlp.meta.ncon))
  reset!(nlp)
end

path = dirname(@__FILE__)
hs33 = AmplModel(joinpath(path, "hs033.nl"))
rosenbrock = AmplModel(joinpath(path, "rosenbr.nl"))
exercise_ampl_model(rosenbrock)
hs9 = AmplModel(joinpath(path, "hs009.nl"))
exercise_ampl_model(hs33)
exercise_ampl_model(hs9)
amplmodel_finalize(hs33)
amplmodel_finalize(rosenbrock)
amplmodel_finalize(hs9)
@test amplmodel_finalize(hs9) === nothing
@test_throws AmplException obj(hs9, hs9.meta.x0)
@test_throws AmplException AmplModel("this_file_does_not_exist")
@test_throws AmplException AmplModel("this_file_does_not_exist.nl")

# check maximization problems are handled correctly, i.e.,
# we do not explicitly change the sign of the objective
# hs6max is the same as hs6 except we maximize the objective
hs6min = AmplModel(joinpath(path, "hs6.nl"))
hs6max = AmplModel(joinpath(path, "hs6max.nl"))
x = hs6min.meta.x0
@test obj(hs6min, x) ≈ obj(hs6max, x)
@test all(grad(hs6min, x) .≈ grad(hs6max, x))
@test all(Matrix(hess(hs6min, x)) .≈ Matrix(hess(hs6max, x)))

