# Test ampl.jl.
include("ampl.jl")

function exercise_ampl_model(nlp :: AmplModel)
  print(nlp)

  f = obj( nlp, nlp.x0)
  g = grad(nlp, nlp.x0)
  c = cons(nlp, nlp.x0)
  J = jac( nlp, nlp.x0)
  H = hess(nlp, nlp.x0, y=ones(nlp.ncon,))

  @printf("f(x0) = %f\n", f)
  @printf("∇f(x0) = "); display(g'); @printf("\n")
  @printf("c(x0) = ");  display(c'); @printf("\n")
  for j = 1 : nlp.ncon
    @printf("∇c_%d(x0) =", j)
    display(jth_congrad(nlp, nlp.x0, j)); @printf("\n")
    @printf("sparse ∇c_%d(x0) =", j)
    display(jth_sparse_congrad(nlp, nlp.x0, j)); @printf("\n")
  end
  @printf "J(x0) = \n";      display(J); @printf("\n")
  @printf "∇²L(x0,y0) = \n"; display(H); @printf("\n")

  e = ones(nlp.nvar)
  for j = 1 : nlp.ncon
    Hje = jth_hprod(nlp, nlp.x0, e, j)
    @printf("∇²c_%d(x0) * e = ", j); display(Hje'); @printf("\n")
  end

  ghje = ghjvprod(nlp, nlp.x0, g, e)
  @printf "(∇f(x0), ∇²c_j(x0) * e) = "; display(ghje'); @printf("\n")
end

hs33 = AmplModel("hs033.nl")
rosenbrock = AmplModel("rosenbr.nl")
exercise_ampl_model(rosenbrock)
hs9 = AmplModel("hs009.nl")
exercise_ampl_model(hs33)
exercise_ampl_model(hs9)
