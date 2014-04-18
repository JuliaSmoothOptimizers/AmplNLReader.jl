# Test ampl.jl.

include("ampl.jl")
include("ampl_utils.jl")

stub = "hs033.nl"
nlp = AmplModel(stub)

@printf "problem %s\n" nlp.name
@printf "nvar = %d, ncon = %d\n" nlp.nvar nlp.ncon
@printf "lvar = "; print_array(nlp.lvar)
@printf "uvar = "; print_array(nlp.uvar)
@printf "lcon = "; print_array(nlp.lcon)
@printf "ucon = "; print_array(nlp.ucon)
@printf "x0 = "; print_array(nlp.x0)
@printf "y0 = "; print_array(nlp.y0)

f = obj( nlp, nlp.x0)
g = grad(nlp, nlp.x0)
c = cons(nlp, nlp.x0)
J = jac( nlp, nlp.x0)
H = hess(nlp, nlp.x0, y=ones(nlp.ncon,))

@printf "f(x0) = %f\n" f
@printf "∇f(x0) = "; print_array(g)
@printf "c(x0) = "; print_array(c)
for j = 1 : nlp.ncon
  println("∇c_$j(x0) =")
  println(jth_sparse_congrad(nlp, nlp.x0, j))
end
@printf "J(x0) = \n"; println(J)
@printf "∇²L(x0,y0) = \n"; println(H)

e = ones(nlp.nvar)
for j = 1 : nlp.ncon
  print_formatted("∇²c_%d(x0) * e =", j);
  print_array(jth_hprod(nlp, nlp.x0, e, j))
end

ghje = ghjvprod(nlp, nlp.x0, g, e)
@printf "(∇f(x0), ∇²c_j(x0) * e) = "; print_array(ghje)
