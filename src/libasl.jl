function asl_init(stub)
  @ccall libasl.asl_init(stub::Ptr{Cchar})::Ptr{Cvoid}
end

function asl_finalize(asl)
  @ccall libasl.asl_finalize(asl::Ptr{Cvoid})::Cvoid
end

function asl_write_sol(asl, msg, x, y)
  @ccall libasl.asl_write_sol(asl::Ptr{Cvoid}, msg::Ptr{Cchar}, x::Ptr{Cdouble},
                              y::Ptr{Cdouble})::Cvoid
end

function asl_objtype(asl)
  @ccall libasl.asl_objtype(asl::Ptr{Cvoid})::Cint
end

function asl_nlo(asl)
  @ccall libasl.asl_nlo(asl::Ptr{Cvoid})::Cint
end

function asl_nzo(asl)
  @ccall libasl.asl_nzo(asl::Ptr{Cvoid})::Cint
end

function asl_nvar(asl)
  @ccall libasl.asl_nvar(asl::Ptr{Cvoid})::Cint
end

function asl_nbv(asl)
  @ccall libasl.asl_nbv(asl::Ptr{Cvoid})::Cint
end

function asl_niv(asl)
  @ccall libasl.asl_niv(asl::Ptr{Cvoid})::Cint
end

function asl_nlvb(asl)
  @ccall libasl.asl_nlvb(asl::Ptr{Cvoid})::Cint
end

function asl_nlvo(asl)
  @ccall libasl.asl_nlvo(asl::Ptr{Cvoid})::Cint
end

function asl_nlvc(asl)
  @ccall libasl.asl_nlvc(asl::Ptr{Cvoid})::Cint
end

function asl_nlvbi(asl)
  @ccall libasl.asl_nlvbi(asl::Ptr{Cvoid})::Cint
end

function asl_nlvci(asl)
  @ccall libasl.asl_nlvci(asl::Ptr{Cvoid})::Cint
end

function asl_nlvoi(asl)
  @ccall libasl.asl_nlvoi(asl::Ptr{Cvoid})::Cint
end

function asl_nwv(asl)
  @ccall libasl.asl_nwv(asl::Ptr{Cvoid})::Cint
end

function asl_ncon(asl)
  @ccall libasl.asl_ncon(asl::Ptr{Cvoid})::Cint
end

function asl_nlc(asl)
  @ccall libasl.asl_nlc(asl::Ptr{Cvoid})::Cint
end

function asl_lnc(asl)
  @ccall libasl.asl_lnc(asl::Ptr{Cvoid})::Cint
end

function asl_nlnc(asl)
  @ccall libasl.asl_nlnc(asl::Ptr{Cvoid})::Cint
end

function asl_nnzj(asl)
  @ccall libasl.asl_nnzj(asl::Ptr{Cvoid})::Cint
end

function asl_nnzh(asl)
  @ccall libasl.asl_nnzh(asl::Ptr{Cvoid})::Cint
end

function asl_islp(asl)
  @ccall libasl.asl_islp(asl::Ptr{Cvoid})::Cint
end

function asl_n_cc(asl)
  @ccall libasl.asl_n_cc(asl::Ptr{Cvoid})::Cint
end

function asl_x0(asl)
  @ccall libasl.asl_x0(asl::Ptr{Cvoid})::Ptr{Cdouble}
end

function asl_y0(asl)
  @ccall libasl.asl_y0(asl::Ptr{Cvoid})::Ptr{Cdouble}
end

function asl_lvar(asl)
  @ccall libasl.asl_lvar(asl::Ptr{Cvoid})::Ptr{Cdouble}
end

function asl_uvar(asl)
  @ccall libasl.asl_uvar(asl::Ptr{Cvoid})::Ptr{Cdouble}
end

function asl_lcon(asl)
  @ccall libasl.asl_lcon(asl::Ptr{Cvoid})::Ptr{Cdouble}
end

function asl_ucon(asl)
  @ccall libasl.asl_ucon(asl::Ptr{Cvoid})::Ptr{Cdouble}
end

function asl_cvar(asl)
  @ccall libasl.asl_cvar(asl::Ptr{Cvoid})::Ptr{Cint}
end

function asl_varscale(asl, s, err)
  @ccall libasl.asl_varscale(asl::Ptr{Cvoid}, s::Ptr{Cdouble}, err::Ptr{Cint})::Cvoid
end

function asl_lagscale(asl, s, err)
  @ccall libasl.asl_lagscale(asl::Ptr{Cvoid}, s::Cdouble, err::Ptr{Cint})::Cvoid
end

function asl_conscale(asl, s, err)
  @ccall libasl.asl_conscale(asl::Ptr{Cvoid}, s::Ptr{Cdouble}, err::Ptr{Cint})::Cvoid
end

function asl_obj(asl, x, err)
  @ccall libasl.asl_obj(asl::Ptr{Cvoid}, x::Ptr{Cdouble}, err::Ptr{Cint})::Cdouble
end

function asl_grad(asl, x, g, err)
  @ccall libasl.asl_grad(asl::Ptr{Cvoid}, x::Ptr{Cdouble}, g::Ptr{Cdouble}, err::Ptr{Cint})::Cvoid
end

function asl_cons(asl, x, c, err)
  @ccall libasl.asl_cons(asl::Ptr{Cvoid}, x::Ptr{Cdouble}, c::Ptr{Cdouble}, err::Ptr{Cint})::Cvoid
end

function asl_jcon(asl, x, j, err)
  @ccall libasl.asl_jcon(asl::Ptr{Cvoid}, x::Ptr{Cdouble}, j::Cint, err::Ptr{Cint})::Cdouble
end

function asl_jcongrad(asl, x, g, j, err)
  @ccall libasl.asl_jcongrad(asl::Ptr{Cvoid}, x::Ptr{Cdouble}, g::Ptr{Cdouble}, j::Cint,
                             err::Ptr{Cint})::Cvoid
end

function asl_hprod(asl, y, v, hv, w)
  @ccall libasl.asl_hprod(asl::Ptr{Cvoid}, y::Ptr{Cdouble}, v::Ptr{Cdouble}, hv::Ptr{Cdouble},
                          w::Cdouble)::Cvoid
end

function asl_hvcompd(asl, v, hv, nobj)
  @ccall libasl.asl_hvcompd(asl::Ptr{Cvoid}, v::Ptr{Cdouble}, hv::Ptr{Cdouble}, nobj::Cint)::Cvoid
end

function asl_ghjvprod(asl, g, v, ghjv)
  @ccall libasl.asl_ghjvprod(asl::Ptr{Cvoid}, g::Ptr{Cdouble}, v::Ptr{Cdouble},
                             ghjv::Ptr{Cdouble})::Cvoid
end

function asl_sparse_congrad_nnz(asl, j)
  @ccall libasl.asl_sparse_congrad_nnz(asl::Ptr{Cvoid}, j::Cint)::Csize_t
end

function asl_sparse_congrad(asl, x, j, inds, vals, err)
  @ccall libasl.asl_sparse_congrad(asl::Ptr{Cvoid}, x::Ptr{Cdouble}, j::Cint, inds::Ptr{Cint},
                                   vals::Ptr{Cdouble}, err::Ptr{Cint})::Cvoid
end

function asl_jac(asl, x, rows, cols, vals, err)
  @ccall libasl.asl_jac(asl::Ptr{Cvoid}, x::Ptr{Cdouble}, rows::Ptr{Cint}, cols::Ptr{Cint},
                        vals::Ptr{Cdouble}, err::Ptr{Cint})::Cvoid
end

function asl_jac_structure(asl, rows, cols)
  @ccall libasl.asl_jac_structure(asl::Ptr{Cvoid}, rows::Ptr{Cint}, cols::Ptr{Cint})::Cvoid
end

function asl_jacval(asl, x, vals, err)
  @ccall libasl.asl_jacval(asl::Ptr{Cvoid}, x::Ptr{Cdouble}, vals::Ptr{Cdouble},
                           err::Ptr{Cint})::Cvoid
end

function asl_hess(asl, y, w, rows, cols, vals)
  @ccall libasl.asl_hess(asl::Ptr{Cvoid}, y::Ptr{Cdouble}, w::Cdouble, rows::Ptr{Cint},
                         cols::Ptr{Cint}, vals::Ptr{Cdouble})::Cvoid
end

function asl_hess_structure(asl, rows, cols)
  @ccall libasl.asl_hess_structure(asl::Ptr{Cvoid}, rows::Ptr{Cint}, cols::Ptr{Cint})::Cvoid
end

function asl_hessval(asl, y, w, vals)
  @ccall libasl.asl_hessval(asl::Ptr{Cvoid}, y::Ptr{Cdouble}, w::Cdouble, vals::Ptr{Cdouble})::Cvoid
end
