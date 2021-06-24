function _checked_dlaev2(a::T,b::T,c::T) where {T <: Float64}
  blam1,blam2,bc,bs = _dlaev2_blas(a,b,c)
  jlam1,jlam2,jc,js = _dlaev2(a,b,c)
  if (blam1,blam2,bc,bs) != (jlam1,jlam2,jc,js)
    @show a,b,c
    @show blam1,blam2,bc,bs
    @show jlam1,jlam2,jc,js
    error("Difference detected")
  end
  return jlam1,jlam2,jc,js
end



function _checked_dlasr_right_side_variable_pivot_backward!(m::Integer, n::Integer,
    c::AbstractVector{T}, s::AbstractVector{T}, a::AbstractVecOrMat{T}; extra=similar(a)) where T
  copyto!(extra, a)
  savecopy = copyto!(similar(a), a)
  _dlasr_right_side_variable_pivot_backward!(m, n, c, s, a)
  _dlasr_rvb_blas!(m, n, c, s, extra.parent, 1)
  ablas = extra
  if norm(extra - a,Inf) >= 1e-8
    @show "Difference on"
    @show m, n
    @show c
    @show s
    @show savecopy
    @show a
    @show ablas
    error("Difference detected")
  end
  return a
end
