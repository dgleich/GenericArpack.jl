function _eps23(::Type{T}) where {T <: AbstractFloat}
  return (eps(T)/2)^((2one(T))/(3one(T)))
end

function _eps2(::Type{T}) where {T <: AbstractFloat}
  return eps(T)*eps(T)
end


""" Print out the fields for a specific type in Markdown format.

This is similar to what is in the Base.Docs module to automatically
summarize a type without any documentation, but allows me to use it
and lightly customize it. Also works to document parametric types.

See <https://github.com/JuliaLang/julia/blob/6f71de4193717f9267e98b4188aa6e4547146375/stdlib/REPL/src/docview.jl>
(Line 282)
"""
function _output_fields_in_markdown(T::Type)
  pad = maximum(length(string(f)) for f in fieldnames(T))
  return string("# Fields\n",
  "```\n",
    map( (Z) -> string(rpad(Z[1], pad), " :: ", Z[2], "\n"), zip(fieldnames(T), T.types) )...,
  "```\n")
end

#=
This file has simple macros to replicate the stats and timing functionality
from the ARPACK macros debug.h and stats.h. These define _common_ variables
in Fortran, which work like Globals. But they are not thread-safe, so we want
to improve the design in Julia. Consequently, we use julia's keyword arguments,
see note-03-stats to explain

c
c\SCCS Information: @(#)
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer  logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
c\SCCS Information: @(#)
c FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2
c
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
c
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/
     &           nopx, nbx, nrorth, nitref, nrstrt,
     &           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &           tmvopx, tmvbx, tgetv0, titref, trvec

=#

""" A debugging structure for Arpack computations.

Usage:
=====
- `debug=ArpackDebug()` creates a structure to print output to stdout, but does not enable
  any printing.
- `debug=ArpackDebug(); debug.maupd = 1` prints minimal information for the aupd routine.
- `debug=ArpackDebug(); debug.mapps = 4` prints maximal information for the deflation routine.

$(_output_fields_in_markdown(ArpackDebug{IO}))

There is help on each field with `? ArpackDebug.maupd` for example.
"""
Base.@kwdef mutable struct ArpackDebug{IOT <: IO}
  #= original FORTRAN
  integer  logfil, ndigit, mgetv0,
 &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
 &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
 &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
 =#
  """ The output logfile, defaults to stdout and must be created at start. """
  logfile::IOT = Base.stdout
  """ if ndigit < 0, then we assume 80-cols, otherwise 132-cols """
  ndigit::Int = (Base.displaysize(Base.stdout)[2] > 80 ? 1 : -1)

  """ TODO """
  mgetv0::Int = 0
  maupd::Int = 0
  maup2::Int = 0
  maitr::Int = 0
  meigt::Int = 0
  mapps::Int = 0
  mgets::Int = 0
  meupd::Int = 0
end

function set_debug_high!(debug::ArpackDebug)
  debug.mgetv0 = 4
  debug.maupd = 4
  debug.maup2 = 4
  debug.maitr = 4 
  debug.meigt = 4
  debug.mapps = 4
  debug.mgets = 4
  debug.meupd = 4
  return debug
end 

#ArpackDebug(;) = ArpackDebug(;logfile=stdout,ndigit=(displaysize(stdout)[2] > 80) ? 1 : -1)

ArpackTime = Float64

macro jl_arpack_time()
  return esc(:( stats === nothing ? zero(ArpackTime) : time() ))
end

macro jl_arpack_set_stat(field, value)
  return esc( quote
    if stats !== nothing
      stats.$field = $value
    end
  end )
end

macro jl_arpack_increment_stat(field)
  return esc( quote
    if stats !== nothing
      stats.$field += 1
    end
  end )
end

macro jl_update_time(field, t0 )
  return esc( quote
    if stats !== nothing
      stats.$field += time() - $t0
    end
  end )
end

macro jl_arpack_debug(field,default)
  return esc(:( debug === nothing ? $default : debug.$field ))
end

""" A timing structure for Arpack computations

$(_output_fields_in_markdown(ArpackStats))
"""
Base.@kwdef mutable struct ArpackStats
  #= Here are the original fortran
  integer    nopx, nbx, nrorth, nitref, nrstrt
  real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
 &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
 &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
 &           tmvopx, tmvbx, tgetv0, titref, trvec
 =#
  nopx::Int = 0
  nbx::Int = 0
  nrorth::Int = 0
  nitref::Int = 0
  nrstrt::Int = 0
  taupd::ArpackTime = 0.0  # removed the type name, so no "s/n/c"
  taup2::ArpackTime = 0.0
  taitr::ArpackTime = 0.0
  teigt::ArpackTime = 0.0
  tgets::ArpackTime = 0.0
  tapps::ArpackTime = 0.0
  tconv::ArpackTime = 0.0
  tmvopx::ArpackTime = 0.0
  tmvbx::ArpackTime = 0.0
  tgetv0::ArpackTime = 0.0
  titref::ArpackTime = 0.0
  trvec::ArpackTime = 0.0
end

macro jl_arpack_check_bmat(var)
  return esc( quote
    if ($var != :I) && ($var != :G)
      throw(ArgumentError(string("bmat is ", $var, " which is not :I or :G")))
    end 
  end)
end

macro jl_arpack_check_length(var,len)
  varstr = string(var)
  lenstr = string(len)
  errstr1 = "range 1:$lenstr="
  errstr2 = " out of bounds for $varstr of length "
  return esc( quote
    if $len > length($var)
      throw(ArgumentError(string($errstr1, $len, $errstr2, length($var))))
    end
  end)
end

macro jl_arpack_check_size(var,m,n)
  varstr = string(var)
  mstr = string(m)
  nstr = string(n)
  errstr1 = "size $mstr x $nstr = ("
  errstr2 = ") out of bounds for $varstr of size "
  return esc( quote
    if $m > size($var,1) || $n > size($var, 2)
      throw(ArgumentError(string($errstr1, $m, ", ", $n, $errstr2, size($var))))
    end
  end)
end

## State saving variables

function _attach_state(statevar,vars...)
  # we want to build up a block of expressions.
  block = Expr(:block)
  for f in vars
      # each expression in the block consists of
      # the fieldname = struct.fieldname
      e = :($f = $statevar.$f)
      # add this new expression to our block
      push!(block.args, e)
  end
  # now escape the evaled block so that the
  # new variable declarations get declared in the surrounding scope.
  return block
end

##
Base.@kwdef struct AitrState{T}
  orth1::Bool = false
  orth2::Bool = false
  rstart::Bool = false
  step3::Bool = false
  step4::Bool = false
  ierr::Int = 0 # TODO check and remove... 
  ipj::Int = 0 # TODO remove
  irj::Int = 0 # TODO remove
  ivj::Int = 0 # TODO remove 
  iter::Int = 0
  itry::Int = 0 
  j::Int = 0 
  rnorm1::T = zero(T) # TODO remove 
  wnorm::T = zero(T)
  t0::ArpackTime = zero(ArpackTime)
  t2::ArpackTime = zero(ArpackTime)
  t4::ArpackTime = zero(ArpackTime)
end

const _aitr_state_vars = (
        :orth1, :orth2, :rstart, :step3, :step4,   # logical vars
        :ierr, :ipj, :irj, :ivj, :iter, :itry, :j, # integer vars
        :rnorm1, :wnorm, # double vars
        :t0, :t2, :t4
        )

macro attach_aitr_state(statevar)
  expr = _attach_state(:($statevar.aitr), _aitr_state_vars...)
  return esc(:($expr))
end

macro aitr_state_vars()
  e = Expr(:parameters, _aitr_state_vars...)
  return esc( :($e) )
end

##

Base.@kwdef struct Getv0State{T}
  iseed::Base.RefValue{NTuple{4,Int64}} = Base.Ref(tuple(1,3,5,7))
  first::Bool = false
  orth::Bool = false
  iter::Int = 0
  rnorm0::T = zero(T)
  t0::ArpackTime = zero(ArpackTime)
  t2::ArpackTime = zero(ArpackTime) 
end

const _getv0_state_vars = (
        :iseed,
        #:iseed1, :iseed2, :iseed3, :iseed4, # seed vars
        :first, :orth, # logical state vars
        :iter, :rnorm0, # status vars
        :t0, :t2, # time vars
        )

macro attach_getv0_state(statevar)
  expr = _attach_state(:($statevar.getv0), _getv0_state_vars...)
  return esc(:($expr))
end

macro getv0_state_vars()
  e = Expr(:parameters, _getv0_state_vars...)
  return esc( :($e) )
end


## aup2 state
Base.@kwdef struct Saup2State{T}
  cnorm::Bool = false
  getv0::Bool = false
  initv::Bool = false 
  update::Bool = false
  ushift::Bool = false
  iter::Int = 0
  kplusp::Int = 0
  nconv::Int = 0
  nev0::Int = 0 
  np0::Int = 0 
  #rnorm::T = zero(T) # needs to be a ref as it's mutated by other funcs
  t0::ArpackTime = zero(ArpackTime)
  t2::ArpackTime = zero(ArpackTime) 
end

const _saupd2_state_vars = (
        :cnorm, :getv0, :initv, :update, :ushift, 
        :iter, :kplusp, :nconv, :nev0, :np0, 
        :rnorm, 
        :t0, :t2, # time vars
        )

macro attach_saup2_state(statevar)
  expr = _attach_state(:($statevar.saup2), _saupd2_state_vars...)
  return esc(:($expr))
end

macro saup2_state_vars()
  e = Expr(:parameters, _saupd2_state_vars...)
  return esc( :($e) )
end


## aupd state 
## overall state

abstract type AbstractArpackState{T} end

Base.@kwdef mutable struct ArpackState{T} <: AbstractArpackState{T}
  aitr::AitrState{T} = AitrState{T}()
  getv0::Getv0State{T} = Getv0State{T}()
  saup2::Saup2State{T} = Saup2State{T}()
  aupd_nev0 = Ref{Int}(0)
  aupd_np = Ref{Int}(0)
  aupd_mxiter = Ref{Int}(0)
  aup2_rnorm = Ref{T}(zero(T))
end



