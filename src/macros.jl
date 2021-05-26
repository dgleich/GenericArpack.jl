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


ArpackTime = Float64

macro jl_arpack_time()
  return esc(:( stats == nothing ? zero(ArpackTime) : time() ))
end

macro jl_update_time(field, t0 )
  return esc( quote
    if stats != nothing
      stats.$field += time() - $t0
    end
  end )
end

macro jl_arpack_debug(field,default)
  return esc(:( debug == nothing ? default : debug.$field ))
end

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
