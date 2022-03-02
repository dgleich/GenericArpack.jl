""" This function differs from the Arpack original
so that we can cut and paste the resulting output into
a julia vector. """
function _arpack_vout(io::IOT,
  n::Int, idigit::Int, v, msg::String) where {IOT}
  println(io)
  println(io, " ", msg)
  println(io, " ", "-"^(min(80, length(msg))))

  ncols = idigit < 0 ? 132 : 80
  ndigit = abs(idigit)  # TODO actually support ndigit

  nvals = idigit < 0 ? 5 : 2

  ilen=ceil(Int, log10(n+1)) # need log10(n+1) to handle n=100, say
  dlen=24

  for i=1:nvals:n
    k1=i
    k2=min(n,i+nvals-1)
    print(io, " ")
    for k=k1:k2
      print(io, rpad(v[k], dlen), ", ")
    end
    if k2==n
      # on the last entry, check if we need to pad 
      # for right alignment of comments 
      if k2-k1 < nvals
        for extra in nvals - (k2-k1)
          print(io, rpad("", dlen), "  ")
        end
      end
    end
    print(io, " # ", rpad(k1,ilen), "-", lpad(k2,ilen))
    println(io)
  end
end
function _arpack_vout(debug::ArpackDebug, msg::String, v)
  _arpack_vout(debug.logfile, length(v), debug.ndigit, v, msg)
end
