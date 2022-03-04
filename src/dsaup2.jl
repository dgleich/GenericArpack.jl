function dsaup2()
  successexit = false # label 1100
  errorexit = false # label 1200
  startiter = false # label 1000
  while true
    if dgetv0
      # execution will fall through to the "else" block
      # below to more closely mirror what happens
      # in arpack. 
    elseif startiter 
      startiter = false
      iter += 1
      if msglvl > 0
      end
      if msglvl > 1
      end
      ido[] = 0
      update = true 
    elseif update
      # no need to set update to true...
      @debug "label 20"
      dsaitr()
      if ido[] != 99
        break # go to 9000
      end
      update = false 
      # lots of stuff...
      # setup next steps in ushift Block
      # in Julia, we always set ushift = true
      # to get to the next block of code.
      ushift = true
      if ishift == 0 
        ido[] = 3
        break
      end
    elseif ushift
      ushift = false 

      cnorm = true
      # setup cnorm call... 
    elseif cnorm
      # c        | Back from reverse communication; |
      # c        | WORKD(1:N) := B*RESID            |
      cnorm = false 

      startiter = true # next iteration! 
    else
      # this is  
      # init lanczos... 
      @debug "initalizing the iteration"
      startiter = true 
    end
  end

  if successexit || errorexit
    if successexit
      mxiter = iter
      nev = nconv
    end
    ido[] = 99
    @jl_update_time(taup2, t0)
  end
  
  return info
  
end