function check_file_allocation(filename::AbstractString; alines=Vector{Pair{String, String}}())
  lines = readlines(filename)
  for line in lines
    line = lstrip(line)
    if startswith(line, "- ")
      # not run...
    elseif startswith(line, "0 ", )
      # not allocating
    else
      push!(alines, (filename => String(line)))
    end 
  end 
  return alines
end 
""" Add files matching pidstring to the allfiles list """
function _make_filelist!(maindir, allfiles, pidstring)
  if isdir(maindir)
    allfileinfo = walkdir(maindir)
    for (root, dirs, files) in allfileinfo
      for file in files
        if endswith(file, pidstring)
          push!(allfiles, joinpath(root, file))
        end
      end
    end
  end
end 
# TODO, add .julia/julia main dir list
function show_allocations(maindir::AbstractString; pid=nothing, depotpaths=false, cleanup=true)
  # make a list of all files
  
  allfiles = Vector{String}()
  pidstring = ".mem"
  if pid !== nothing
    pidstring = "$(pid).mem"
  end 
  _make_filelist!(maindir, allfiles, pidstring)
  if depotpaths
    for depotdir in DEPOT_PATH
      _make_filelist!(depotdir, allfiles, pidstring)
    end 
  end
  alines = Vector{Pair{String, String}}()
  for afile in allfiles
    check_file_allocation(afile; alines)
  end 
  if cleanup
    for afile in allfiles
      rm(afile)
    end 
  end 
  alines 
end 
#show_allocations("."; depotpaths=true, pid=76133)

function report_allocations(file_to_run; depotpaths=true, cleanup=true)
  # check if we are running from the file we are reporting on! 
  if haskey(ENV, "JULIA_REPORT_ALLOCS")
    exit(0)
  end 
  curenv = dirname(Base.active_project())
  jlpath = joinpath(Sys.BINDIR, "julia")
  cmd = Cmd(`$jlpath --project=$(curenv) --track-allocation=user $file`, env=("JULIA_REPORT_ALLOCS"=>"1",))
  # ripped from run(...) and adaopted
  ps = Base._spawn(cmd, Base.spawn_opts_inherit())
  pid = getpid(ps)
  if success(ps)
    alines = show_allocations("."; depotpaths, cleanup, pid)
  else
    Base.pipeline_error(ps)
  end 
end 