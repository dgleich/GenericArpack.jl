#=
usage:

write the code you want to track allocations for into a single file, say "myfile.jl"

    using LinearAlgebra
    using Profile # to allow reseting allocation to get your target function 

    function myfun(n::Int)
      return  randn(n,n)
    end 
    Z = myfun(100)
    Profile.clear_malloc_data()
    Z = myfun(100)

then run 

    include("allocations.jl")
    lines = report_allocations("myfile.jl")
    println.(lines); 

Or the code is smart enough to be able to run from the same file.

    using LinearAlgebra
    using Profile # to allow reseting allocation to get your target function 

    function myfun(n::Int)
      return  randn(n,n)
    end 
    Z = myfun(100)
    Profile.clear_malloc_data()
    Z = myfun(100)

    ## 
    include("allocations.jl")
    lines = report_allocations(@__FILE__) # this will auto-exit 
    println.(lines);
=#  

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

function report_allocations(file_to_run; depotpaths=true, cleanup=true)
  # check if we are running from the file we are reporting on! 
  if haskey(ENV, "JULIA_REPORT_ALLOCS")
    exit(0)
  end 
  curenv = dirname(Base.active_project())
  jlpath = joinpath(Sys.BINDIR, "julia")
  cmd = Cmd(`$jlpath --project=$(curenv) --track-allocation=user $file_to_run`, env=("JULIA_REPORT_ALLOCS"=>"1",))
  # ripped from run(...) and adaopted
  ps = Base._spawn(cmd, Base.spawn_opts_inherit())
  pid = getpid(ps)
  if success(ps)
    alines = show_allocations("."; depotpaths, cleanup, pid)
    if length(alines) == 0
      println("No allocations found")
    end
  else
    # cleanup if error
    alines = show_allocations("."; depotpaths, cleanup, pid)
    Base.pipeline_error(ps)
  end 
  return alines
end 