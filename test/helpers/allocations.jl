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
function show_allocations(maindir::AbstractString; pid=nothing, depotpaths=false)
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
  alines 
end 
#show_allocations("."; depotpaths=true, pid=76133)
