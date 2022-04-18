push!(LOAD_PATH,"../src/")
using Documenter, GenericArpack, LinearAlgebra, SparseArrays, StableRNGs
makedocs(
  sitename="GenericArpack.jl",
  doctest=false,
  authors="David F. Gleich and the ARPACK authors", 
  pages = Any[
    "Home" => "index.md", 
    "Library" => "library.md"
  ]
)

deploydocs(
  repo   = "github.com/dgleich/GenericArpack.jl.git"
)