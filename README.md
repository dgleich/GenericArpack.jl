ArpackInJulia
=============

This is a "me-project" to work on to I can learn something and see how various ideas work. 

The goal of this exercise is to port the double-precision ARPACK
for symmetric matrices in Julia. Including all ARPACK stuff. So this should
give "exactly" what ARPACK does but be a pure Julia implementation.
(Where exactly is ... it should be executing roughly the same sequence of
floating point operations and can differ on levels that would be expected
for different compilers compiling the same code.)

- not a goal to "Julia-ize" the package; I want to keep as close to the FORTRAN
  as possible so that I might be able to replace calls to Julia's Arpack.saupd /
  Arpack.seupd (which call the Fortran library) with this code.
- small internal function changes are okay, e.g. ARPACK has various debugging
  and timing stuff that would need to be done differently in Julia.
- small simplifications, e.g. if a function computes a single Int, we can
  rewrite that to return the Int rather than writing it into an array like in
  FORTRAN.
- Why? Why not implement my own ImplicitRestart/Eigensolver/Etc.? Simple: I trust
  ARPACK. Also, I want to understand exactly what the symmetric ARPACK solver is doing.
- Why not use a Fortran->Julia compiler? Well, I could. But I could also do 
  this and learn intimate details of how it works to help out in teaching :)
- I want to document some of the cool stuff in ARPACK!

Directories
-----------
- `notes` my notes as I'm going along / draft blog posts
- `dev` codes I hack together along the way. Usually my workflow is ... `dev -> src` where
  `src` is the final stuff that might be usable by others.
- `test` unit tests of the stuff in src
