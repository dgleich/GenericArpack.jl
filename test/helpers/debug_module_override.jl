#=
I had problems getting the dgetv0_override.jl functions working.
So I developed this little test to explore what happens.

It turns out this doesn't work with keyword arguments. <sigh>
This example illustrates that.

This motivates a refactor of the interface involving ArpackState 
where we don't make that a kw parameter and make it a fixed
parameter, which it really was anyway, except for dsaupd... (i.e.
an other call... even though that may also be a good idea to add
a forwarding function)
=#
module MyMod
  abstract type AbstractMyType end
  mutable struct MyType1 <: AbstractMyType
    val::Int
  end 
  function myfun(x::AbstractMyType)
    @show "In myfun with AbstractMyType $(x.val)"
    @show typeof(x.val) 
    return nothing 
  end 

  function mykwfun(z::Int; x::AbstractMyType)
    @show "In mykwfun with z=$(z) AbstractMyType $(x.val)"
    @show typeof(x.val) 
    return nothing 
  end 

  function myoptfun(z::Int, x::Union{AbstractMyType,Nothing} = nothing)
    @show "In myoptfun with z=$(z) AbstractMyType $(typeof(x))"
    if x !== nothing 
      @show typeof(x.val) 
    end
    return nothing 
  end 
end
@show methods(MyMod.myfun)
@show methods(MyMod.mykwfun)
MyMod.myfun(MyMod.MyType1(5))
MyMod.mykwfun(1; x=MyMod.MyType1(5))

module MyNewMod
  import ..MyMod
  mutable struct MyType2 <: MyMod.AbstractMyType
    val::Float64
  end 

  function MyMod.myfun(x::MyType2)
    @show "In myfun with MyType2"
    @show typeof(x.val) 
    MyMod.myfun(MyMod.MyType1(6))
    return nothing 
  end 
  
  function MyMod.mykwfun(z::Int; x::MyType2)
    @show "In mykwfun with z=$(z) MyType2 $(x.val)"
    @show typeof(x.val) 
    MyMod.mykwfun(2, MyMod.MyType1(6))
    return nothing 
  end 

  function MyMod.myoptfun(z::Int, x::MyType2)
    @show "In myoptfun with z=$(z) MyType2 $(x.val)"
    @show typeof(x.val) 
    @show "Forwaring..."
    MyMod.myoptfun(2, MyMod.MyType1(6))
    return nothing 
  end 
end 

@show methods(MyMod.myfun)



MyMod.myfun(MyNewMod.MyType2(5.0))
MyMod.myfun(MyMod.MyType1(5))

# This should still work, but it doesn't
@show methods(MyMod.mykwfun)
MyMod.mykwfun(1; x=MyMod.MyType1(5))

## It does work with opt parameters though
@show methods(MyMod.myoptfun)
MyMod.myoptfun(1, MyMod.MyType1(5))
MyMod.myoptfun(1, MyNewMod.MyType2(5))

