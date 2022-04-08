module ArpackInJulia

import LinearAlgebra

include("macros.jl")
export ArpackDebug, ArpackStats, ArpackTime, ArpackState, set_debug_high!

include("idonow_ops.jl")
export ArpackOp, ArpackSimpleOp, ArpackSimpleFunctionOp, ArpackSymmetricGeneralizedOp

include("output.jl")
include("simple.jl")

include("arpack-blas.jl")
#include("arpack-blas-direct-temp.jl")
include("arpack-blas-qr.jl")
include("dstqrb.jl")

include("dgetv0.jl")
include("dsaitr.jl")

include("dsapps.jl")
include("dsaup2.jl")
include("dsaupd.jl")
include("dseupd.jl")

#include("arpack_jl_interface.jl")
include("interface.jl")
export eigs, symeigs

end # module
