module GenericArpack

import LinearAlgebra

include("macros.jl")
export ArpackDebug, ArpackStats, ArpackTime, ArpackState, set_debug_high!

include("idonow_ops.jl")
export ArpackOp, ArpackSVDOp
export ArpackSimpleOp, ArpackSimpleFunctionOp
export ArpackSymmetricGeneralizedOp, ArpackSymmetricGeneralizedFunctionOp
export ArpackShiftInvertOp, ArpackBucklingOp
export ArpackNormalOp, ArpackNormalFunctionOp

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
export eigs, symeigs, hermeigs, svds, complexsvds

include("fixes.jl")
#export @fix_doublefloats, @fix_multifloats
# we don't export these now, so it's clear where they come from in downstream code.

end # module
