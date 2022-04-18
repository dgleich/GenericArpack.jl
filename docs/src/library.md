# Library

## User Friendly Interfaces
```@docs
symeigs
svds
```

## ArpackOp types

### User-facing ArpackOp types
```@docs
ArpackSimpleOp
ArpackSimpleFunctionOp
ArpackSymmetricGeneralizedOp
#ArpackSymmetricGeneralizedFunctionOp
ArpackNormalOp
ArpackNormalFunctionOp
```

### In-development ArpackOp types
```@docs
ArpackShiftInvertOp
ArpackBucklingOp
```

### Abstract ArpackOp Types
```@docs
ArpackOp
```

```@docs
ArpackSVDOp
```

## Helpers
```@docs
GenericArpack.svd_residuals
```

## Arpack drivers
```@docs
GenericArpack.dsaupd!
GenericArpack.simple_dseupd!
```