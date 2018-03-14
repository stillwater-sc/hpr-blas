# HPR-BLAS
High-Performance Reproducible BLAS using posit arithmetic

# Introduction

Round-off makes floating point addition nonassociative, and different orders of summation often produce different results.
On a parallel machine, the order of summation will vary from run to run, or even subroutine-call to subroutine-call, depending on scheduling of available resources, so results will different from run to run.
HPR-BLAS offers high-performance reproducible BLAS independent of concurrency levels of the underlying hardware.
