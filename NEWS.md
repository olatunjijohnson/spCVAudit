# spCVAudit 0.1.0

Initial release of a focused spatial cross-validation audit package.

- Declare equally or unequally weighted prediction targets.
- Construct and validate random, spatial-block, buffered and optional CAST kNNDM folds.
- Calculate exact weighted empirical A, W, J1, J2, U1 and U2 contrasts.
- Compare candidate schemes with common covariance and explicit tie probabilities.
- Estimate exponential covariance using a constrained variogram profile search.
- Calculate Gaussian reference MSE, the additional-neighbour decomposition,
  its positive-nugget bound, and conditional target MSE.
- Supply four worked tutorials and a complete function reference.

This is a separate package from spLeakage. The paper's retained computations keep
their original package snapshot; the new interface does not change those results.
