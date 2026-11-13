# Validated swordfish diagnostic snapshot

Generated from the canonical Quadra fit and the eight source files recorded in
`../source_manifest.json`. All source hashes were rechecked unchanged.

See `../../FIT_DIAGNOSTICS.md` for methods, limitations, and commands that
regenerate these outputs in `build/`. This snapshot preserves reviewable numerical
results without tracking executable binaries or the full covariance matrices.

Validation: three perturbed starts converge to the same objective; two Hessian
steps give stable positive definite curvature; 39 selected standard errors and
ten management reference points agree with SS3 at their documented tolerances.
