# Style guide

This outlines recommended practices for contributors to this package.

## Naming

- Type parameters: use `F` for `F <: Fixed`, `N` for `N <: Normed`,
  and `X` for `X <: FixedPoint`. Use `f` for the number of fractional bits.
- Use `Ti` for `Ti <: Integer`, `Tf` for `Tf <: AbstractFloat`, and `Tw`
  for `widen`ed types.
- `T` should refer to the underlying "raw" type.
