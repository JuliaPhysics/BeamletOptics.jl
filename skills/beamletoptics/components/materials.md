# Refractive index / materials

A refractive index is anything of type
`RefractiveIndex = Union{Function, DiscreteRefractiveIndex, SellmeierEquation}`. It must be callable as
`n(λ)` with **λ in meters** and return a real number. Plain numbers are **not** accepted.

```julia
n_const = λ -> 1.5                                    # constant index

NBK7 = DiscreteRefractiveIndex([532e-9, 1064e-9],     # exact lookup table,
                               [1.5195, 1.5066])      # NO interpolation: any other λ -> KeyError

N_BK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945,     # B1, B2, B3
                          0.00600069867, 0.0200179144, 103.560653)  # C1, C2, C3 in µm²
N_BK7(532e-9)                                         # still called with meters
```

- Coefficients from refractiveindex.info or glass catalogs (µm-based) go straight into `SellmeierEquation`.
- The surrounding medium is always n = 1.
- Constructors validate the index by calling it with `1`, `1f0` and `1.0`, so closures must accept any `Real`.
- Materials are used by lenses, prisms, plate/cube beamsplitters, compensator plates and linear polarizers.
