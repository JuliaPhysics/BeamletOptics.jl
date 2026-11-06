# Refractive index / materials

A refractive index is a subtype of `BeamletOptics.AbstractRefractiveIndex`: `ConstantRefractiveIndex`,
`DiscreteRefractiveIndex` or `SellmeierEquation`. It is callable as `n(λ)` with **λ in meters** and returns
a real number. Constructors take only these types: a plain number (`1.5`) or a function (`λ -> 1.5`)
throws a `MethodError`.

```julia
n_const = ConstantRefractiveIndex(1.5)                # constant index, no dispersion

NBK7 = DiscreteRefractiveIndex([532e-9, 1064e-9],     # exact lookup table,
                               [1.5195, 1.5066])      # NO interpolation: any other λ -> KeyError

N_BK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945,     # B1, B2, B3
                          0.00600069867, 0.0200179144, 103.560653)  # C1, C2, C3 in µm²
N_BK7(532e-9)                                         # still called with meters
```

- Coefficients from refractiveindex.info or glass catalogs (µm-based) go straight into `SellmeierEquation`.
- The surrounding medium is always n = 1.
- A constant index ignores dispersion. Use `SellmeierEquation` or `DiscreteRefractiveIndex` whenever the
  result depends on λ (several wavelengths, chromatic focal shift, prisms for dispersion).
- Another dispersion formula is a new subtype with a functor:
  `struct Cauchy <: BeamletOptics.AbstractRefractiveIndex; A::Float64; B::Float64 end` and
  `(c::Cauchy)(λ) = c.A + c.B / (λ * 1e6)^2`. It must accept any `Real` λ and return a `Real`.
- An absorbing material, e.g. the metal layer of a `ThinFilmCoating`, is
  `AbsorbingRefractiveIndex(n, κ)` with the extinction coefficient κ (a number or an index type of κ(λ));
  `extinction_coefficient(n, λ)` returns it, 0 for the other types.
- Materials are used by lenses, prisms, plate/cube beamsplitters, compensator plates and linear polarizers.
