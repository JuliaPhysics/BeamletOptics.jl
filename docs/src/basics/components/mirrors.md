```@setup mirrors
mirror_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "mirror_renders")

Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "plano_mirror_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "spherical_mirror_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "parabolic_mirror_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "oap_mirror_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "ellipsoidal_mirror_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "hyperbolic_mirror_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "conic_frame_convention.jl"))
```

# Mirrors

A common optical element with a straight-forward optical interaction. This kind of component is in general defined as a [`BeamletOptics.AbstractReflectiveOptic`](@ref). For a basic [`Ray`](@ref) the interaction is simply defined by the [`BeamletOptics.reflection3d`](@ref) function. A more complex algorithm is required when when a [`PolarizedRay`](@ref) interacts with a reflecting surface. The polarization calculus that is performed is explained in the [Polarized rays](@ref) section. Below, some of the concrete implemented mirror types are shown. In general, the [`Mirror`](@ref) is used as a concrete type to represent an arbitrary reflecting shape.

```@docs; canonical=false
Mirror
```

The following constructors can be used to generate flat reflecting shapes. Additional types are explained below.

- [`SquarePlanoMirror2D`](@ref)
- [`SquarePlanoMirror`](@ref)
- [`RectangularPlanoMirror`](@ref)
- [`Retroreflector`](@ref)


## Plano Mirrors

A category of mirrors with a flat reflecting surface. A round version of this mirror can be easily generated using the [`RoundPlanoMirror`](@ref) or [`RightAnglePrismMirror`](@ref) types. An optional central through-hole can be added to [`RoundPlanoMirror`](@ref) via the `hole_diameter` keyword argument:

```@docs; canonical=false
RoundPlanoMirror(::Real, ::Real)
```

Below, a trivial example of a beam path propagating through a system of Ø1"-mirrors mounted in [KM100CP/M](https://www.thorlabs.de/thorproduct.cfm?partnumber=KM100CP/M#ad-image-0) kinematic mounts is shown (e.g. [PF10-03-P01](https://www.thorlabs.com/thorproduct.cfm?partnumber=PF10-03-P01)). Note that the mounts are modeled as [`NonInteractableObject`](@ref)s.

![Plano mirror showcase](plano_mirror_showcase.png)

## Spherical Mirrors

The [`SphericalMirror`](@ref) represents an ideal optical element with a spherical concave reflective surface, commonly used for non-dispersive focusing applications. Its geometry is modeled using a combination of a concave spherical surface and a plano substrate, represented internally by a [`BeamletOptics.UnionSDF`](@ref) (refer also to the [SDF-based spherical lenses](@ref) section). An optional central through-hole can be added via the `hole_diameter` keyword argument.

![Spherical mirror multipass showcase](spherical_mirror_showcase.png)

The following constructor allows the spawning of spherical mirrors.

```@docs; canonical=false
SphericalMirror(::Real, ::Real, ::Real)
```

## Conic Mirrors

All rotationally-symmetric conic mirrors (paraboloids, ellipsoids and hyperboloids) share the same underlying surface of revolution, i.e. the [`BeamletOptics.ConicSDF`](@ref) which is defined by the governing equation

$$Z(r) = \frac{r^2}{R\left(1 + \sqrt{1 - (1+k)\,r^2/R^2}\right)},$$

where the surface is concave for $R > 0$ (i.e. opens towards $-y$) and convex for $R < 0$. This equation is the ISO 10110-12 aspheric surface description without the power series [ISO10110-12:2019, Sasian:2019](@cite) and also used for the description of aspherical lenses within this package (albeit in a different form, see [Aspherical lenses](@ref)). The $k$-factor determines the surface shape. Depending on its value, the surface type will be either one of the ones listed in the table below.

```@raw html
<div class="bmo-table">
```

| $k$ | surface family |
| :---: | --- |
| $k < -1$ | hyperboloid |
| $k = -1$ | paraboloid |
| $-1 < k < 0$ | prolate ellipsoid |
| $k = 0$ | sphere |
| $k > 0$ | oblate ellipsoid |

```@raw html
</div>
```

For $k > -1$ the parent conic is only defined for $r < \lvert R \rvert / \sqrt{1+k}$; for $k \le -1$ every aperture is admissible.

![Conic mirror frame convention](conic_frame_convention.png)

Ellipsoidal and hyperbolic mirrors are parameterized by the object/image conjugate distances $s, s'$ (measured from the parent vertex, positive towards $-y$) rather than $R, k$ directly:

$$R = \frac{2ss'}{s+s'}, \qquad k = -\left(\frac{s'-s}{s'+s}\right)^2$$

Same-sign $s, s'$ give a real second focus (ellipsoid, $-1 < k \le 0$); opposite signs give a virtual second focus (hyperboloid, $k < -1$).

The following constructors allow the spawning of on-axis and off-axis conic, ellipsoidal and hyperbolic mirrors. All on-axis constructors accept an optional `hole_diameter` keyword argument to subtract an axial cylindrical through-hole from the substrate (e.g. for Cassegrain, Gregorian, Ritchey-Chrétien, or Dall-Kirkham telescope primaries). Parabolic mirrors are covered in the [Parabolic Mirrors](@ref) section below.

```@docs; canonical=false
ConicMirror(::Real, ::Real, ::Real)
OffAxisConicMirror(::Real, ::Real, ::Real, ::Real)
```

### Parabolic Mirrors

Parabolic mirrors are the $k = -1$ special case of the conic mirrors above, with $R = 2f$. In contrast to the [`SphericalMirror`](@ref), a paraboloid focuses a collimated beam that is parallel to its optical axis into a single point without spherical aberration. BMO provides an on-axis and an off-axis variant.

#### On-axis parabolic mirrors

The [`ParabolicMirror`](@ref) represents an on-axis parabolic mirror. Its surface $y = -\frac{x^2 + z^2}{4f}$ opens towards $-y$, such that the focus lies at $(0, -f, 0)$. Optionally, a central through-hole can be added via `hole_diameter`, e.g. to model a Cassegrain primary.

![Parabolic mirror showcase](parabolic_mirror_showcase.png)

```@docs; canonical=false
ParabolicMirror(::Real, ::Real)
```

#### Off-axis parabolic mirrors

The [`OffAxisParabolicMirror`](@ref) represents an off-axis parabolic (OAP) mirror used for achromatic focusing and beam deflection without introducing spherical aberration.

Its geometry is constructed from a parent paraboloid with focal length $f$ and off-axis distance $x_{\text{off}}$, parameterized by the Reflected Focal Length ($RFL$) and deflection angle $\theta_d$ (default 90°):

$$f = RFL \cdot \cos^2\left(\frac{\theta_d}{2}\right), \quad x_{\text{off}} = RFL \cdot \sin(\theta_d)$$

##### Through-holes (Thorlabs POH style)

An optional through-hole can be specified via the `hole_diameter` keyword argument. The orientation of the bore is controlled by `hole_axis`:
- `:collimated` (default): A cylindrical bore parallel to the incident collimated beam (local $y$-axis / substrate normal).
- `:focused`: A cylindrical bore oriented towards the focal point $(-x_{\text{off}}, -RFL\cos\theta_d, 0)$, enabling collinear pump-probe or THz transmission through the mirror substrate directly onto the focus.

![Off-Axis Parabolic mirror showcase](oap_mirror_showcase.png)

```@docs; canonical=false
OffAxisParabolicMirror(::Real, ::Real)
```

### Ellipsoidal mirrors

Ellipsoidal mirrors cover the $-1 < k \le 0$ range of the conic mirrors above, i.e. prolate ellipsoids with the sphere ($s = s'$) as limiting case. Both conjugate foci are real and lie in front of the mirror, such that a point source placed in one focus is imaged into the other one without spherical aberration. Typical applications are refocusing/relay mirrors, e.g. to couple a source into a fiber, or the concave secondary of a Gregorian telescope. The [`EllipsoidalMirror`](@ref) places both foci on its optical axis, while the [`OffAxisEllipsoidalMirror`](@ref) uses an off-axis segment of the parent ellipsoid in order to separate the reflected from the incoming beam.

Below, a fan of rays is emitted from the far focus of an [`EllipsoidalMirror`](@ref) (left). After the reflection, all rays are focused into the near focus of the mirror.

![Ellipsoidal mirror showcase](ellipsoidal_mirror_showcase.png)

```@docs; canonical=false
EllipsoidalMirror(::Real, ::Real, ::Real)
OffAxisEllipsoidalMirror(::Real, ::Real, ::Real, ::Real)
```

### Hyperbolic mirrors

Hyperbolic mirrors cover the $k < -1$ range of the conic mirrors above. Exactly one of the two conjugate foci is real, the other one is virtual and lies behind the mirror, i.e. $s$ and $s'$ have opposite signs. A beam converging towards the virtual focus is therefore reflected into the real focus without spherical aberration. The classic application is the convex secondary of a Cassegrain telescope, which reimages the prime focus of a parabolic primary (see [Parabolic Mirrors](@ref)). Since this prime focus lies behind the secondary, it is passed as a negative $s$. As before, the [`HyperbolicMirror`](@ref) is the on-axis variant and the [`OffAxisHyperbolicMirror`](@ref) the off-axis segment of the parent hyperboloid.

Below, a classical Cassegrain telescope is shown. A collimated beam (left) is focused by a [`ParabolicMirror`](@ref) with a central bore towards its prime focus. Before reaching it, the rays are intercepted by a convex [`HyperbolicMirror`](@ref), which reimages the prime focus through the bore into the final focus behind the primary (black).

![Hyperbolic mirror showcase](hyperbolic_mirror_showcase.png)

```@docs; canonical=false
HyperbolicMirror(::Real, ::Real, ::Real)
OffAxisHyperbolicMirror(::Real, ::Real, ::Real, ::Real)
```