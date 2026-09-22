```@setup rays
dir = joinpath(@__DIR__, "..", "assets")

Main.DocUtils.conditional_include(joinpath(dir, "api_assets", "raymarching.jl"))
Main.DocUtils.conditional_include(joinpath(dir, "api_assets", "sdf_showcase.jl"))
```

# Signed Distance Functions (SDFs)

Signed distance functions in general represent geometry by providing a function that e.g. maps a point in R³ onto the closest scalar distance to the surface of the geometry [Osher:2003](@cite). This representation technique has already been used extensively for ray tracing applications in the field of computer graphics. For a quick introduction into SDFs the [website of Inigo Quilez](https://iquilezles.org/articles/distfunctions/) is referred to.

In the field of ray-based optical simulations, mathematical surface definitions, e.g. in the form of analytical equations for the surface sag or NURBs, are common [Sasian:2019; p. 44 ff](@cite). However, they are also costly to implement since no "one-size-fits-all" algorithm exists which could be implemented for the `intersect3d` function. This is why for this package, SDFs have been chosen to represent curved surfaces. This is due to the fact that for exact SDFs

1. the point of intersection can be calculated exactly with reasonable convergence criteria
2. the normal vector at the point of intersection can be calculated exactly using automatic differentiation
3. the ray marching algorithm works for any shape that provides a SDF
4. the ray marching procedure can be computed very efficiently

The viability of SDFs for radiation transfer simulations has already been demonstrated in the literature [McMillan:2022](@cite). **For optical ray tracing this is a new application to the best of our knowledge**. Representing lens shapes using the exact shape SDFs provided by Inigo Quilez is straight forward by combining cylinder and cut-sphere SDFs using boolean addition and subtraction. Further, custom iterative algorithms have been developed by O. Kliebisch in order to tie in (rotationally symmetrical) aspherical lenses into the ray marching formalism.

In order to implement new SDFs, the `AbstractSDF` interface should be used.

```@docs; canonical=false
BeamletOptics.AbstractSDF
```

The following shapes are currently implemented:

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractSDF);
```

## Ray marching algorithm

The algorithm works by iteratively propagating a ray from the point of origin through a scene. At each step, the SDF/s is/are evaluated at the current point along the ray. The scaral value tells the algorithm how far the current positions is from the nearest surface. This distance is guaranteed to be a safe step size: the ray can advance by that amount without “tunneling” through any geometry. The process repeats until one of two conditions is met: (1) the distance becomes smaller than a threshold (the ray has hit a surface), or (2) a maximum number of steps or maximum distance is exceeded (the ray escapes without hitting anything). A visualization of this process is provided below. A ray propagates from the left to the right end of the scene. The step size is visualized as the radius of a sphere for each iteration.

![Ray marching algorithm showcase](raymarching.png)

Surface normals can be estimated numerically by sampling the SDF gradient around the hit point or by generating the automatic derivative of the SDF for exact descriptions. BMO uses both procedures, depending on the SDF type.

## Composite SDFs

In order to make use of the dispatch capability of Julia, boolean composite SDFs allow users to easily combine SDFs using the `+` and `-` operators. Both share a common field layout and pivot-aware kinematics, defined by the `AbstractCompositeSDF` interface:

```@docs; canonical=false
BeamletOptics.AbstractCompositeSDF
```

![SDF composition examples](composite_sdf.png)

### Union

```@docs; canonical=false
BeamletOptics.UnionSDF
```

### Difference

```@docs; canonical=false
BeamletOptics.DifferenceSDF
```

!!! info "Exactness"
    Boolean addition of non-overlapping, exact SDFs stays exact. Boolean subtraction, by
    contrast, only yields a *bound*: `max(-a, b)` under-estimates the true distance near
    the seam between operands. This under-estimate is the *safe* direction for sphere
    tracing — it never causes tunneling, only a few extra ray marching iterations near
    concave creases. For information refer to the website of [I. Quilez](https://iquilezles.org/articles/distfunctions/)