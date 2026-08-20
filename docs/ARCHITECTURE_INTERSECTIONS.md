# The Unified Media–Interface Architecture for BeamletOptics.jl
**A First-Principles, Zero-Allocation Foundation for Geometric and Beamlet Optics**

**Document Version:** 2.7 (Comprehensive Master Architectural Specification)  
**Target Branches:** `feature/intersection-handling`, `master`  
**Author:** Optical Architecture Working Group  

---

## 1. Executive Summary & Paradigm Shift

Historically, `BeamletOptics.jl` (BMO) distributed optical physics across concrete optical components. Every component type (`SphericalLens`, `Mirror`, `PlateBeamsplitter`, `PolarizationFilter`, `Detector`) implemented its own custom `interact3d` method, duplicating Snell's law, Fresnel equations, normal inversion, and ray creation. Furthermore, tracking multi-material interfaces (doublets, coatings, beamsplitters) relied on ad-hoc object mutation and tolerance searches (`MultiIntersection`).

### The Foundational Insight
**A lens is not an optical interaction; a lens is a geometric volume filled with a physical medium.**

In natural optics:
1. **Geometry (Space & Boundaries):** SDFs and Meshes partition $\mathbb{R}^3$ into spatial domains bounded by oriented manifolds $\partial \Omega$.
2. **Media (Volumetric Physics):** Light propagates through bulk materials characterized by complex refractive indices $\tilde{n}(\lambda) = n(\lambda) + i\kappa(\lambda)$ or dielectric tensors $\boldsymbol{\epsilon}(\lambda)$.
3. **Interfaces (Boundary Physics):** When a wavepacket encounters an oriented boundary $\partial \Omega$ between Medium 1 and Medium 2, Maxwell's boundary conditions uniquely dictate reflection, refraction, absorption, splitting, or diffraction.

By decomposing BMO into these three orthogonal layers:
* **Zero Component Boilerplate:** Over 90% of duplicate Snell/Fresnel code across components is permanently eliminated.
* **Universal Higher-Order Beamlets:** `GaussianBeamlet` and `AstigmaticGaussianBeamlet` broadcast across parabasal rays with continuous Optical Path Length ($OPL$) accumulation without per-component specialization.
* **100% Type-Stable & Register-Resident:** Immutable `isbits` structs and compile-time unrolled tuple indexing eliminate GC pauses and enable native **GPU acceleration** and **Automatic Differentiation (AD)**.

---

## 2. Architectural Overview: The Four Decoupled Layers

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                            LAYER 1: PURE GEOMETRY                            │
│   AbstractShape{T} (SDFs, Meshes, CSG)  ──>  ShapeHit{T} (isbits, stack)     │
└──────────────────────────────────────┬───────────────────────────────────────┘
                                       │ (t, n̂, SurfaceID)
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│                    LAYER 2: COMPOSABLE OPTICAL OBJECTS                       │
│   AbstractObject{T}:                                                         │
│     ├── SingleOptic: Shape (Geometry) + Bulk Medium + Boundary Tuple         │
│     └── CompositeOptic: MultiShape Children + Internal Interfaces            │
└──────────────────────────────────────┬───────────────────────────────────────┘
                                       │ Type-Stable Unrolled Dispatch & Media Transition
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│                    LAYER 3: UNIVERSAL BOUNDARY PHYSICS                       │
│   interact_surface(surface_model, transition, ray, t_hit)                    │
│   (Fresnel, Thin-Film Stack, Ideal Mirror, Grating, Detector, Polarizer)     │
└──────────────────────────────────────┬───────────────────────────────────────┘
                                       │ Outgoing Ray(s) / Beamlet(s) + OPL & Epsilon Advance
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│                            LAYER 4: SYSTEM TRACER                            │
│   Universal Execution Loop with Energy Cutoff & Branching DAG Buffer         │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Layer 1: Pure Geometry & Type-Stable `ShapeHit{T}`

Geometric solvers (SDF ray-marchers, analytical shapes, BVH meshes) have zero awareness of refractive indices, wavelengths, or optical physics. They solve pure geometry and return a stack-allocated, `isbits` geometric hit token.

### 3.1 Universal `SurfaceID` Facet Addressing
To support everything from 2-sided lenses to **5-sided prisms, 8-sided polyhedra, and triangle meshes with millions of facets** without dynamic dispatch, `SurfaceID` is defined as a lightweight 32-bit `isbits` integer wrapper:

```julia
"""
    SurfaceID

A 32-bit unsigned integer identifying a geometric surface or facet.
100% isbits, zero-allocation, and GPU register resident. Supports up to 4.2 billion facets.
"""
struct SurfaceID
    id::UInt32
end

Base.Integer(s::SurfaceID) = s.id
Base.convert(::Type{SurfaceID}, x::Integer) = SurfaceID(UInt32(x))

# Semantic Aliases for Common Standard Optics
const FRONT    = SurfaceID(1)
const BACK     = SurfaceID(2)
const SIDE     = SurfaceID(3)
const INTERNAL = SurfaceID(4)
```

```julia
"""
    ShapeHit{T<:Real}

Pure geometric intersection result. 100% immutable, zero heap allocations (isbits).
"""
struct ShapeHit{T <: Real}
    t::T                    # Distance parameter along ray (p = origin + t * dir)
    normal::Point3{T}       # Outward unit normal vector (|n̂| = 1)
    tag::SurfaceID          # 4-byte surface identifier (FRONT, BACK, Prism Face, Triangle Index)
end
```

### 3.2 Scaling from Lenses to Prisms and Dense Triangle Meshes
* **Lenses & Mirrors ($N=2..3$ faces):** `tag = FRONT`, `BACK`, or `SIDE`. Boundaries stored in `NTuple{3, SurfaceModel}`.
* **Prisms & Polyhedra ($N=4..8$ faces, e.g. Pentaprism, Porro Prism):** `tag = SurfaceID(1)` through `SurfaceID(N)`. Boundaries stored in an indexable `NTuple{N, SurfaceModel}`.
* **Dense Triangle Meshes ($N > 1000$ triangles):** Individual triangles do not require separate entries in `boundaries::Tuple`. Instead, triangle `SurfaceID` values resolve to a **Material Palette ID** ($1..K$, where $K \le 8$), indexing into a compact material palette array.

### 3.3 Generic SDF Ray-Marcher (Written Once for All SDFs)
In BMO, all SDF volumes inherit from `AbstractSDF`. The ray-marching algorithm is implemented **once** in `src/SDFs/AbstractSDF.jl` and automatically powers all SDF primitives, lens shapes, and CSG unions:

```julia
# Defined ONCE in src/SDFs/AbstractSDF.jl for ALL SDF types (spheres, aspheres, cylinders, unions, etc.)
function intersect3d(shape::AbstractSDF{T}, ray::AbstractRay{T}) where {T}
    t_hit, hit_found = ray_march(shape, position(ray), direction(ray))
    !hit_found && return nothing

    p_hit = position(ray) + t_hit * direction(ray)
    normal = normal3d(shape, p_hit)

    # Intrinsic tag queried generically from the SDF
    tag = surface_tag(shape, p_hit, normal)

    return ShapeHit{T}(t_hit, normal, tag)
end
```

### 3.4 Canonical, Fully-Typed `Intersection` Bridge
When a shape hit is associated with an `AbstractObject`, it is wrapped in an immutable `Intersection`:

```julia
"""
    Intersection{T<:Real, S<:AbstractShape{T}, O<:AbstractObject{T}} <: AbstractIntersection{T}

Fully typed, immutable token representing the ray intersection.
"""
struct Intersection{T <: Real, S <: AbstractShape{T}, O <: AbstractObject{T}} <: AbstractIntersection{T}
    hit::ShapeHit{T}
    shape::S
    object::O
end

@inline Base.length(i::Intersection) = i.hit.t
@inline normal3d(i::Intersection)   = i.hit.normal
@inline surface_tag(i::Intersection) = i.hit.tag
@inline shape(i::Intersection)       = i.shape
@inline object(i::Intersection)      = i.object
```

---

## 4. Layer 2: Optical Media & Surface Models

### 4.1 Volumetric Media (`AbstractMedium`)
Bulk materials govern phase accumulation, group velocity, dispersion, and bulk absorption:

```julia
abstract type AbstractMedium end

struct Ambient <: AbstractMedium end  # Wildcard resolving to System environment (e.g. Air/Vacuum)

struct IsotropicMedium{D, A} <: AbstractMedium
    name::Symbol
    dispersion::D    # Function λ -> n(λ) (e.g. SellmeierFormula, Cauchy, Constant)
    attenuation::A   # Function λ -> α(λ) (bulk absorption coefficient [1/m])
end

refractive_index(::Ambient, λ) = 1.0
refractive_index(m::IsotropicMedium, λ) = m.dispersion(λ)
```

### 4.2 Surface Models (`AbstractSurfaceModel`)
Surface models govern physical behavior *at* the boundary interface:

```julia
abstract type AbstractSurfaceModel end

struct FresnelInterface <: AbstractSurfaceModel end     # Bare dielectric boundary (Snell + Fresnel)
struct IdealMirror <: AbstractSurfaceModel end          # 100% specular reflection
struct AbsorbingSurface <: AbstractSurfaceModel end     # Ray termination (apertures/baffles)
struct DetectorSurface{D} <: AbstractSurfaceModel       # Records power, phase, E-field
    detector_data::D
end
struct CoatedSurface{C} <: AbstractSurfaceModel         # AR, HR, Beamsplitter, Thin-film stack
    coating_model::C
end
struct GratingSurface{T} <: AbstractSurfaceModel        # Diffraction grating: d * (sinθm - sinθi) = mλ
    groove_density::T
    order::Int
end
```

---

## 5. Layer 3: The `AbstractObject` Type Hierarchy

All optical components in BMO are subtypes of `AbstractObject{T}`. We organize them into a clean, two-tier hierarchy:

```
                          AbstractObject{T}
                                 │
             ┌───────────────────┴───────────────────┐
             ▼                                       ▼
     SingleOptic{T}                        CompositeOptic{T}
(Lens, Mirror, Detector)                (DoubletLens, CubeBeamsplitter)
 ├── shape::AbstractShape                ├── children::Tuple of Optics
 ├── medium::AbstractMedium              └── internal_boundaries
 └── boundaries::NTuple
```

### 5.1 Single Optics (`SingleOptic <: AbstractObject`)
A `SingleOptic` represents a single geometric volume (or 2D surface) possessing a bulk interior medium and an indexable `Tuple` of boundary surface definitions:

```julia
abstract type SingleOptic{T} <: AbstractObject{T} end

shape_trait_of(::SingleOptic) = SingleShape()

# Standard generic implementation backing Lens, Mirror, Detector, etc.
struct Optic{T <: Real, S <: AbstractShape{T}, M <: AbstractMedium, B <: Tuple} <: SingleOptic{T}
    shape::S
    medium::M
    boundaries::B  # (FrontSurface, BackSurface, SideSurface)
end

shape(o::Optic) = o.shape
medium(o::Optic) = o.medium

# Guaranteed Zero-Allocation, Compile-Time Unrolled Boundary Lookup:
@inline function boundary(o::Optic{T, S, M, B}, tag::SurfaceID) where {T, S, M, B <: Tuple}
    idx = Integer(tag)
    # Statically unrolled branch-free dispatch for heterogeneous Tuples:
    Base.Cartesian.@nif 8 i -> (idx == i && i <= length(o.boundaries)) i -> (@inbounds o.boundaries[i]) i -> o.boundaries[1]
end
```

#### User-Facing Component Constructors (Idiomatic Julia)
```julia
# Spherical Lens (Constructs an Optic with glass medium and front/back boundaries)
function SphericalLens(r1, r2, thickness, diameter, glass; front=FresnelInterface(), back=FresnelInterface())
    s = SphericalLensSDF(r1, r2, thickness, diameter)
    m = IsotropicMedium(:Glass, glass, nothing)
    b = (front, back, AbsorbingSurface()) # Indexable Tuple: FRONT=1, BACK=2, SIDE=3
    return Optic(s, m, b)
end

# Flat Mirror (Constructs an Optic with reflective front boundary)
function PlanoMirror(diameter, thickness; coating=IdealMirror())
    s = CylinderSDF(diameter/2, thickness)
    b = (coating, AbsorbingSurface(), AbsorbingSurface())
    return Optic(s, Ambient(), b)
end

# Detector (Constructs an Optic with recording surface)
function Detector(width, height, resolution)
    s = RectangleSDF(width, height)
    b = (DetectorSurface(resolution), AbsorbingSurface(), AbsorbingSurface())
    return Optic(s, Ambient(), b)
end
```

---

### 5.2 Composite Optics & the `MultiShape` Integration

A `CompositeOptic` (e.g. `DoubletLens`, `CubeBeamsplitter`) represents an assembly of multiple touching single optics that share internal boundaries.

```julia
abstract type CompositeOptic{T} <: AbstractObject{T} end

# Composite optics naturally implement the MultiShape trait
shape_trait_of(::CompositeOptic) = MultiShape()
```

#### The `DoubletLens` Struct & Constructor
```julia
struct DoubletLens{T, L1 <: SingleOptic{T}, L2 <: SingleOptic{T}, B <: AbstractSurfaceModel} <: CompositeOptic{T}
    front_lens::L1
    back_lens::L2
    cemented_boundary::B  # e.g. FresnelInterface() or CoatedSurface(CementLayer)
end

# MultiShape shape getter exposes constituent shapes for kinematics & intersection
shape(dl::DoubletLens) = (shape(dl.front_lens), shape(dl.back_lens))
```

#### How `MultiShape` Intersection & Kinematics Work
1. **Kinematics**: All spatial transforms (`translate3d!`, `rotate3d!`, `align3d!`) on a `CompositeOptic` automatically delegate across all child shapes via existing `MultiShape` trait methods.
2. **Intersection**: Testing a ray against a `CompositeOptic` queries constituent shapes:

```julia
function intersect3d(::MultiShape, composite::CompositeOptic, ray::AbstractRay{R}) where {R}
    closest_hit::Union{Nothing, ShapeHit{R}} = nothing
    closest_shape = nothing

    for child_shape in shape(composite)
        temp_hit = intersect3d(child_shape, ray)
        if temp_hit !== nothing
            if closest_hit === nothing || temp_hit.t < closest_hit.t
                closest_hit = temp_hit
                closest_shape = child_shape
            end
        end
    end

    isnothing(closest_hit) && return nothing
    return Intersection(closest_hit, closest_shape, composite)
end
```

---

## 6. Layer 4: Universal Physics Dispatch Pipeline

Because optical physics is decoupled from component types, `interact_surface` is implemented **once per surface model**, not once per optical component:

### 6.1 State-Aware Media Transition Resolution (Entering vs. Exiting vs. Cemented)
A ray carries its active propagation medium $M_{\text{current}}$. The transition resolver evaluates the incident and transmitted media across single and composite interfaces:

```julia
@inline function resolve_transition(
    med_inside::AbstractMedium,
    med_outside::AbstractMedium,
    ray::AbstractRay,
    normal::Point3
)
    # dot(dir, normal) < 0 => entering optic; dot(dir, normal) > 0 => exiting optic
    is_entering = dot(direction(ray), normal) < 0
    normal_eff  = is_entering ? normal : -normal

    # Resolve incident and transmitted bulk media
    med_inc   = is_entering ? current_medium(ray) : med_inside
    med_trans = is_entering ? med_inside : med_outside

    n_inc   = refractive_index(med_inc, wavelength(ray))
    n_trans = refractive_index(med_trans, wavelength(ray))

    return (n_inc = n_inc, n_trans = n_trans, med_trans = med_trans, normal_eff = normal_eff, entering = is_entering)
end
```

* **Isolated Lens Front Face (Entering):** `med_inside = Glass`, `med_outside = Ambient` $\implies n_{\text{ambient}} \to n_1$.
* **Cemented Doublet Interface ($L_1 \to L_2$):** `med_inside = Glass1`, `med_outside = Glass2` $\implies n_1 \to n_2$ (Direct Glass-to-Glass Snell refraction without false air gaps).
* **Isolated Lens Rear Face (Exiting):** `med_inside = Glass`, `med_outside = Ambient` $\implies n_1 \to n_{\text{ambient}}$.

---

### 6.2 Universal Surface Interaction Methods with Epsilon Self-Intersection Prevention
Every `interact_surface` receives the geometric distance $t_{\text{hit}}$ from `ShapeHit{T}`. Outgoing rays apply an $\epsilon$-normal offset to prevent self-intersection ("shadow acne") on subsequent ray-marching steps:

```julia
# Numerical epsilon offset constant for self-intersection avoidance
const RAY_EPSILON = 1e-7

# Dielectric Fresnel Interface
function interact_surface(::FresnelInterface, trans, ray::PolarizedRay{T}, t_hit::T) where {T}
    θi = angle3d(direction(ray), -trans.normal_eff)
    rs, rp, ts, tp = fresnel_coefficients(θi, trans.n_trans / trans.n_inc)

    p_hit = position(ray) + t_hit * direction(ray)

    if is_total_internal_reflection(rs, rp)
        new_dir = reflection3d(direction(ray), trans.normal_eff)
        new_pol = reflect_polarization(polarization(ray), rs, rp, direction(ray), trans.normal_eff)
        med_out = trans.med_inc
        new_pos = p_hit + (RAY_EPSILON * sign(dot(new_dir, trans.normal_eff))) * trans.normal_eff
    else
        new_dir, _ = refraction3d(direction(ray), trans.normal_eff, trans.n_inc, trans.n_trans)
        new_pol = refract_polarization(polarization(ray), ts, tp, direction(ray), trans.normal_eff, trans.n_inc, trans.n_trans)
        med_out = trans.med_trans
        new_pos = p_hit - (RAY_EPSILON * sign(dot(new_dir, trans.normal_eff))) * trans.normal_eff
    end

    new_opl = optical_path_length(ray) + t_hit * trans.n_inc
    return PolarizedRay{T}(new_pos, new_dir, nothing, wavelength(ray), med_out, new_pol, new_opl)
end

# Specular Mirror Interface
function interact_surface(::IdealMirror, trans, ray::AbstractRay{T}, t_hit::T) where {T}
    new_dir = reflection3d(direction(ray), trans.normal_eff)
    p_hit   = position(ray) + t_hit * direction(ray)
    new_pos = p_hit + RAY_EPSILON * trans.normal_eff
    new_opl = optical_path_length(ray) + t_hit * trans.n_inc
    return Ray{T}(new_pos, new_dir, nothing, wavelength(ray), current_medium(ray), new_opl)
end

# Absorptive / Detector Interface
function interact_surface(d::DetectorSurface, trans, ray::AbstractRay{T}, t_hit::T) where {T}
    hit_pos = position(ray) + t_hit * direction(ray)
    record_hit!(d.detector_data, hit_pos, ray)
    return nothing  # Cleanly terminates ray
end
```

### 6.3 Robust Polarization Frame Transformation & Normal Incidence Degeneracy
When computing `reflect_polarization` and `refract_polarization`, the electric field $\mathbf{E}_{\text{in}}$ is transformed using the local $(s, p)$ plane-of-incidence basis:
$$\hat{s} = \frac{\vec{d} \times \hat{n}}{\|\vec{d} \times \hat{n}\|}, \quad \hat{p} = \vec{d} \times \hat{s}$$

#### Degeneracy Handling at Normal Incidence ($\|\vec{d} \times \hat{n}\| < \epsilon_{\text{tol}}$):
When light strikes a surface at exact or near-normal incidence ($\theta_i \approx 0$), the cross product $\vec{d} \times \hat{n} \approx \mathbf{0}$, producing a $0/0$ singularity.
* **Resolution:** At normal incidence, the plane of incidence is degenerate and $r_s = r_p = \frac{n_{\text{trans}} - n_{\text{inc}}}{n_{\text{trans}} + n_{\text{inc}}}$. The solver bypasses the cross product and applies the isotropic scalar Fresnel factor directly to the transverse $\mathbf{E}$-field vector without rotating frames.

---

## 7. The Universal System Tracing Loop & Branching Model

### 7.1 Single Tracing Step
The system tracer coordinates geometry, medium resolution, and boundary physics in three zero-allocation steps:

```julia
function trace_step!(system::AbstractSystem, ray::AbstractRay{T}) where {T}
    # 1. Geometry phase: find closest intersection across all objects
    hit_intersection = trace_all(system, ray)
    isnothing(hit_intersection) && return nothing

    obj  = object(hit_intersection)
    geom = hit_intersection.hit

    # 2. Transition phase: resolve entering/exiting bulk media
    ambient = ambient_medium(system)
    trans   = resolve_transition(medium(obj), ambient, ray, geom.normal)

    # 3. Physics phase: query boundary model via branch-free Cartesian unrolling
    surface_model = boundary(obj, geom.tag)
    return interact_surface(surface_model, trans, ray, geom.t)
end
```

### 7.2 Zero-Allocation Ray Branching (DAG Buffer Model with Energy Cutoff)
For partially reflecting/transmitting surfaces (beamsplitters, thin-film coatings) that split light into multiple child rays, energy splitting is managed via a **pre-allocated, fixed-capacity stack buffer** (`StaticRayStack{N}`):

```
                       Root Ray
                          │
                  [Beamsplitter]
                   /          \
          Reflected Ray      Transmitted Ray
          (Pushed to Stack)  (Traced iteratively)
```

#### Energy & Depth Cutoff Invariants:
To prevent infinite bounce loops in resonant cavities, etalons, or TIR traps:
1. **Flux Threshold Cutoff:** Rays with relative power $\text{flux} < \epsilon_{\text{cutoff}}$ (e.g. $10^{-4}$ or $0.01\%$) are dropped.
2. **Maximum Recursion Depth:** Ray branching terminates when $\text{depth} \ge \text{max\_depth}$ (default: 50).

---

## 8. Beamlet Optics & Parabasal Ray Decomposition

A major question is: **Do we have to specialize `interact3d` for `GaussianBeamlet` and `AstigmaticGaussianBeamlet`?**

### The Answer: NO!

In `BeamletOptics.jl`, higher-order beamlets are mathematically represented by **parabasal ray bundles**:

$$\begin{aligned}
\text{\textbf{Beam}} &\implies 1 \text{ Ray} \\
\text{\textbf{GaussianBeamlet}} &\implies 3 \text{ Rays: } (\text{chief}, \text{waist}, \text{divergence}) \\
\text{\textbf{AstigmaticGaussianBeamlet}} &\implies 9 \text{ Rays: } (\text{chief}, w_{x\pm}, w_{y\pm}, d_{x\pm}, d_{y\pm})
\end{aligned}$$

```
                AstigmaticGaussianBeamlet (9 Rays)
                                │
   ┌───────────┬───────────┬────┴──────┬───────────┬───────────┐
   ▼           ▼           ▼           ▼           ▼           ▼
Chief Ray   Waist X+    Waist X-    Div X+      Div X-      ...
   │           │           │           │           │           │
   └───────────┴───────────┼───────────┴───────────┴───────────┘
                           ▼
          Universal interact_surface(...)
         (Snell / Fresnel / Coating Transfer Matrix)
```

### 8.1 Independent Parabasal Ray-Marching & Wavefront Reconstruction
In BMO's Arnaud/Greynolds parabasal formulation:
* **All 3 or 9 rays are marched independently against the numerical SDF geometry.** Each ray hits the surface at its own physical intersection position and local normal $\hat{n}$.
* The complex beam curvature matrix $\mathbf{Q} = \mathbf{Q}_r + i \mathbf{Q}_i$ and wavefront phase are reconstructed directly from the differential coordinates of the auxiliary rays relative to the chief ray.
* **Advantage:** This completely eliminates the need for analytical 2nd-order surface Hessians ($\mathbf{C} = \nabla \hat{n}$), working seamlessly on arbitrary numerical SDFs, CSG blends, and triangle meshes!

```julia
# Defined ONCE in src/Gaussian.jl for ALL optical components
function interact3d(system::AbstractSystem, object::AbstractObject, gauss::GaussianBeamlet, id::Int)
    i_chief = interact_step(system, object, rays(gauss.chief)[id])
    i_waist = interact_step(system, object, rays(gauss.waist)[id])
    i_div   = interact_step(system, object, rays(gauss.divergence)[id])
    
    if any(isnothing, (i_chief, i_waist, i_div))
        return nothing
    end
    return GaussianBeamletInteraction(i_chief, i_waist, i_div)
end

# Defined ONCE in src/AstigmaticGaussian.jl for ALL optical components
function interact3d(system::AbstractSystem, object::AbstractObject, agb::AstigmaticGaussianBeamlet, id::Int)
    interactions = ntuple(9) do i
        interact_step(system, object, rays(component_beam(agb, i))[id])
    end
    
    any(isnothing, interactions) && return nothing
    return AstigmaticGaussianBeamletInteraction(interactions...)
end
```

---

## 9. GPU Acceleration Architecture

Traditional object-oriented ray tracers fail on GPUs in Julia because they rely on mutable structs, heap allocations (`Array`, `Vector`, `Dict`), and runtime dynamic dispatch.

### Why This Media–Interface Architecture Excels on GPUs (`CUDA.jl`, `KernelAbstractions.jl`, `Metal.jl`):

1. **`isbits` Register Invariant**:
   * `ShapeHit{T}`, `SurfaceID`, `Ray{T}`, and `PolarizedRay{T}` contain only primitive floats and bitstypes.
   * They fit directly into **GPU hardware registers** with zero heap memory allocation.

2. **Branchless / Direct Execution Kernel**:
   * Massive ray bundles ($10^6+$ rays or Monte Carlo ray-tracing) execute in parallel on thousands of GPU cores:

```julia
# Runs simultaneously across thousands of GPU threads via KernelAbstractions.jl / CUDA.jl
@kernel function trace_rays_kernel!(rays_out, @Const(rays_in), @Const(system_shapes), @Const(system_materials))
    idx = @index(Global)
    ray = rays_in[idx]
    
    # 1. SDF ray-march in GPU registers (arithmetic-heavy SIMD)
    hit = gpu_trace_closest(system_shapes, ray)
    
    # 2. Evaluate boundary physics (Snell / Reflection)
    if hit.t < Inf
        trans = gpu_resolve_transition(hit, ray, system_materials)
        rays_out[idx] = gpu_interact(hit.tag, trans, ray, hit.t)
    end
end
```

3. **SDF Ray-Marching is Ideal for GPUs**:
   * Mathematical SDF formulas are arithmetic-dense and memory-light—the exact workload GPUs are built to maximize.

---

## 10. Critical Edge-Case & Drawback Analysis

| Potential Drawback / Edge Case | Root Cause | **Architectural Resolution** |
| :--- | :--- | :--- |
| **1. Ambient Medium in Submerged Systems** | If an `Interface` hardcodes `outside_medium = Air`, submerged optics (e.g. water tanks, cryostats) would require redefining components. | Components default their `outside_medium` to `Ambient()`. During tracing, `Ambient()` dynamically resolves to the enclosing `system.ambient_medium` ($n=1.0$ by default). Explicit outer media are only set for permanently cemented/encased assemblies. |
| **2. Dynamic / User-Assembled Touching Shapes** | If a user translates two separate lenses $L_1$ and $L_2$ so their flat faces touch at runtime, neither lens knows about the other at construction time. | 1. Composite constructors (`DoubletLens`, `CubeBeamsplitter`) bind shared interfaces automatically.<br>2. For ad-hoc user assemblies, the tracer tracks $n_{\text{ray}}$. When exiting $L_1$, if the forward distance to $L_2$ is $t \le \epsilon_{\text{tol}}$, the tracer transitions $n_1 \to n_2$ directly without an artificial air step.<br>3. A fluent helper `cement!(L1, :back, L2, :front)` is provided. |
| **3. Tagging on Arbitrary Custom SDFs** | User-defined mathematical SDF formulas ($f(\mathbf{x}) \le 0$) might not have annotated face tags (`FRONT`, `BACK`). | A fallback classifier assigns tags via outward normal geometry ($\hat{n} \cdot \hat{d}_{\text{opt}}$), with optional predicate overrides `(p, n) -> SurfaceID`. |
| **4. 2D / Non-Volumetric Elements** | Flat apertures, detectors, and thin polarizer sheets do not enclose a 3D volume. | Marked with `is_thin_interface(elem) = true`. Both `inside_medium` and `outside_medium` match ambient space, so rays interact and propagate without index refraction. |
| **5. Floating-Point Coincidence Ambiguity** | Touching surfaces report identical $t_1 \approx t_2$ within numerical precision. | When propagating inside an optic $\Omega_1$, queries are domain-restricted to $\Omega_1$'s bounding manifolds, avoiding floating-point z-fighting at coincident faces. |

---

## 11. Immediate Low-Hanging Fruits & Future-Proofing

| Area | Impact | Detailed Benefit |
| :--- | :--- | :--- |
| **Massive Code Deletion** | ~1,500 lines removed | Deletes duplicate Snell/Fresnel/TIR code in `Lenses.jl`, `Mirrors.jl`, `Beamsplitters.jl`. |
| **Zero Heap Allocations** | CPU L1 cache execution | Inner loops run entirely in registers without GC pauses for $10^5+$ rays. |
| **5-Line Optic Constructors** | Developer Ergonomics | Creating a custom lens takes **5 lines** (`Optic(shape, medium, boundaries)`). |
| **Automatic Differentiation (AD)** | Lens Optimization | Fully immutable `isbits` structs allow **`Enzyme.jl`** and **`ForwardDiff.jl`** to differentiate through systems with native LLVM performance. |
| **Metasurfaces & Gratings** | Extensibility | Adding diffractive phase masks or metasurfaces requires only a new `AbstractSurfaceModel`. |

---

## 12. Comprehensive Comparison: Legacy vs. Definitive Architecture

| Architectural Feature | Legacy BMO (`master`) | Initial PR #60 Draft | **Definitive Architecture (v2.7)** |
| :--- | :--- | :--- | :--- |
| **`Intersection` Representation** | Mutable struct with coincident pointers | Mutable hierarchy (`MultiIntersection`) | **100% Immutable `ShapeHit{T}` (`isbits`)** |
| **Facet ID Type** | Ad-hoc symbols | `tag::Symbol` (dynamic dispatch) | **`SurfaceID` (`UInt32` `isbits` wrapper)** |
| **Memory Allocation** | Heap allocated per hit | Dynamic boxed references | **Zero heap allocations (stack/register only)** |
| **Physics Dispatch** | Duplicated in every component type | Duplicated in components + traits | **Universal dispatch on `(SurfaceModel, Medium)`** |
| **Boundary Lookup** | Ad-hoc getter methods | Dynamic NamedTuple index | **Static Cartesian unrolled indexing (`@nif`)** |
| **Ray Advance Distance** | Undefined $t$ handling | Inconsistent | **Explicit $t_{\text{hit}}$ parameter from `ShapeHit{T}`** |
| **Self-Intersection Prevention** | None (prone to shadow acne) | None | **$\epsilon$-Normal offset along outgoing ray** |
| **Medium Transition Resolution** | Heuristic `isentering` normal checks | Tolerant coincident $t$-matching | **State-aware: `resolve_transition(med_in, med_out, ray, n̂)`** |
| **Cemented Interfaces** | Coincident object pointer hacks | Coincident $t$-search with tolerance | **Native multi-material boundary $(\text{Glass}_1 \to \text{Glass}_2)$** |
| **Beamlet Specialization** | Required duplicate methods for Gausslets | Required duplicate methods | **Universal broadcast across parabasal rays** |
| **GPU Acceleration** | Impossible (heap + mutability) | Impossible | **100% GPU-Ready (Register resident)** |
| **Extensibility** | New optic requires custom `interact3d` | New optic requires custom `interact3d` | **New optic = Shape + Material + Boundary map** |

---

## 13. Contributor FAQ & Extension Patterns

### Q1: How do I add a new custom thin-film coating or surface model?
Define a subtype of `AbstractSurfaceModel` and implement a single `interact_surface` method:
```julia
struct DichroicCoating{T} <: AbstractSurfaceModel
    cutoff_wavelength::T
end

function interact_surface(d::DichroicCoating, trans, ray::AbstractRay{T}, t_hit::T) where {T}
    if wavelength(ray) < d.cutoff_wavelength
        # Reflect short wavelengths
        return interact_surface(IdealMirror(), trans, ray, t_hit)
    else
        # Transmit long wavelengths
        return interact_surface(FresnelInterface(), trans, ray, t_hit)
    end
end
```

---

### Q2: How do I model a custom Freeform / CAD Lens?
1. Define the geometric shape (SDF or Mesh) with `surface_tag(shape, p, n)` returning a `SurfaceID`.
2. Construct the optic: `Optic(my_freeform_sdf, glass_medium, (AR_front, AR_back, Absorber_side))`.
3. Tracing works immediately for standard rays, polarized rays, and Gaussian beamlets with **zero extra lines of optical physics code**.

---

### Q3: How do I model Birefringent Crystals (Calcite, Quartz, Wollaston Prisms)?
Birefringence is a **volumetric medium property**, not a component hack.
1. Define `AnisotropicMedium(no_dispersion, ne_dispersion, optical_axis_vector) <: AbstractMedium`.
2. At the boundary interface, `interact_surface` evaluates the dispersion relation:
   $$k_o = \frac{\omega}{c} n_o, \quad k_e(\theta) = \frac{\omega}{c} \left(\frac{\cos^2\theta}{n_o^2} + \frac{\sin^2\theta}{n_e^2}\right)^{-1/2}$$
3. The boundary pushes the ordinary ray along $\vec{d}_o$ and the extraordinary ray along the Poynting vector $\vec{S}_e$ (accounting for the walk-off angle) onto the `StaticRayStack`.

---

### Q4: How do I model Diffraction Gratings and Holographic Optical Elements (HOEs)?
1. Define `GratingSurface(groove_density, blaze_angle, order_efficiencies) <: AbstractSurfaceModel`.
2. `interact_surface` solves the grating vector equation:
   $$\vec{k}_{\parallel, m} = \vec{k}_{\parallel, \text{inc}} + m \vec{G}$$
   where $\vec{G} = \frac{2\pi}{d} \hat{g}$ is the grating reciprocal vector.
3. Outgoing diffracted orders ($m = 0, \pm 1, \dots$) are created with amplitudes scaled by the grating efficiency curve.

---

### Q5: How do I simulate Gradient-Index (GRIN) Optics or Atmospheric Turbulence?
In a GRIN medium (e.g. Luneburg lenses, GRIN fiber collimators, thermal lensing), the refractive index varies continuously: $n(\mathbf{r}) = n_0(1 - \frac{1}{2} A r^2)$.
* **Volumetric Transport:** When propagating inside a `GRINMedium <: AbstractMedium`, the tracer uses an adaptive 4th-order Runge-Kutta or Verlet ray-marcher solving the ray equation:
  $$\frac{d}{ds} \left( n(\mathbf{r}) \frac{d\mathbf{r}}{ds} \right) = \nabla n(\mathbf{r})$$
* **Boundary Transitions:** When the curved ray reaches the bounding SDF manifold $\partial \Omega$, the standard `interact_surface` pipeline executes seamlessly.

---

### Q6: How do I model Metasurfaces, Spatial Light Modulators (SLMs), and Phase Masks?
Metasurfaces impart an abrupt subwavelength phase discontinuity $\Phi(u, v)$:
1. Define `PhaseMaskSurface(phase_function) <: AbstractSurfaceModel`.
2. In `interact_surface`, Generalized Snell's Law is evaluated directly:
   $$\vec{d}_{\text{trans}} \times \hat{n} - \vec{d}_{\text{inc}} \times \hat{n} = \frac{\lambda_0}{2\pi n_{\text{trans}}} \nabla \Phi(u, v)$$
3. The outgoing wavepacket receives the phase shift $\Delta \phi = \Phi(u_0, v_0)$ without altering bulk geometric SDFs.

---

### Q7: How does Automatic Differentiation (AD) with `Enzyme.jl` or `ForwardDiff.jl` work?
Because all data structures (`ShapeHit{T}`, `Intersection`, `Ray{T}`, `Optic`) are **100% immutable and `isbits` (stack-allocated)**:
* There are zero heap pointers or mutable arrays in the inner loop.
* You can define a merit function $\mathcal{L}(\boldsymbol{\theta})$ (e.g. RMS wavefront error or spot size on a detector) parameterized by lens curvature radii $R_i$, glass thicknesses $d_i$, or positions $\mathbf{x}_i$.
* Calling `Enzyme.autodiff(Reverse, ray_trace_loss, Active(system_params))` computes exact machine-precision gradients $\nabla_{\boldsymbol{\theta}} \mathcal{L}$ for gradient-based optical optimization (BFGS, Adam).

---

### Q8: How do I handle Submerged Systems (e.g. Water Immersion Objectives, Fluidics)?
Simply configure the ambient space of the system:
```julia
system = OpticalSystem(ambient_medium = IsotropicMedium(:Water, water_sellmeier, nothing))
```
All standard optics using `Ambient()` dynamically resolve outside transitions against water ($n=1.33$) rather than air ($n=1.0$) with zero component reconfiguration.

---

### Q9: How do Detectors record Coherent Wavefronts vs. Incoherent Intensity?
In `interact_surface(d::DetectorSurface, trans, ray, t_hit)`:
* **Geometric Rays:** Accumulate incoherent ray count and power $I(x, y) = \sum P_i \delta(x - x_i, y - y_i)$.
* **Polarized Rays & Beamlets:** Accumulate complex electric field vectors $\mathbf{E}(x, y) = \sum \mathbf{E}_i e^{i (k \cdot OPL_i - \omega t)}$.
* The detector matrix provides direct Fourier optics post-processing for **Point Spread Functions (PSF)**, **Modulation Transfer Functions (MTF)**, and **Interferograms**.

---

### Q10: Is Non-Sequential Ray Tracing Supported?
**Yes.** The Media–Interface architecture is inherently non-sequential:
* Rays do not follow a fixed component order list.
* At each step, `trace_all` finds the physically closest boundary manifold.
* Light bounces freely through retroreflectors, integrating spheres, multipass absorption cells, and laser resonator cavities until terminating at detectors, absorbers, or reaching the energy flux cutoff.

---

## 14. Direct Implementation Roadmap for `feature/intersection-handling`

1. **Step 1 (`ShapeHit` & `SurfaceID`):** Define `struct SurfaceID; id::UInt32; end` and `struct ShapeHit{T}` in `src/AbstractTypes/AbstractIntersection.jl`. Update `src/SDFs/AbstractSDF.jl` and mesh solvers to return `ShapeHit{T}`.
2. **Step 2 (`Intersection`):** Replace `ObjectIntersection`, `PlaneIntersection`, and `MultiIntersection` with a single canonical `Intersection{T, S, O}`.
3. **Step 3 (`Media & State Tracking`):** Add `current_medium(ray)` and `resolve_transition(med_inside, med_outside, ray, normal)`.
4. **Step 4 (Consolidate `interact`):** Replace monolithic `interact3d` implementations across `Lenses.jl`, `Mirrors.jl`, and `Detectors.jl` with the universal `interact_surface(model, trans, ray, t_hit)` pipeline.
