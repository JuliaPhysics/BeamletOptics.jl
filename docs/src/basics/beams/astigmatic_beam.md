```@setup astigmatic_beamlet
include(joinpath(@__DIR__, "..", "..", "assets", "cond_save.jl"))

beam_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

conditional_include(joinpath(beam_showcase_dir, "agb_showcase.jl"), use_placeholder=false)
```

# Astigmatic polarized beamlets

Standard stigmatic beamlets are limited to symmetric systems. For the modeling of arbitrary optical systems, including off-axis components, tilted lenses, or complex media like atmospheric turbulence, the package provides the [`AstigmaticGaussianBeamlet`](@ref). This type implements the astigmatic Gaussian beamlet ray tracing formalism initially proposed by Greynolds (1986) [Greynolds:1986_1, Greynolds:1986_2](@cite).

However, this idea has been referred to in multiple ways over the last few decades, including but not limited to: Gaussian Beamlet Tracing (GBT), Fat Ray Tracing (FRT), Parabasal Ray Tracing (PRT), Complex Ray Tracing (CRT), Gauss Beamlet Propagation (GBP) and Beam Synthesis Propagation (BSP).

## Formalism

Similar to the Gauss model presented in the [Stigmatic beamlets](@ref) chapter, an astigmatic beamlet can represented by a cluster of **9 rays**:

!!! tip "Ray representation"
    1.  **Chief Ray**: a [`PolarizedRay`](@ref) defining the central path and polarization state.
    2.  **Waist Rays**: four "positional" rays (waist $x+, x-, y+, y-$)
    3.  **Divergence Rays**: four "directional" rays (divergence $x+, x-, y+, y-$)

These rays can be thought of tracking the complex curvature matrix $\mathbf{Q}(z)$ through the system. The key assumptions made in order for this formalism to hold are:

!!! info "Key assumptions"
    1. **Paraxiality**: the auxiliary rays must remain within the paraxial regime relative to the chief ray.
    2. **Parabolic interaction**: since only astigmatism can be captured, each surface interaction must be approximately parabolical
    3. **Homogeneous polarization**: the polarization state of the traced field is assumed to be homogeneous over each beamlet
    4. **Matched scale**: the beamlet must be smaller than the optical element it interacts with, see also 2.

The initial ordering of the geometric beams/rays is shown in the figure below. The colors represent the chief (red), waist (blue) and divergence (green) beams.

![Astigmatic ray tracing I](agbtest1.png)

As with the [`GaussianBeamlet`](@ref), tracing these rays through a system allows the reconstruction of the waist envelope and electric field. The reduced complex amplitude of the Gaussian beamlet is computed as

```math
\psi(\mathbf{r}) = \frac{E_0}{\sqrt{\mathbf{h}_1 \times \mathbf{h}_2}} \cdot \exp \left( i k \frac{(\mathbf{h}_1 \times \mathbf{r})(\mathbf{u}_2 \cdot \mathbf{r}) - (\mathbf{h}_2 \times \mathbf{r})(\mathbf{u}_1 \cdot \mathbf{r})}{2 \mathbf{h}_1 \times \mathbf{h}_2} \right) \,,
```
where $\mathbf{h}_{1,2}$ and $\mathbf{u}_{1,2}$ are the two-dimensional complex ray height and angle vectors on the plane perpendicular to the chief ray [Greynolds:1986_1, Wilhelm:2001, Greynolds:2014](@cite). By also tracing a 3D-field vector along the chief ray based on the formalism introduced in the [Polarized rays](@ref) section, polarization effects can be considered when calculating $\psi$ [Worku:2017](@cite). For a simple lens setup, the resulting waist is illustrated in the image below.

![Astigmatic ray tracing II](agbtest2.png)

While the above beam closely matches the example Gaussian given in the previous section ([Obtaining the beam parameters](@ref)), the real advantage of this extended method lies in the tracing of non-symmetric systems. For a tilted two cylinder lens system, the following result with strong astigmatic effects is obtained. Multiple waist slices are marked with red dots.

![Astigmatic ray tracing III](agbtest3.png)

The phase factor due to the optical path length of the central ray ($e^{i k (z + \Delta L)}$) is added to the reduced field $\psi$ during the final field calculation for each beamlet. More information on the mathemathics behind this formalism can be found below in [The Curvature Matrix $\mathbf{Q}$](@ref) section.

!!! info "Optical invariant"
    In order to ensure the correctness of the traced beamlet, the complex ray vectors must satisfy the **vanishing complex optical invariant**:
    ```math
    \mathbf{h}_1 \times \mathbf{u}_2 - \mathbf{h}_2 \times \mathbf{u}_1 = 0
    ```
    This ensures that the complex curvature matrix $\mathbf{Q}$ is symmetric, since the beamlet formalism can not capture higher-order abberations.

## Astigmatic beamlets

You can construct a single astigmatic beamlet with a specified waist and polarization:

```@docs; canonical=false
AstigmaticGaussianBeamlet(::AbstractArray, ::AbstractArray, ::Real, ::Real)
```

The constructor will spawn an `AstigmaticGaussianBeamlet` which is implemented as follows:

```@docs; canonical=false
AstigmaticGaussianBeamlet
```

### Calculating astigmatic beam parameters

To obtain the physical vector electric field [V/m] at any point:
```julia
# Field at world position (x, y, z)
E_vec = polarized_field(agb, [0.01, 0, 0], 10.0)
```

```julia
# Decompose a field from a detector into a new AstigmaticBeamGroup
beams = WavefrontBeamletDecomposition(x, z, amplitude, phase, dir, λ)

# Solve system (propagate all beamlets to next segment)
solve_system!(system, beams)
```

## The Curvature Matrix $\mathbf{Q}$

The relationship between the auxiliary rays and the complex curvature of the beam is defined by the matrix equation $\mathbf{Q} = \mathbf{U}\mathbf{H}^{-1}$. By arranging the transverse components of the complex rays into $2 \times 2$ matrices $\mathbf{H} = [\mathbf{h}_1, \mathbf{h}_2]$ and $\mathbf{U} = [\mathbf{u}_1, \mathbf{u}_2]$, we can solve for the individual components of $\mathbf{Q}$:

```math
\mathbf{Q} = \frac{1}{\mathbf{h}_1 \times \mathbf{h}_2} \begin{bmatrix} 
u_{1x} h_{2y} - u_{2x} h_{1y} & u_{2x} h_{1x} - u_{1x} h_{2x} \\
u_{1y} h_{2y} - u_{2y} h_{1y} & u_{2y} h_{1x} - u_{1y} h_{2x}
\end{bmatrix}
```

The diagonal elements $Q_{xx}$ and $Q_{yy}$ describe the phase curvature (and beam width) along the primary axes, while the off-diagonal elements $Q_{xy}$ (which must equal $Q_{yx}$ for a paraxial beam) describe the **general astigmatism** or twist of the wavefront. This representation allows the formalism to track complex, rotating beam profiles as they propagate through non-orthogonal optical systems. It is important to note that this matrix is not directly tracked during beamlet propagation in **BMO**. More information can be found in the works of Kochkina (2013) and Ashcraft et al. (2024) [Kochkina:2013, Ashcraft:2024](@cite).

!!! note "Analytic Bilinear Form"
    All vector operations in the formulas above ($\cdot, \times$) and the matrix product $\mathbf{r}^T \mathbf{Q} \mathbf{r}$ use the **analytic bilinear form** rather than the conjugated dot product. This preservation of analytic continuity is essential for the stability of Gaussian beamlets in the complex domain.
