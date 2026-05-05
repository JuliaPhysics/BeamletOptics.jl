```@setup astigmatic_beamlet
include(joinpath(@__DIR__, "..", "..", "assets", "cond_save.jl"))

beam_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

conditional_include(joinpath(beam_showcase_dir, "agb_ds_showcase.jl"), use_placeholder=true)
```

# Astigmatic polarized beamlets

Standard stigmatic beamlets are limited to symmetric systems. For modeling arbitrary optical systems—including off-axis components, tilted lenses, or complex media like atmospheric turbulence—the package provides the [`AstigmaticGaussianBeamlet`](@ref). This type implements the **Parabasal Gaussian Beamlet** formalism TODO: Add citations

## Formalism

An astigmatic beamlet is represented by a cluster of **9 rays**:
1.  **Chief Ray**: A [`PolarizedRay`](@ref) defining the central path and polarization state.
2.  **Eight Auxiliary Rays**: Four "positional" rays (waist $x+, x-, y+, y-$) and four "directional" rays (divergence $x+, x-, y+, y-$) that track the complex curvature matrix $\mathbf{Q}(z)$ through the system.

The reduced complex amplitude of the Gaussian beamlet is computed as:
```math
\psi(\mathbf{r}) = \frac{E_0}{\sqrt{\mathbf{h}_1 \times \mathbf{h}_2}} \cdot \exp \left( i k \frac{(\mathbf{h}_1 \times \mathbf{r})(\mathbf{u}_2 \cdot \mathbf{r}) - (\mathbf{h}_2 \times \mathbf{r})(\mathbf{u}_1 \cdot \mathbf{r})}{2 \mathbf{h}_1 \times \mathbf{h}_2} \right)
```
where $\mathbf{h}_{1,2}$ and $\mathbf{u}_{1,2}$ are the two-dimensional complex ray height and angle vectors on the plane perpendicular to the chief ray. The phase factor due to the optical path length of the central ray ($e^{i k (z + \Delta L)}$) is added to this reduced field during the final superposition of multiple beamlets. 

For the paraxial Gaussian beam to be unique, the complex ray vectors must satisfy the **vanishing complex optical invariant**:
```math
\mathbf{h}_1 \times \mathbf{u}_2 - \mathbf{h}_2 \times \mathbf{u}_1 = 0
```
This ensures that the complex curvature matrix $\mathbf{Q}$ is symmetric, a fundamental requirement for physical Gaussian beams.

## The Curvature Matrix $\mathbf{Q}$

The relationship between the auxiliary rays and the complex curvature of the beam is defined by the matrix equation $\mathbf{Q} = \mathbf{U}\mathbf{H}^{-1}$. By arranging the transverse components of the complex rays into $2 \times 2$ matrices $\mathbf{H} = [\mathbf{h}_1, \mathbf{h}_2]$ and $\mathbf{U} = [\mathbf{u}_1, \mathbf{u}_2]$, we can solve for the individual components of $\mathbf{Q}$:

```math
\mathbf{Q} = \frac{1}{\mathbf{h}_1 \times \mathbf{h}_2} \begin{bmatrix} 
u_{1x} h_{2y} - u_{2x} h_{1y} & u_{2x} h_{1x} - u_{1x} h_{2x} \\
u_{1y} h_{2y} - u_{2y} h_{1y} & u_{2y} h_{1x} - u_{1y} h_{2x}
\end{bmatrix}
```

The diagonal elements $Q_{xx}$ and $Q_{yy}$ describe the phase curvature (and beam width) along the primary axes, while the off-diagonal elements $Q_{xy}$ (which must equal $Q_{yx}$ for a paraxial beam) describe the **general astigmatism** or twist of the wavefront. This representation allows the package to track complex, rotating beam profiles as they propagate through non-orthogonal optical systems.

!!! note "Analytic Bilinear Form"
    All vector operations in the formulas above ($\cdot, \times$) and the matrix product $\mathbf{r}^T \mathbf{Q} \mathbf{r}$ use the **analytic bilinear form** rather than the conjugated dot product. This preservation of analytic continuity is essential for the stability of Gaussian beamlets in the complex domain.

## API and Usage

### Construction

You can construct a single astigmatic beamlet with a specified waist and polarization:
```julia
using BeamletOptics
# A 10mm waist beamlet propagating along Y, polarized along X
agb = AstigmaticGaussianBeamlet([0,0,0], [0,1,0], 1064e-9, 10e-3; E0 = [1,0,0])
```

### Field Calculation

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

## Key Assumptions
- **Paraxiality**: The auxiliary rays must remain within the paraxial regime relative to the chief ray.
- **Analytic Continuity**: Quadratic phase calculations use the bilinear dot product to avoid unphysical conjugation of complex ray parameters.
- **Power Conservation**: The decomposition methods automatically scale amplitudes to preserve total field energy regardless of grid resolution.

## Example: Double Slit Diffraction

One of the most powerful applications of the AGB framework is the simulation of coherent diffraction from complex apertures. By using [`WavefrontBeamletDecomposition`](@ref), any arbitrary aperture mask (including phase-varying masks) can be decomposed into a set of coherent beamlets.

In the example below, a classic double-slit mask (10 µm width, 100 µm separation) is decomposed and propagated 0.2 meters. The simulation perfectly reproduces the textbook Fraunhofer diffraction pattern, showcasing both the individual slit diffraction (envelope) and the coherent interference fringes.

![Double Slit Experiment Visualization](agb_doubleslit_experiment.png)

This example demonstrates the package's ability to maintain phase coherence across thousands of beamlets—a prerequisite for high-fidelity modeling of systems like Atmospheric Turbulence or Phase-Modulating SLMs.