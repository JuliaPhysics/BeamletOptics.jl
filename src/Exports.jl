# kinematic export
export translate3d!, translate_to3d!, rotate3d!, xrotate3d!, yrotate3d!, zrotate3d!,
       align3d!, reset_translation3d!, reset_rotation3d!, set_pivot3d!
export position, direction, orientation

# ray and beam type export
export Ray, PolarizedRay, Beam, PointSource, CollimatedSource, UniformDiscSource, UniformPointSource, set_num_rays!,
       GaussianBeamlet, AstigmaticGaussianBeamlet, rayleigh_range, rays, point_on_beam,
       normal3d
export CollimatedGaussianBeamletSource, GaussianBeamletDecomposition,
       SphericalGaussianBeamletSource, EllipticalGaussianBeamletSource, WavefrontBeamletDecomposition,
       GaussianModeDecomposition, AstigmaticBeamGroup

# system
export System, StaticSystem, solve_system!

# object group
export ObjectGroup

# display
export properties, default_properties

# additional
export DiscreteRefractiveIndex, SellmeierEquation

#=
components
=#

# mirrors
export Mirror, SquarePlanoMirror2D, RectangularPlanoMirror, SquarePlanoMirror,
       RoundPlanoMirror, SphericalMirror, RightAnglePrismMirror,
       ConicMirror, OffAxisConicMirror,
       ParabolicMirror, OffAxisParabolicMirror,
       EllipsoidalMirror, OffAxisEllipsoidalMirror,
       HyperbolicMirror, OffAxisHyperbolicMirror

# lenses
export Lens, DoubletLens, ThinLens, SphericalLens, SphericalDoubletLens, thickness,
       TripletLens, SphericalTripletLens

# surfaces
export CircularFlatSurface, RectangularFlatSurface, SphericalSurface, EvenAsphericalSurface,
       CylindricalSurface, AcylindricalSurface

# prisms
export Prism, RightAnglePrism

# detectors
export Detector, electric_field, intensity, spot_diagram, optical_power, gauss_parameters,
       waist_parameters, Centroid, MinMax

# splitters
export ThinBeamsplitter, RoundThinBeamsplitter, RectangularPlateBeamsplitter,
       RoundPlateBeamsplitter, CubeBeamsplitter, RectangularCompensatorPlate

# polarizing components
export PolarizationFilter, RoundPolarizationFilter, LinearPolarizer, RoundLinearPolarizer,
       transmission_axis

# dummies
export NonInteractableObject, MeshDummy, IntersectableObject

# misc
export Retroreflector, get_invariant_threshold, set_invariant_threshold!,
    get_sdf_surface_threshold, get_sdf_raymarch_eps, get_sdf_inside_step,
    get_internal_reflection_threshold, get_line_plane_intersection_threshold,
    get_orthogonality_threshold, get_default_r_max, get_default_depth_max,
    get_default_wavelength, get_default_waist, get_default_power,
    get_progress_threshold, set_progress_threshold!

# render
export render!, get_view, set_view, hide_axis, set_orthographic, arrow!, render_lcs!, look_at!
export live_render!, update_render!, remove_render!, pick_object, studio_lighting!, set_render_look

# render handle protocol and developer API for packages built on BeamletOptics: public, not exported
public AbstractRenderHandle, AbstractObjectRenderHandle, AbstractSystemRenderHandle,
       AbstractBeamRenderHandle, rendered, render_plots, render_children, render_parent,
       render_settings, render_settings!, pickable_plots, look_colors
public is_static, hit_count, wavelength, min_num_rays, ProgressSink, PROGRESS_SINK,
       progress_state, is_cancelled
# re-emitting components (e.g. solvers coupled through a field)
public relaunch!, beamlet_hit_field, GaussianBeamletHit, AstigmaticGaussianBeamletHit
public AbstractSampling, NoSampling, DiscRings, DiscSunflower, ConeRings, ConeSunflower,
       source_beams, sampling_basis
