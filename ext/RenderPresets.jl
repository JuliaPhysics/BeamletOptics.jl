#=
Material classes of the components, see `_MATERIALS`. The parts of MultiShape objects (e.g. the
prisms and the coating of a `CubeBeamsplitter`) are rendered with the class of each part.
=#

_material_class(::BMO.AbstractObject) = :mechanics
_material_class(::BMO.AbstractRefractiveOptic) = :refractive
_material_class(::BMO.AbstractReflectiveOptic) = :reflective
_material_class(::ThinBeamsplitter) = :coating
_material_class(::BMO.AbstractJonesPolarizer) = :polarizer
_material_class(::BMO.AbstractDetector) = :detector
