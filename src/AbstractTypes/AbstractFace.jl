"""
    AbstractFace

A face of a body, e.g. the front surface of a lens, that carries a surface model, see
[`surface_model!`](@ref). A face is a question about a hit rather than a piece of geometry: whether the
hit lies on it, answered in the local frame of the body. Hence a face moves with its body, and it does
not depend on how the shape of the body is represented.

# Implementation reqs.

## Fields

None required.

## Getters/setters

None required.

## Functions

- `isonface(face, part, hit, point, normal)`: returns `true` if the `hit` of the `part` at the
  `point` with its outward `normal`, both in the local frame of the `part`, lies on the `face`. By
  default, the face is the native face of the hit, see [`face`](@ref): a shape that knows its faces,
  e.g. a boundary representation, sets it in its `intersect3d` method, and its faces are coated
  directly. A face that is defined by geometry, e.g. [`FrontFace`](@ref), implements its own method.
"""
abstract type AbstractFace end

