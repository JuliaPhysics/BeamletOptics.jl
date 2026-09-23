"""
    DifferenceSDF{T, S, TT <: Tuple} <: AbstractCompositeSDF{T}

Represents the boolean subtraction of one or more `tools` from a `base` [`AbstractSDF`](@ref),
i.e. `base \\ (tool_1 ∪ tool_2 ∪ …)`. Unlike [`UnionSDF`](@ref), the constituent SDFs are
**required** to overlap: the tools must intersect the base for anything to be removed.

The intended way to construct these is not explicitly but by subtracting `AbstractSDF`s
using the regular `-` operator:

```julia
s1 = SphereSDF(1.0)
s2 = SphereSDF(0.5)
translate3d!(s2, Point3(0.5, 0.0, 0.0))

# will result in a sphere with a smaller, off-center sphere carved out of it
s_diff = s1 - s2
```

# Non-commutativity and evaluation order

Subtraction is neither commutative nor associative, so `a - b` and `b - a` differ, and a
chain like `a + b - c + d` must keep track of *when* each operand was combined. Repeated
subtraction is flattened along its own chain via the identity

```math
(A \\setminus B) \\setminus C = A \\setminus (B \\cup C)
```

so `a - b - c` produces a single `DifferenceSDF` with one `base` and a flat `tools` tuple,
rather than deepening the type on every chained `-`. Mixed `+`/`-` chains still *nest*
(e.g. `a + b - c` is `(a + b) - c`, and `... + d` afterwards wraps the whole difference in a
[`UnionSDF`](@ref)) — this is what keeps their meaning, since material added after a
subtraction must not be carved away by it.

# Exactness

Per [iquilezles.org/articles/distfunctions](https://iquilezles.org/articles/distfunctions/),
`opSubtraction = max(-a, b)` is only a *bound*, not an exact SDF: it under-estimates the true
distance near the seam between `base` and `tool`. This is the *safe* direction for sphere
tracing — an under-estimate never causes tunneling, it only costs a few extra ray marching
iterations near concave creases. Correctness of ray marching hinges on [`normal3d`](@ref)
being resolved by explicit operand selection (as implemented here) rather than by automatic
differentiation, which picks the wrong sub-shape at the seam.
"""
mutable struct DifferenceSDF{T <: Number, S <: AbstractSDF{T},
        TT <: Tuple{Vararg{AbstractSDF{T}}}} <: AbstractCompositeSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    base::S
    tools::TT
end

operands(d::DifferenceSDF) = (d.base, d.tools...)

function DifferenceSDF{T}(base::AbstractSDF{T}, tools::Tuple{Vararg{AbstractSDF{T}}}) where T
    _check_kinematic_members((base, tools...))
    DifferenceSDF{T, typeof(base), typeof(tools)}(
        SMatrix{3,3}(one(T)*I),
        SMatrix{3,3}(one(T)*I),
        Point3{T}(zero(T)),
        base,
        tools
        )
end

"""
    thickness(difference)

Calculates the thickness of a [`DifferenceSDF`](@ref) as the thickness of its `base` —
removing material cannot increase the axial extent.
"""
function thickness(d::DifferenceSDF{T}) where T
    if hasmethod(thickness, Tuple{typeof(d.base)})
        return thickness(d.base)
    end
    return zero(T)
end

"""
    bounding_sphere(d::DifferenceSDF)

Returns the [`bounding_sphere`](@ref) of `d.base` — the result of a subtraction is always a
subset of the base, so the base's bounding sphere (or `nothing`, if it has none) also bounds
`d`.
"""
bounding_sphere(d::DifferenceSDF) = bounding_sphere(d.base)

function sdf(d::DifferenceSDF, pos)
    # sdf to world transform handled by sub-SDFs
    return max(sdf(d.base, pos), maximum(-sdf(t, pos) for t in d.tools))
end

Base.:-(s1::AbstractSDF{T}, s2::AbstractSDF{T}) where T = DifferenceSDF{T}(s1, (s2,))
Base.:-(d::DifferenceSDF{T}, s::AbstractSDF{T}) where T = DifferenceSDF{T}(d.base, (d.tools..., s))

# NOTE: no unary `Base.:-(::AbstractSDF)` is defined on purpose. The complement of a
# bounded solid is unbounded, which would break both ray marchers (they assume a finite
# object to step towards/away from) and the `bounding_box` fallback (which probes a fixed
# ±1000 m range around the shape).

# Without this function it is not possible for SDFs encapsulated in a DifferenceSDF
# to specialize normal3d as always the generic normal3d function is called. Do NOT fall
# through to normal_fd/AD: AD through a composite max/min picks the wrong sub-shape at the
# seam (see MeniscusLensSDF.jl), and the numeric_gradient fallback only triggers on NaN
# with a stencil far wider than the hit epsilon, so it straddles the seam too.
function normal3d(d::DifferenceSDF, pos)
    best = sdf(d.base, pos)
    idx = 0
    for (i, t) in enumerate(d.tools)
        v = -sdf(t, pos)
        if v > best
            best = v
            idx = i
        end
    end
    # On a carved wall, the outward normal of the result is the inward normal of the
    # active tool, hence the sign flip. Returned directly in world coordinates, per the
    # AbstractCompositeSDF convention (no orientation(d) * factor).
    return idx == 0 ? normal3d(d.base, pos) : -normal3d(d.tools[idx], pos)
end
