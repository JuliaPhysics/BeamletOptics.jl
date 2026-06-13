# Thin Beamsplitter Component
# Reuses the generic Coating system under the hood with a SimpleBeamsplitterCoating.

"""
    ThinBeamsplitter(shape; reflectance=0.5)
    ThinBeamsplitter(width, height; reflectance=0.5)
    RoundThinBeamsplitter(diameter; reflectance=0.5)

Creates a zero-thickness, lossless, non-polarizing 2D rectangular or circular beamsplitter.
"""
function ThinBeamsplitter(shape::AbstractShape{T}; reflectance::Real = 0.5) where {T}
    if reflectance >= 1 || reflectance <= 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    r = sqrt(reflectance)
    t = sqrt(1 - reflectance)
    model = SimpleBeamsplitterCoating(r, r, t, t) # rs, rp, ts, tp
    return Coating(shape, model)
end

function ThinBeamsplitter(width::Real, height::Real; reflectance::Real = 0.5)
    return ThinBeamsplitter(RectangularFlatMesh(width, height); reflectance)
end

function ThinBeamsplitter(width::Real; reflectance::Real = 0.5)
    return ThinBeamsplitter(width, width; reflectance)
end

"""
    RoundThinBeamsplitter(diameter; reflectance=0.5)

Creates a circular [`ThinBeamsplitter`](@ref).
"""
function RoundThinBeamsplitter(diameter::Real; reflectance::Real = 0.5)
    return ThinBeamsplitter(CircularFlatMesh(diameter / 2); reflectance)
end

# Check validity of the split coefficients
function Base.isvalid(c::Coating{T, S, SimpleBeamsplitterCoating}) where {T, S}
    model = c.model
    return abs2(model.rs) + abs2(model.ts) ≈ 1
end
