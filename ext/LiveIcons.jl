# Icons and flat icon buttons of the app layout of the live view
#
# Icon paths not marked "own design" are from Material Symbols (outlined, weight 400, grade 0,
# optical size 24) by Google, licensed under the Apache License 2.0:
#   https://github.com/google/material-design-icons (symbols/web/<name>/materialsymbolsoutlined)
#   https://www.apache.org/licenses/LICENSE-2.0
# The path data is used unchanged. The icons marked "own design" are drawn on the same grid.

using Makie: Observable, Point2f, Point2d, RGBAf, BezierPath, MoveTo, LineTo, CurveTo, ClosePath,
             Consume, lift, onany, on, translate!

# SVG path data on the Material Symbols grid: viewBox "0 -960 960 960", y down
const _ICON_SVG = Dict{Symbol, String}(
    # keep (outlined and filled): pin a card
    :pin => "m640-480 80 80v80H520v240l-40 40-40-40v-240H240v-80l80-80v-280h-40v-80h400v80h-40v280Zm-286 80h252l-46-46v-314H400v314l-46 46Zm126 0Z",
    :pinned => "m640-480 80 80v80H520v240l-40 40-40-40v-240H240v-80l80-80v-280h-40v-80h400v80h-40v280Z",
    # refresh
    :trace => "M480-160q-134 0-227-93t-93-227q0-134 93-227t227-93q69 0 132 28.5T720-690v-110h80v280H520v-80h168q-32-56-87.5-88T480-720q-100 0-170 70t-70 170q0 100 70 170t170 70q77 0 139-44t87-116h84q-28 106-114 173t-196 67Z",
    # autoplay
    :auto_trace => "M380-300v-360l280 180-280 180ZM480-40q-108 0-202.5-49.5T120-228v108H40v-240h240v80h-98q51 75 129.5 117.5T480-120q115 0 208.5-66T820-361l78 18q-45 136-160 219.5T480-40ZM42-520q7-67 32-128.5T143-762l57 57q-32 41-52 87.5T123-520H42Zm214-241-57-57q53-44 114-69.5T440-918v80q-51 5-97 25t-87 52Zm449 0q-41-32-87.5-52T520-838v-80q67 6 128.5 31T762-818l-57 57Zm133 241q-5-51-25-97.5T761-705l57-57q44 52 69 113.5T918-520h-80Z",
    # home
    :home => "M240-200h120v-240h240v240h120v-360L480-740 240-560v360Zm-80 80v-480l320-240 320 240v480H520v-240h-80v240H160Zm320-350Z",
    # fit_screen
    :fit => "M800-600v-120H680v-80h120q33 0 56.5 23.5T880-720v120h-80Zm-720 0v-120q0-33 23.5-56.5T160-800h120v80H160v120H80Zm600 440v-80h120v-120h80v120q0 33-23.5 56.5T800-160H680Zm-520 0q-33 0-56.5-23.5T80-240v-120h80v120h120v80H160Zm80-160v-320h480v320H240Zm80-80h320v-160H320v160Zm0 0v-160 160Z",
    # photo_camera
    :views => "M480-260q75 0 127.5-52.5T660-440q0-75-52.5-127.5T480-620q-75 0-127.5 52.5T300-440q0 75 52.5 127.5T480-260Zm0-80q-42 0-71-29t-29-71q0-42 29-71t71-29q42 0 71 29t29 71q0 42-29 71t-71 29ZM160-120q-33 0-56.5-23.5T80-200v-480q0-33 23.5-56.5T160-760h126l74-80h240l74 80h126q33 0 56.5 23.5T880-680v480q0 33-23.5 56.5T800-120H160Zm0-80h640v-480H638l-73-80H395l-73 80H160v480Zm320-240Z",
    # add_a_photo
    :save_view => "M440-440ZM120-120q-33 0-56.5-23.5T40-200v-480q0-33 23.5-56.5T120-760h126l74-80h240v80H355l-73 80H120v480h640v-360h80v360q0 33-23.5 56.5T760-120H120Zm640-560v-80h-80v-80h80v-80h80v80h80v80h-80v80h-80ZM440-260q75 0 127.5-52.5T620-440q0-75-52.5-127.5T440-620q-75 0-127.5 52.5T260-440q0 75 52.5 127.5T440-260Zm0-80q-42 0-71-29t-29-71q0-42 29-71t71-29q42 0 71 29t29 71q0 42-29 71t-71 29Z",
    # view_in_ar
    :orthographic => "M440-181 240-296q-19-11-29.5-29T200-365v-230q0-22 10.5-40t29.5-29l200-115q19-11 40-11t40 11l200 115q19 11 29.5 29t10.5 40v230q0 22-10.5 40T720-296L520-181q-19 11-40 11t-40-11Zm0-92v-184l-160-93v185l160 92Zm80 0 160-92v-185l-160 93v184ZM80-680v-120q0-33 23.5-56.5T160-880h120v80H160v120H80ZM280-80H160q-33 0-56.5-23.5T80-160v-120h80v120h120v80Zm400 0v-80h120v-120h80v120q0 33-23.5 56.5T800-80H680Zm120-600v-120H680v-80h120q33 0 56.5 23.5T880-800v120h-80ZM480-526l158-93-158-91-158 91 158 93Zm0 45Zm0-45Zm40 69Zm-80 0Z",
    # content_cut
    :clip => "M760-120 480-400l-94 94q8 15 11 32t3 34q0 66-47 113T240-80q-66 0-113-47T80-240q0-66 47-113t113-47q17 0 34 3t32 11l94-94-94-94q-15 8-32 11t-34 3q-66 0-113-47T80-720q0-66 47-113t113-47q66 0 113 47t47 113q0 17-3 34t-11 32l494 494v40H760ZM600-520l-80-80 240-240h120v40L600-520ZM240-640q33 0 56.5-23.5T320-720q0-33-23.5-56.5T240-800q-33 0-56.5 23.5T160-720q0 33 23.5 56.5T240-640Zm240 180q8 0 14-6t6-14q0-8-6-14t-14-6q-8 0-14 6t-6 14q0 8 6 14t14 6ZM240-160q33 0 56.5-23.5T320-240q0-33-23.5-56.5T240-320q-33 0-56.5 23.5T160-240q0 33 23.5 56.5T240-160Z",
    # flare
    :sources => "M40-440v-80h240v80H40Zm270-154-84-84 56-56 84 84-56 56Zm130-86v-240h80v240h-80Zm210 86-56-56 84-84 56 56-84 84Zm30 154v-80h240v80H680Zm-200 80q-50 0-85-35t-35-85q0-50 35-85t85-35q50 0 85 35t35 85q0 50-35 85t-85 35Zm198 134-84-84 56-56 84 84-56 56Zm-396 0-56-56 84-84 56 56-84 84ZM440-40v-240h80v240h-80Z",
    # straighten
    :measure => "M160-240q-33 0-56.5-23.5T80-320v-320q0-33 23.5-56.5T160-720h640q33 0 56.5 23.5T880-640v320q0 33-23.5 56.5T800-240H160Zm0-80h640v-320H680v160h-80v-160h-80v160h-80v-160h-80v160h-80v-160H160v320Zm120-160h80-80Zm160 0h80-80Zm160 0h80-80Zm-120 0Z",
    # ios_share
    :export => "M240-40q-33 0-56.5-23.5T160-120v-440q0-33 23.5-56.5T240-640h120v80H240v440h480v-440H600v-80h120q33 0 56.5 23.5T800-560v440q0 33-23.5 56.5T720-40H240Zm200-280v-447l-64 64-56-57 160-160 160 160-56 57-64-64v447h-80Z",
    # dock_to_right
    :panel_left => "M200-120q-33 0-56.5-23.5T120-200v-560q0-33 23.5-56.5T200-840h560q33 0 56.5 23.5T840-760v560q0 33-23.5 56.5T760-120H200Zm120-80v-560H200v560h120Zm80 0h360v-560H400v560Zm-80 0H200h120Z",
    # dock_to_left
    :panel_right => "M200-120q-33 0-56.5-23.5T120-200v-560q0-33 23.5-56.5T200-840h560q33 0 56.5 23.5T840-760v560q0 33-23.5 56.5T760-120H200Zm440-80h120v-560H640v560Zm-80 0v-560H200v560h360Zm80 0h120-120Z",
    # dock_to_bottom
    :panel_bottom => "M200-120q-33 0-56.5-23.5T120-200v-560q0-33 23.5-56.5T200-840h560q33 0 56.5 23.5T840-760v560q0 33-23.5 56.5T760-120H200Zm0-200v120h560v-120H200Zm0-80h560v-360H200v360Zm0 80v120-120Z",
    # help
    :help => "M478-240q21 0 35.5-14.5T528-290q0-21-14.5-35.5T478-340q-21 0-35.5 14.5T428-290q0 21 14.5 35.5T478-240Zm-36-154h74q0-33 7.5-52t42.5-52q26-26 41-49.5t15-56.5q0-56-41-86t-97-30q-57 0-92.5 30T342-618l66 26q5-18 22.5-39t53.5-21q32 0 48 17.5t16 38.5q0 20-12 37.5T506-526q-44 39-54 59t-10 73Zm38 314q-83 0-156-31.5T197-197q-54-54-85.5-127T80-480q0-83 31.5-156T197-763q54-54 127-85.5T480-880q83 0 156 31.5T763-763q54 54 85.5 127T880-480q0 83-31.5 156T763-197q-54 54-127 85.5T480-80Zm0-80q134 0 227-93t93-227q0-134-93-227t-227-93q-134 0-227 93t-93 227q0 134 93 227t227 93Zm0-320Z",
    # visibility
    :eye => "M480-320q75 0 127.5-52.5T660-500q0-75-52.5-127.5T480-680q-75 0-127.5 52.5T300-500q0 75 52.5 127.5T480-320Zm0-72q-45 0-76.5-31.5T372-500q0-45 31.5-76.5T480-608q45 0 76.5 31.5T588-500q0 45-31.5 76.5T480-392Zm0 192q-146 0-266-81.5T40-500q54-137 174-218.5T480-800q146 0 266 81.5T920-500q-54 137-174 218.5T480-200Zm0-300Zm0 220q113 0 207.5-59.5T832-500q-50-101-144.5-160.5T480-720q-113 0-207.5 59.5T128-500q50 101 144.5 160.5T480-280Z",
    # visibility_off
    :eye_off => "m644-428-58-58q9-47-27-88t-93-32l-58-58q17-8 34.5-12t37.5-4q75 0 127.5 52.5T660-500q0 20-4 37.5T644-428Zm128 126-58-56q38-29 67.5-63.5T832-500q-50-101-143.5-160.5T480-720q-29 0-57 4t-55 12l-62-62q41-17 84-25.5t90-8.5q151 0 269 83.5T920-500q-23 59-60.5 109.5T772-302Zm20 246L624-222q-35 11-70.5 16.5T480-200q-151 0-269-83.5T40-500q21-53 53-98.5t73-81.5L56-792l56-56 736 736-56 56ZM222-624q-29 26-53 57t-41 67q50 101 143.5 160.5T480-280q20 0 39-2.5t39-5.5l-36-38q-11 3-21 4.5t-21 1.5q-75 0-127.5-52.5T300-500q0-11 1.5-21t4.5-21l-84-82Zm319 93Zm-151 75Z",
    # chevron_right
    :expand => "M504-480 320-664l56-56 240 240-240 240-56-56 184-184Z",
    # expand_more
    :collapse => "M480-345 240-585l56-56 184 184 184-184 56 56-240 240Z",
    # folder
    :group => "M160-160q-33 0-56.5-23.5T80-240v-480q0-33 23.5-56.5T160-800h240l80 80h320q33 0 56.5 23.5T880-640v400q0 33-23.5 56.5T800-160H160Zm0-80h640v-400H447l-80-80H160v480Zm0 0v-480 480Z",
    # deployed_code
    :mesh => "M440-183v-274L200-596v274l240 139Zm80 0 240-139v-274L520-457v274Zm-40-343 237-137-237-137-237 137 237 137ZM160-252q-19-11-29.5-29T120-321v-318q0-22 10.5-40t29.5-29l280-161q19-11 40-11t40 11l280 161q19 11 29.5 29t10.5 40v318q0 22-10.5 40T800-252L520-91q-19 11-40 11t-40-11L160-252Zm320-228Z",
    # category
    :object => "m260-520 220-360 220 360H260ZM700-80q-75 0-127.5-52.5T520-260q0-75 52.5-127.5T700-440q75 0 127.5 52.5T880-260q0 75-52.5 127.5T700-80Zm-580-20v-320h320v320H120Zm580-60q42 0 71-29t29-71q0-42-29-71t-71-29q-42 0-71 29t-29 71q0 42 29 71t71 29Zm-500-20h160v-160H200v160Zm202-420h156l-78-126-78 126Zm78 0ZM360-340Zm340 80Z",
    # Own designs, in the outlined style of Material Symbols: strokes of 80 units (2 px at 24 px),
    # solid parts clockwise, holes counterclockwise (nonzero fill rule of FreeType)
    # own design: a ray that ends at a plane in perspective, dashed behind it
    :clip_beams => "M500-720L640-880L640-240L500-80ZM80-520L500-520L500-440L80-440ZM700-520L780-520L780-440L700-440ZM820-520L880-520L880-440L820-440Z",
    # own design: biconvex lens with flat edges, the optical axis on both sides
    :lens => "M400-880H560Q840-480 560-80H400Q120-480 400-880ZM440-800Q160-480 440-160H520Q720-480 520-800ZM80-520L200-520L200-440L80-440ZM760-520L880-520L880-440L760-440Z",
    # own design: thick mirror with a ray reflected on it
    :mirror => "M80-240L880-240L880-120L80-120ZM182.7-853L480-429.6L668.1-697.5L733.6-651.5L480-290.4L117.3-807ZM810-830L784.2-567.1L571.4-716.5Z",
    # own design: sensor with 2 × 2 pixels
    :detector => "M120-840L840-840L840-120L120-120ZM200-200L760-200L760-760L200-760ZM280-680L440-680L440-520L280-520ZM520-680L680-680L680-520L520-520ZM280-440L440-440L440-280L280-280ZM520-440L680-440L680-280L520-280Z",
    # own design: laser head emitting a beam
    :source => "M80-700L500-700L500-260L80-260ZM160-340L420-340L420-620L160-620ZM500-520L720-520L720-440L500-440ZM880-480L680-310L680-650Z",
    # own design: beamsplitter cube with its diagonal coating, an incoming ray from the left, the
    # transmitted ray to the right and the reflected ray upwards
    :beamsplitter => "M200-760L760-760L760-200L200-200ZM280-280L680-280L680-680L280-680ZM280-336.6L623.4-680L680-680L680-623.4L336.6-280L280-280ZM40-520L200-520L200-440L40-440ZM760-520L920-520L920-440L760-440ZM440-920L520-920L520-760L440-760Z",
    # own design: polarizer disc with its transmission axis
    :polarizer => "M480-840C678.8-840 840-678.8 840-480C840-281.2 678.8-120 480-120C281.2-120 120-281.2 120-480C120-678.8 281.2-840 480-840ZM480-760C325.4-760 200-634.6 200-480C200-325.4 325.4-200 480-200C634.6-200 760-325.4 760-480C760-634.6 634.6-760 480-760ZM440-920L520-920L520-40L440-40Z",
    # own design: plane in perspective with its normal
    :clip_plane => "M80-120L300-400L880-400L660-120ZM625-180L745-340L335-340L215-180ZM440-670L520-670L520-260L440-260ZM480-880L630-650L330-650Z",
    # own design: bar chart, the icon of the tabs of `add_panel!`
    :chart => "M120-120V-840H200V-200H840V-120ZM280-280V-560H400V-280ZM480-280V-760H600V-280ZM680-280V-480H800V-280Z",
    # own design: a window open at its top right corner and an arrow out of it, i.e. a card that
    # leaves the sidebar and floats in the 3D view
    :float => "M120-840L200-840L200-120L120-120ZM120-200L840-200L840-120L120-120ZM120-840L480-840L480-760L120-760ZM760-480L840-480L840-120L760-120ZM560-840L840-840L840-760L560-760ZM760-840L840-840L840-560L760-560ZM371.7-428.3L771.7-828.3L828.3-771.7L428.3-371.7Z",
    # own design: a window with a sidebar at the right and an arrow into it, i.e. a card that goes
    # back into the sidebar
    :dock => "M120-840L840-840L840-120L120-120ZM200-200L760-200L760-760L200-760ZM560-760L640-760L640-200L560-200ZM220-520L400-520L400-440L220-440ZM360-620L500-480L360-340Z",
    # own design: breadboard with mounting holes, i.e. an optical system
    :system => "M80-760H880V-200H80ZM160-680V-280H800V-680ZM355-560C355-529.6 330.4-505 300-505C269.6-505 245-529.6 245-560C245-590.4 269.6-615 300-615C330.4-615 355-590.4 355-560ZM535-560C535-529.6 510.4-505 480-505C449.6-505 425-529.6 425-560C425-590.4 449.6-615 480-615C510.4-615 535-590.4 535-560ZM715-560C715-529.6 690.4-505 660-505C629.6-505 605-529.6 605-560C605-590.4 629.6-615 660-615C690.4-615 715-590.4 715-560ZM355-400C355-369.6 330.4-345 300-345C269.6-345 245-369.6 245-400C245-430.4 269.6-455 300-455C330.4-455 355-430.4 355-400ZM535-400C535-369.6 510.4-345 480-345C449.6-345 425-369.6 425-400C425-430.4 449.6-455 480-455C510.4-455 535-430.4 535-400ZM715-400C715-369.6 690.4-345 660-345C629.6-345 605-369.6 605-400C605-430.4 629.6-455 660-455C690.4-455 715-430.4 715-400Z",
)

const _SVG_TOKEN = r"[A-Za-z]|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"

"""
    _svg_path(d::AbstractString; grid = 960)

Parses the SVG path data `d` on a `grid × grid` view box with the origin at the bottom left and y
down, i.e. `viewBox = "0 -grid grid grid"` as used by Material Symbols, into a `BezierPath` that
spans the unit square centered at the origin, y up. The view box, not the bounding box of the
path, is mapped to the unit square, such that all icons share the same scale and padding.

Supports the commands `M L H V Q T C S Z`, absolute and relative (no arcs); quadratic curves are
converted to cubic ones. Subpaths without a segment (`m…Z`, which Material Symbols use as a
placeholder) are dropped.
"""
function _svg_path(d::AbstractString; grid::Real = 960)
    tokens = [m.match for m in eachmatch(_SVG_TOKEN, d)]
    commands = Makie.PathCommand[]
    g = Float64(grid)
    tf(p) = Point2d(p[1] / g - 0.5, -p[2] / g - 0.5)
    cur = start = ctrl = Point2d(0, 0)
    cmd = prev = 'M'
    i = 1
    num() = (i += 1; parse(Float64, tokens[i - 1]))
    pt(rel) = (p = Point2d(num(), num()); rel ? cur + p : p)
    while i <= length(tokens)
        if isletter(tokens[i][1])
            cmd = tokens[i][1]
            i += 1
        end
        rel = islowercase(cmd)
        c = uppercase(cmd)
        if c == 'M'
            cur = start = pt(rel)
            push!(commands, MoveTo(tf(cur)))
            # further coordinate pairs are lines
            cmd = rel ? 'l' : 'L'
        elseif c == 'L'
            cur = pt(rel)
            push!(commands, LineTo(tf(cur)))
        elseif c == 'H'
            x = num()
            cur = Point2d(rel ? cur[1] + x : x, cur[2])
            push!(commands, LineTo(tf(cur)))
        elseif c == 'V'
            y = num()
            cur = Point2d(cur[1], rel ? cur[2] + y : y)
            push!(commands, LineTo(tf(cur)))
        elseif c == 'Q' || c == 'T'
            # T reflects the control point of a preceding quadratic curve
            q = c == 'Q' ? pt(rel) : (prev in ('Q', 'T') ? 2cur - ctrl : cur)
            p = pt(rel)
            push!(commands, CurveTo(tf(cur + 2 / 3 * (q - cur)), tf(p + 2 / 3 * (q - p)), tf(p)))
            ctrl, cur = q, p
        elseif c == 'C' || c == 'S'
            # S reflects the second control point of a preceding cubic curve
            c1 = c == 'C' ? pt(rel) : (prev in ('C', 'S') ? 2cur - ctrl : cur)
            c2 = pt(rel)
            p = pt(rel)
            push!(commands, CurveTo(tf(c1), tf(c2), tf(p)))
            ctrl, cur = c2, p
        elseif c == 'Z'
            if last(commands) isa MoveTo
                pop!(commands)
            else
                # An explicit line back to the start: Makie's bounding box of a path, which
                # scales the marker, misses the start of a subpath that is only reached by
                # `ClosePath`, e.g. the tip of an arrow head drawn from its tip
                cur == start || push!(commands, LineTo(tf(start)))
                push!(commands, ClosePath())
            end
            cur = start
        else
            throw(ArgumentError("SVG path command $cmd is not supported"))
        end
        prev = c
    end
    return BezierPath(commands)
end

# Parsed once, at precompile time
const _ICONS = Dict{Symbol, BezierPath}(name => _svg_path(d) for (name, d) in _ICON_SVG)

"""
    _icon(name::Symbol)
    _icon(path::BezierPath)

Returns the icon `name` as a `BezierPath` in the unit square centered at the origin, for use as a
`scatter` marker (`markersize` = size of the icon in pixels, with the usual padding of Material
Symbols). See `keys(_ICONS)` for the available names; `:object` is the generic icon of objects
without an own one. An own `path`, e.g. the icon of a tool of [`add_tool!`](@ref), is returned
as it is.
"""
function _icon(name::Symbol)
    haskey(_ICONS, name) && return _ICONS[name]
    throw(ArgumentError("Unknown icon :$name, available: $(join(sort!(collect(keys(_ICONS))), ", "))"))
end
_icon(path::BezierPath) = path

# Colors of the icon buttons for a light theme; the app layout passes the tokens of its theme
const _ICON_COLOR = RGBAf(0.23, 0.25, 0.29, 1)
const _ICON_HOVER_COLOR = RGBAf(0.1, 0.12, 0.16, 0.08)
const _ICON_ACTIVE_COLOR = RGBAf(0.16, 0.44, 0.9, 0.16)
const _ICON_ACTIVE_ICON_COLOR = RGBAf(0.09, 0.35, 0.78, 1)
const _TOOLTIP_COLOR = RGBAf(0.17, 0.18, 0.21, 0.96)
const _TOOLTIP_TEXT_COLOR = RGBAf(0.97, 0.97, 0.97, 1)
const _TRANSPARENT = RGBAf(0, 0, 0, 0)
# GLMakie draws the plots in the order of their z translation (clip range ±10000 of the pixel
# camera): the tooltip comes after the 3D view and the progress window (`_PROGRESS_Z`)
const _TOOLTIP_Z = 9000.0f0

"""
    _IconButton(parent; icon::Union{Symbol, BezierPath}, tooltip::String = "", size = 28, kwargs...)

A flat icon button of the app layout of the live view at the grid position `parent`, e.g.
`fig[1, 2]`: the icon [`_icon(icon)`](@ref _icon) on a transparent square of `size` pixels, which
gets a rounded background while the mouse is over it. A left click increments `clicks`, like
`Makie.Button`. The tooltip appears below the button after the mouse rests on it for
`tooltip_delay` seconds.

# Keyword arguments

- `icon_size = round(0.72 * size)`: size of the icon in pixels
- `icon_color`, `hover_color`, `active_color`, `active_icon_color`: colors of the theme, see
  `_ICON_COLOR` etc.; `active_*` are used by [`_IconToggle`](@ref) only
- `tooltip_color`, `tooltip_text_color`, `tooltip_placement = :below` (`:above`, `:left`,
  `:right`), `tooltip_delay = 0.5`
- further keyword arguments, e.g. `halign` or `tellwidth`, go to the `Makie.Box`

# Fields

- `box`: the `Makie.Box` that holds the layout and draws the background
- `clicks`: number of clicks
- `hovered`: whether the mouse is over the button
- `tooltip`: the text of the tooltip, may be changed
- `icon`, `icon_color`, `background`: the marker and the colors as drawn
- `plots`: background, icon and tooltip

Hover and clicks are handled by `Makie.addmouseevents!` on the bounding box: observables change
only when the mouse enters or leaves the button, there is no per-frame work.
"""
struct _IconButton
    box::Makie.Box
    clicks::Observable{Int}
    hovered::Observable{Bool}
    tooltip::Observable{String}
    icon::Observable{BezierPath}
    icon_color::Observable{RGBAf}
    background::Observable{RGBAf}
    plots::Vector{Makie.AbstractPlot}
end

"""
    _IconToggle(parent; icon::Union{Symbol, BezierPath}, icon_off = icon, tooltip::String = "",
                active = false, size = 28, kwargs...)

A flat icon toggle of the app layout of the live view, see [`_IconButton`](@ref) for the look and
the keyword arguments. A left click flips `active`, which can also be set from code, like
`Makie.Toggle.active`; `active` may be passed as an `Observable{Bool}`, which the toggle then uses.
While active, the toggle shows `icon` in `active_icon_color` on an `active_color` background,
otherwise `icon_off` in `icon_color`.

# Fields

As [`_IconButton`](@ref), with `active::Observable{Bool}` in place of `clicks`.
"""
struct _IconToggle
    box::Makie.Box
    active::Observable{Bool}
    hovered::Observable{Bool}
    tooltip::Observable{String}
    icon::Observable{BezierPath}
    icon_color::Observable{RGBAf}
    background::Observable{RGBAf}
    plots::Vector{Makie.AbstractPlot}
end

function _IconButton(parent; icon::Union{Symbol, BezierPath}, tooltip::String = "", kwargs...)
    clicks = Observable(0)
    w = _icon_widget(parent, Observable(_icon(icon)), Observable(false), tooltip,
        () -> (clicks[] += 1); kwargs...)
    return _IconButton(w.box, clicks, w.hovered, w.tooltip, w.icon, w.icon_color, w.background,
        w.plots)
end

function _IconToggle(parent; icon::Union{Symbol, BezierPath}, icon_off::Union{Symbol, BezierPath} = icon,
        tooltip::String = "", active = false, kwargs...)
    active = convert(Observable{Bool}, active)
    on_icon, off_icon = _icon(icon), _icon(icon_off)
    marker = lift(a -> a ? on_icon : off_icon, active)
    w = _icon_widget(parent, marker, active, tooltip, () -> (active[] = !active[]); kwargs...)
    return _IconToggle(w.box, active, w.hovered, w.tooltip, w.icon, w.icon_color, w.background,
        w.plots)
end

Base.delete!(b::Union{_IconButton, _IconToggle}) = delete!(b.box)

"""
    _icon_widget(parent, icon, active, tooltip, onclick; kwargs...)

Builds the plots and the mouse handling shared by [`_IconButton`](@ref) and
[`_IconToggle`](@ref); `onclick()` runs on a left click. Returns a `NamedTuple` of the fields.
"""
function _icon_widget(parent, icon::Observable{BezierPath}, active::Observable{Bool},
        tooltip::String, onclick;
        size::Real = 28, icon_size::Real = round(0.72 * size),
        icon_color = _ICON_COLOR, hover_color = _ICON_HOVER_COLOR,
        active_color = _ICON_ACTIVE_COLOR, active_icon_color = _ICON_ACTIVE_ICON_COLOR,
        tooltip_color = _TOOLTIP_COLOR, tooltip_text_color = _TOOLTIP_TEXT_COLOR,
        tooltip_placement::Symbol = :below, tooltip_delay::Real = 0.5, layout_kwargs...)
    c_icon, c_hover = RGBAf(Makie.to_color(icon_color)), RGBAf(Makie.to_color(hover_color))
    c_active = RGBAf(Makie.to_color(active_color))
    c_active_icon = RGBAf(Makie.to_color(active_icon_color))
    hovered = Observable(false)
    background = Observable(_TRANSPARENT)
    fg = Observable(c_icon)
    box = Makie.Box(parent; width = size, height = size, color = background,
        strokevisible = false, cornerradius = max(size / 6, 3), layout_kwargs...)
    scene = box.blockscene
    bg = only(scene.plots)
    bbox = box.layoutobservables.computedbbox
    # Only on changes of the state, not per mouse move
    function update_look!(_...)
        b = active[] ? c_active : hovered[] ? c_hover : _TRANSPARENT
        f = active[] ? c_active_icon : c_icon
        background[] == b || (background[] = b)
        fg[] == f || (fg[] = f)
        return nothing
    end
    onany(update_look!, scene, hovered, active)
    update_look!()
    center = lift(r -> Point2f(Makie.origin(r) .+ Makie.widths(r) ./ 2), scene, bbox)
    marker = scatter!(scene, center; marker = icon, markersize = icon_size, color = fg,
        markerspace = :pixel, visible = box.visible, inspectable = false)
    translate!(marker, 0, 0, 1)

    # Tooltip, placed when shown
    tip_text = Observable(tooltip)
    tip_visible = Observable(false)
    tip_pos = Observable(Point2f(0))
    tip_align = Observable(0.5f0)
    tip = Makie.tooltip!(scene, tip_pos, tip_text; placement = tooltip_placement,
        align = tip_align, visible = tip_visible, backgroundcolor = tooltip_color,
        textcolor = tooltip_text_color, outline_linewidth = 0, triangle_size = 6, offset = 4,
        fontsize = 13, textpadding = (6, 6, 4, 4), overdraw = true, inspectable = false)
    translate!(tip, 0, 0, _TOOLTIP_Z)
    function show_tooltip!()
        (hovered[] && !isempty(tip_text[])) || return nothing
        r = bbox[]
        (x0, y0), (w, h) = Makie.origin(r), Makie.widths(r)
        tip_pos[] = tooltip_placement === :below ? Point2f(x0 + w / 2, y0) :
                    tooltip_placement === :above ? Point2f(x0 + w / 2, y0 + h) :
                    tooltip_placement === :left ? Point2f(x0, y0 + h / 2) :
                    Point2f(x0 + w, y0 + h / 2)
        # Keep the tooltip of buttons near the left or right edge inside the window: `align` is
        # the position of the arrow along the tooltip, from the left
        if tooltip_placement in (:below, :above)
            W = Makie.widths(Makie.viewport(scene)[])[1]
            tip_align[] = x0 + w / 2 < W / 4 ? 0.15f0 : x0 + w / 2 > 3W / 4 ? 0.85f0 : 0.5f0
        end
        tip_visible[] = true
        return nothing
    end
    timer = Ref{Union{Nothing, Timer}}(nothing)
    function hide_tooltip!()
        isnothing(timer[]) || (close(timer[]); timer[] = nothing)
        tip_visible[] && (tip_visible[] = false)
        return nothing
    end

    mouse = Makie.addmouseevents!(scene, bbox)
    on(scene, mouse.obs) do event
        t = event.type
        if t === Makie.MouseEventTypes.enter && box.visible[]
            hovered[] = true
            if tooltip_delay > 0
                timer[] = Timer(_ -> show_tooltip!(), tooltip_delay)
            else
                show_tooltip!()
            end
        elseif t === Makie.MouseEventTypes.out
            hovered[] && (hovered[] = false)
            hide_tooltip!()
        elseif t in (Makie.MouseEventTypes.leftclick, Makie.MouseEventTypes.leftdoubleclick) &&
               box.visible[]
            # The second of two fast clicks is a double click, which counts as a click
            hide_tooltip!()
            onclick()
            return Consume(true)
        end
        return Consume(false)
    end
    plots = Makie.AbstractPlot[bg, marker, tip]
    return (; box, hovered, tooltip = tip_text, icon, icon_color = fg, background, plots)
end
