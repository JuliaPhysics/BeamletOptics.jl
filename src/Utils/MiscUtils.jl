"""
    list_subtypes(T::Type; max_depth::Int=5)

Prints a tree of all subtypes, e.g. `list_subtypes(AbstractObject)`.
Maximum exploration depth can be limited by passing `max_depth`.
Returns the total number of types encountered.
"""
function list_subtypes(parent::Type, is_last=true, prefix="", depth=1; max_depth=10)
    # Determine children of parent type
    children = subtypes(parent)
    num_children = length(children)
    # Print out parent type information
    connector = is_last ? "└── " : "├── "
    msg = prefix * connector * string(parent)
    println(msg)
    # Return to parent function if childless
    if iszero(num_children)
        return 0
    end
    # Extend prefix string based on layer depth
    if is_last
        prefix = prefix * "    "
    else
        prefix = prefix * "│   "
    end
    # Check if max tree depth has been reached
    if depth > max_depth
        msg = prefix * "└── " * "$num_children more ..."
        println(msg)
        return num_children
    end
    # Delve deeper into tree
    for child in children[1:end-1]
        num_children += list_subtypes(child, false, prefix, depth+1; max_depth)
    end
    num_children += list_subtypes(last(children), true, prefix, depth+1; max_depth)
    # If tree/branch has been resolved, return total number of types found.
    if depth == 1
        msg = "\n" * "At least $num_children types have been found."
        println(msg)
    end
    return num_children
end

"""
    countlines_in_dir(dir)

Counts the number of lines of all `.jl` files in `dir`.
"""
function countlines_in_dir(dir::String)
    l = 0
    items = readdir(dir)
    for item in items
        path = joinpath(dir, item)
        if isfile(path) && endswith(path, ".jl")
            li = countlines(path)
            l += li
            println("$li lines of code in $path")
        end
        if isdir(path)
            l += countlines_in_dir(path)
        end
    end
    return l
end