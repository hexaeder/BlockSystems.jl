"""
    _uncouple_algebraic_equations(eqs; verbose=false)

Uncouple a given set of explicit algebraic equations. Returns a tuple

    (rules, keep)

of independent substitution rules (non recursive) and a list of indices
of equations, which could not be reduced.
"""
function _uncouple_algebraic_equations(eqs; verbose=false)
    @check all([eq_type(eq)[1] for eq in eqs] .== :explicit_algebraic) "Can only hanlde expl. alg. eqs!"

    eqs = deepcopy(eqs)
    g = _expl_algebraic_dependency_graph(eqs)

    keep = Int[]
    for set in cyclebreaking_vertices(g)
        if length(set) == 1
            push!(keep, only(set))
        else
            # TODO: cleverly pick which equations to keep
            # there are several possible combinations of equations to keep which all break the
            # circles. Maybe some are better then other given that we later try to solve the
            # new implicit algebraic equations?
           push!(keep, last(set))
        end
    end

    _substitute_along_depgraph!(eqs, g, keep)
    rules = Dict(eqs[i].lhs => eqs[i].rhs for i in eachindex(eqs) if i ∉ keep)

    # if the algorightm could not remove alle of the equations, there have been
    verbose && !isempty(keep) && println("Try to handle the newly created implicit equations")
    transformed = map!(eq->_transform_implicit_algebraic(eq; verbose), eqs[keep], eqs[keep])
    explicit = findall(isequal(:explicit_algebraic), map(eq -> eq_type(eq)[1], transformed))
    if !isempty(explicit)
        verbose && println("This lead to new explicit algebraic equations! Recursivly uncouple")
        subrules, subkeep = _uncouple_algebraic_equations(transformed[explicit]; verbose)
        keep = keep[subkeep]

        for (k, r) in rules
            rules[k] = substitute(r, subrules)
        end
        merge!(rules, subrules)
    end

    return (rules, keep)
end


"""
    _substitute_along_depgraph!(eqs, g, keep)

Remove all connections from from equations in `keep`. The remaining graph
should be cycle free. Start with the outermost equations and substitute along
the dependency graph.
After this, the `eqs` object should only depend on `eqs[keep]` without any further
interdependency between them!
"""
function _substitute_along_depgraph!(eqs, g::SimpleDiGraph, keep::Vector{Int})
    gcut = _cut_substitutions_from(g, keep)

    while ne(gcut) != 0
        for i in 1:nv(gcut)
            # only substitute if there are no dependencies
            !iszero(indegree(gcut, i)) && continue

            rule = eqs[i].lhs => eqs[i].rhs
            for nb in copy(outneighbors(gcut, i))
                eq = eqs[nb]
                eqs[nb] = eq.lhs ~ substitute(eq.rhs, rule)
                t = rem_edge!(gcut, i, nb)
                @assert t "Could not remove edge. Weird."
            end
        end
    end
end

"""
    _cut_substitutions_from(g, vertices)

Cut all connections from and to given vertices in Graph. Returns a new graph.
"""
function _cut_substitutions_from(g, vertices)
    g = copy(g)
    for v in vertices
        for nb in copy(outneighbors(g, v))
            rem_edge!(g, v, nb)
        end
    end
    return g
end

"""
    _expl_algebraic_dependency_graph(eqs)

Get the dependency graph for a given set of explicit algebraic equations.
Edge 1->2 means, that state 1 appears in state 2.
"""
function _expl_algebraic_dependency_graph(eqs)
    g = SimpleDiGraph(length(eqs))
    symbols = [eq.lhs for eq ∈ eqs]
    for (i, eq) in enumerate(eqs)
        rhs_vars = get_variables(eq.rhs)
        for (isym, sym) in enumerate(symbols)
            if Set([sym]) ⊆ Set(rhs_vars)
                add_edge!(g, isym => i)
            end
        end
    end
    return g
end

# This code ist adapted from the Hawick & James implementation in Graphs.jl which is
# licenced as follows: https://github.com/JuliaGraphs/Graphs.jl/blob/master/LICENSE.md
# Copyright (c) 2015: Seth Bromberger and other contributors. Copyright (c) 2012: John Myles White and other contributors.

"""
    cyclebreaking_vertices(g)

Find a minimal set of vertices, which has an intersection with all circles in the Graph,
i. e. by removing those vertices, the Graph becomes cycle free.

Returns a `Vector{Vector{T}}` for a given `SimpleDiGraph{T}`.
Each combination, which picks at least 1 element from each of those subvectors will be sufficient
to break all cycles.

### References
- Hawick & James, "Enumerating Circuits and Loops in Graphs with Self-Arcs and Multiple-Arcs", 2008
- Julia Implementation taken from https://github.com/JuliaGraphs/Graphs.jl/blob/master/src/cycles/hawick-james.jl
"""
function cyclebreaking_vertices(g::SimpleDiGraph{T}) where T
    nvg = nv(g)
    B = Vector{T}[Vector{T}() for i in vertices(g)]
    blocked = zeros(Bool, nvg)
    stack = Vector{T}()
    cycle_intersecs = Vector{Set{T}}()
    for v in vertices(g)
        cycle_intersec_recursive!(g, v, v, blocked, B, stack, cycle_intersecs)

        # reset blocked and B
        fill!(blocked, false)
        map!(empty!, B, B)
    end
    return [collect(set) for set in cycle_intersecs]
end

"""
    cycle_intersec_recursive!(g, v1, v2, blocked, B, stack, cycle_intersecs)

Find circuits in `g` recursively starting from v1.
"""
function cycle_intersec_recursive!(g::SimpleDiGraph, v1, v2, blocked, B, stack, cycle_intersecs)
    f = false
    push!(stack, v2)
    blocked[v2] = true

    # if the stack is a superset of on of those cycle intersections
    # treat it as if a cycle was just found. This will prevent the algorithm from
    # further sucht which would not help top narrow down the cycle intersections.
    for c in cycle_intersecs
        if stack ⊇ c
            pop!(stack)
            unblock!(v2, blocked, B)
            return true
        end
    end

    Av = outneighbors(g, v2)
    for w in Av
        (w < v1) && continue
        if w == v1 # Found a circuit
            # first check if this circuit has a cut with an existing cycle_intersec
            matching_intersection = false
            for i in eachindex(cycle_intersecs)
                cut = cycle_intersecs[i] ∩ stack
                if !isempty(cut)
                    cycle_intersecs[i] = cut
                    matching_intersection = true
                    break
                end
            end
            # if the cycle contains only new vertices, create new set
            if !matching_intersection
                push!(cycle_intersecs, Set(stack))
            end
            f = true
        elseif !blocked[w]
            f = cycle_intersec_recursive!(g, v1, w, blocked, B, stack, cycle_intersecs)
        end
    end
    if f
        unblock!(v2, blocked, B)
    else
        for w in Av
            (w < v1) && continue
            if !(v2 in B[w])
                push!(B[w], v2)
            end
        end
    end
    pop!(stack)
    return f
end

"""
    unblock!(v, blocked, B)

Unblock the value `v` from the `blocked` list and remove from `B`.
"""
function unblock!(v, blocked, B)
    blocked[v] = false
    Bv = B[v]

    while !isempty(Bv)
        w = pop!(Bv)
        if blocked[w]
            unblock!(w, blocked, B)
        end
    end

    return nothing
end
