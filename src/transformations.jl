export connect_system

"""
$(SIGNATURES)

Recursively transform `IOSystems` to `IOBlocks`.

- substitute inputs with connected outputs
- try to eliminate equations for internal states which are not used to calculate the specified outputs of the system.
- try to eliminate explicit algebraic equations (i.e. outputs of internal blocks) by substituting each occurrence
  with their rhs. Explicit algebraic states which are marked as system outputs won't be removed.

Parameters:
- `ios`: system to connect
- `verbose=false`: toggle verbosity (show equations at different steps)
- `simplify_eqs=true`: toggle simplification of all equations at the end
"""
function connect_system(ios::IOSystem; verbose=false, simplify_eqs=true)
    # recursive connect all subsystems
    for (i, subsys) in enumerate(ios.systems)
        if subsys isa IOSystem
            ios.systems[i] = connect_system(subsys, verbose=verbose)
        end
    end
    eqs = vcat([namespace_equations(iob.system) for iob in ios.systems]...)

    verbose && @info "Transform IOSystem $(ios.name) to IOBlock" ios.name ios.inputs ios.outputs ios.connections eqs

    # get rid of closed inputs by substituting output states
    connections = ios.connections
    for (i, eq) in enumerate(eqs)
        eqs[i] = eq.lhs ~ substitute(eq.rhs, connections)
    end

    verbose && @info "substitute inputs with outputs" eqs

    # get red of unused states
    # (internal variables which are not used for the outputs)
    reduced_eqs1 = remove_superfluous_states(eqs, keys(ios.outputs_map))
    verbose && @info "without superfluous states" reduced_eqs1

    # reduce algebraic states of the system
    (reduced_eqs2, removed_eqs) = remove_algebraic_states(reduced_eqs1, skip = keys(ios.outputs_map))
    verbose && @info "without explicit algebraic states" reduced_eqs2

    # apply the namespace transformations
    namespace_promotion = merge(ios.inputs_map, ios.iparams_map, ios.istates_map, ios.outputs_map)

    promoted_eqs = deepcopy(reduced_eqs2)
    for (i, eq) in enumerate(promoted_eqs)
        promoted_eqs[i] = eqsubstitute(eq, namespace_promotion)
    end
    for (i, eq) in enumerate(removed_eqs)
        removed_eqs[i] = eqsubstitute(eq, namespace_promotion)
    end
    verbose && @info "promoted namespaces" promoted_eqs removed_eqs

    if simplify_eqs
        promoted_eqs = simplify.(promoted_eqs)
        removed_eqs = simplify.(removed_eqs)
        verbose && @info "simplified equations" promoted_eqs removed_eqs
    end

    try
        IOBlock(promoted_eqs, ios.inputs, ios.outputs, name=ios.name, removed_eq=removed_eqs)
    catch e
        @error "Failed to build IOBlock from System" ios.inputs ios.outputs ios.name eqs reduced_eqs1 reduced_eqs2 promoted_eqs removed_eqs
        throw(e)
    end
end

"""
    remove_superfluous_states(eqs::Vector{Equation}, outputs)

This function removes equations, which are not used in order to generate the given
`outputs`. It recursively looks for attracting components in the dependency graph which
are not connected to the output nodes.
Returns a new, reduced set of equations without these states.

TODO: Change strategy, remove i if there is no path from i to outputs.
"""
function remove_superfluous_states(eqs::Vector{Equation}, outputs)
    neweqs = deepcopy(eqs)
    sys = ODESystem(neweqs) # will be used for the dependency graph
    neweqs = sys.eqs # the ODESystem might reorder the equations
    # generate dependency graph
    graph = eqeq_dependencies(asgraph(sys), variable_dependencies(sys))
    # find 'main' eq for each output
    eq_idx = [findfirst(x->o ∈ Set(ModelingToolkit.vars(x.lhs)), neweqs) for o in outputs]
    # find the attracting components, remove attracting components which do not include eqs
    attr_components = attracting_components(graph)
    removable = []
    for attr in attr_components
        if isempty(Set(eq_idx) ∩ Set(attr))
            append!(removable, attr)
        end
    end
    @assert allunique(removable)
    deleteat!(neweqs, sort(removable))

    # reduction has to be repeated recursively
    if neweqs == eqs
        return neweqs
    else
        return remove_superfluous_states(neweqs, outputs)
    end
end

"""
    remove_algebraic_states(eqs:Vector{Equation}, skip=[])

Reduces the number of equations by substituting explicit algebraic equations.
Returns a tuple containing two `Vector{Equations}`
  - the new equations with substituted states
  - the algebraic states which have been removed

The optional skip argument is used to declare states which should not be removed.

This function checks for cyclic dependencies between algebraic equations by generating
a dependency graph between them. All substitutions have to be pairwise cycle free.
```example
julia> @variables i x y o1 o2;
julia> @derivatives D'~t;
julia> eqs = [D(x) ~ i,
              o1 ~ x + o2,
              D(y) ~ i,
              o2 ~ y + o1];
julia> remove_algebraic_states(eqs)
(Equation[Equation((D'~t)(x), i),
          Equation(o1, x + (o1 + y)),
          Equation((D'~t)(y), i)],
 Equation[Equation(o2, o1 + y)])
```
"""
function remove_algebraic_states(eqs::Vector{Equation}; skip=[])
    reduced_eqs = deepcopy(eqs)

    # only consider states for reduction which are explicit algebraic and not in skip
    condition = eq -> is_explicit_algebraic(eq) && !(Set(get_variables(eq.lhs)) ⊆ Set(skip))
    algebraic_idx = findall(condition, reduced_eqs)

    # symbols of all algebraic eqs
    symbols = [eq.lhs for eq ∈ reduced_eqs[algebraic_idx]]

    # generate dependency graph
    g = SimpleDiGraph(length(algebraic_idx))
    for (i, eq) in enumerate(reduced_eqs[algebraic_idx])
        rhs_vars = vars(eq.rhs)
        for (isym, sym) in enumerate(symbols)
            if Set([sym]) ⊆ Set(rhs_vars)
                add_edge!(g, isym => i)
            end
        end
    end
    removable = algebraic_idx[pairwise_cycle_free(g)]
    @assert allunique(removable)

    rules = Dict(eq.lhs => eq.rhs for eq in reduced_eqs[removable])
    # subsitute all the equations, remove substituted
    for (i, eq) in enumerate(reduced_eqs)
        reduced_eqs[i] = eq.lhs ~ substitute(eq.rhs, rules)
    end

    removed_eqs = reduced_eqs[removable]
    deleteat!(reduced_eqs, sort(removable))

    return reduced_eqs, removed_eqs
end

"""
    pairwise_cycle_free(g:SimpleDiGraph)

Returns an array of vertices, which pairwise do not belong to any cycle in `g`.
Uses `simplecycles` from `LightGraphs`. The algorithm starts with all vertices
and iteratively removes the vertices is part of most cycles.
"""
function pairwise_cycle_free(g::SimpleDiGraph)
    idx = collect(vertices(g))
    cycles = simplecycles(g)
    while true
        cycles_per_idx = zeros(Int, length(idx))
        for (i, id1) ∈ enumerate(idx)
            for (j, id2) ∈ enumerate(@view idx[i+1 : end])
                for c in cycles
                    if id1 ∈ c && id2 ∈ c
                        cycles_per_idx[i] += 1
                        cycles_per_idx[i+j] += 1
                    end
                end
            end
        end
        if sum(cycles_per_idx) > 0
            worst_idx = findmax(cycles_per_idx)[2]
            deleteat!(idx, worst_idx)
        else
            return idx
        end
    end
end
