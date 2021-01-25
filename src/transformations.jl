export connect_system

"""
$(SIGNATURES)

Recursively transform `IOSystems` to `IOBlocks`.

- substitute inputs with connected outputs
- try to eliminate equations for internal states which are not used to calculate the specified outputs of the system.
- try to eliminate explicit algebraic equations (i.e. outputs of internal blocks) by substituting each occurrence
  with their rhs. Explicit algebraic states which are marked as system outputs won't be reduced.

Parameters:
- `ios`: system to connect
- `verbose=false`: toggle verbosity (show equations at different steps)
"""
function connect_system(ios::IOSystem; verbose=false)
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

    verbose && @info "subsitute inputs with outputs" eqs

    # TODO: check assumption
    # - each state is represented by one lhs
    # - lhs only first order or algebraic
    # - no self loop in algebraic (not needed?)
    # - check that every variable found in the equations is referenced by namespace map
    # - dependency graph does not work correctly for implicit algebraic states! forbid them
    #   -> does this have implications for the reduction of depended explicit a states?

    # get red of unused states
    # (internal variables which are not used for the outputs)
    reduced_eqs1 = reduce_superflous_states(eqs, keys(ios.outputs_map))
    verbose && @info "w/o superflous states" reduced_eqs1

    # reduce algebraic states of the system
    reduced_eqs2 = reduce_algebraic_states(reduced_eqs1, skip = keys(ios.outputs_map))
    verbose && @info "with reduced algebraic states" reduced_eqs2

    # TODO redo namespace promotion after deletion of states?
    # apply the namespace transformations
    namespace_promotion = merge(ios.inputs_map, ios.iparams_map, ios.istates_map, ios.outputs_map)

    promoted_eqs = deepcopy(reduced_eqs2)
    for (i, eq) in enumerate(promoted_eqs)
        promoted_eqs[i] = eqsubstitute(eq, namespace_promotion)
    end
    verbose && @info "promoted namespaces" promoted_eqs

    try
        IOBlock(promoted_eqs, ios.inputs, ios.outputs, name=ios.name)
    catch e
        @error "Failed to build IOBlock from System" ios.inputs ios.outputs ios.name eqs reduced_eqs1 reduced_eqs2 promoted_eqs
        throw(e)
    end
end

"""
    reduce_superflous_states(eqs::Vector{Equation}, outputs)

This function reduce equations which are not used in order to generate the given
`outputs`. It recoursively looks for attractic components in the dependency graph which
are not connected to the output nodes.
Returns a new set of equations without these states.

TODO: Change strategy, remove i if there is no path from i to outputs.
"""
function reduce_superflous_states(eqs::Vector{Equation}, outputs)
    neweqs = deepcopy(eqs)
    sys = ODESystem(neweqs) # will be used for the dependency graph
    neweqs = sys.eqs # the ODESystem might reorder the equations
    # generate dependency graph
    graph = eqeq_dependencies(asgraph(sys), variable_dependencies(sys))
    # find 'main' eq for each output
    eq_idx = [findfirst(x->o ∈ Set(ModelingToolkit.vars(x.lhs)), neweqs) for o in outputs]
    # find the attracting components, remove attracting components which do not include eqs
    attr_components = attracting_components(graph)
    deleq = []
    for attr in attr_components
        if isempty(Set(eq_idx) ∩ Set(attr))
            append!(deleq, attr)
        end
    end
    @assert isunique(deleq)
    deleteat!(neweqs, sort(deleq))

    # reduction has to be repeated recursively
    if neweqs == eqs
        return neweqs
    else
        return reduce_superflous_states(neweqs, outputs)
    end
end

"""
    reduce_algebraic_states(eqs:Vector{Equation}, skip=[])

Reduces the number of equations by substituting explicit algebraic equations.
Returns a new set of equations.

The optional skip argument is used to declare states which should not be reduced.

This function checks for cyclic dependencies between algebraic equations by generating
a dependency graph between them. All substitutions have to be pairwise cycle free.
```example
julia> @variables i x y o1 o2;
julia> @derivatives D'~t;
julia> eqs = [D(x) ~ i,
              o1 ~ x + o2,
              D(y) ~ i,
              o2 ~ y + o1];
julia> reduce_algebraic_states(eqs)
3-element Array{Equation,1}:
 Equation(Differential(x), i)
 Equation(o1, x + (y + o1))
 Equation(Differential(y), i)
```
"""
function reduce_algebraic_states(eqs::Vector{Equation}; skip=[])
    neweqs = deepcopy(eqs)

    # only consider states for reduction which are explicit algebraic and not in skip
    condition = eq -> is_explicit_algebraic(eq) && !(Set(get_variables(eq.lhs)) ⊆ Set(skip))
    algebraic_idx = findall(condition, neweqs)

    # symbols of all algebraic eqs
    symbols = [eq.lhs for eq ∈ neweqs[algebraic_idx]]

    # generate dependency graph
    g = SimpleDiGraph(length(algebraic_idx))
    for (i, eq) in enumerate(neweqs[algebraic_idx])
        rhs_vars = vars(eq.rhs)
        for (isym, sym) in enumerate(symbols)
            if Set([sym]) ⊆ Set(rhs_vars)
                add_edge!(g, isym => i)
            end
        end
    end
    reducable = algebraic_idx[pairwise_cycle_free(g)]

    rules = Dict(eq.lhs => eq.rhs for eq in neweqs[reducable])
    # subsitute all the equations, remove substituted
    for (i, eq) in enumerate(neweqs)
        neweqs[i] = eq.lhs ~ substitute(eq.rhs, rules)
    end
    @assert isunique(reducable)
    deleteat!(neweqs, sort(reducable))

    return neweqs
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
