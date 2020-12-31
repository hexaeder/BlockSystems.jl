export connect_system

function connect_system(ios::IOSystem)
    # recursive connect all subsystems
    for (i, subsys) in enumerate(ios.systems)
        if subsys isa IOSystem
            ios.systems[i] = connect_system(subsys)
        end
    end
    eqs = vcat([namespace_equations(iob.system) for iob in ios.systems]...)

    # get rid of closed inputs by substituting output states
    connections = ios.connections
    for (i, eq) in enumerate(eqs)
        eqs[i] = eq.lhs ~ substitute(eq.rhs, connections)
    end

    # TODO: check assumtions
    # - each state is represented by one lhs
    # - lhs only first order or algebraic
    # - no self loop in algeraic (not needed?)
    # - check that every variable found in the equations is referenced by namespacemap
    # - dependency graph does not work correctly for implicit algebraic states! forbid them
    #   -> does this have implications for the reduction of depended explicit a states?

    # get red of unused states
    # (internal variables which are not used for the outputs)
    eqs = reduce_superflous_states(eqs, keys(ios.outputs_map))

    # reduce algebraic states of the system
    eqs = reduce_algebraic_states(eqs, skip = keys(ios.outputs_map))

    # apply the namespace transformations
    namespace_promotion = merge(ios.inputs_map, ios.iparams_map, ios.istates_map, ios.outputs_map)

    for (i, eq) in enumerate(eqs)
        eqs[i] = eqsubstitute(eq, namespace_promotion)
    end

    # eqs
    IOBlock(eqs, ios.inputs, ios.outputs, name=ios.name)
end

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
    deleteat!(neweqs, deleq)

    # reduction has to be repeated recursively
    if neweqs == eqs
        return neweqs
    else
        return reduce_superflous_states(neweqs, outputs)
    end
end

function reduce_algebraic_states(eqs::Vector{Equation}; skip=[])
    eqs = deepcopy(eqs)
    sys = ODESystem(eqs) # will be used for the dependency graph
    eqs = sys.eqs # the ODESystem might reorder the equations

    # only consider states for reduction which are explicit algebraic and not in skip
    condition = eq -> is_explicit_algebraic(eq) && !(Set(get_variables(eq.lhs)) ⊆ Set(skip))
    algebraic_idx = findall(condition, eqs)

    # generate dependency graph
    graph = eqeq_dependencies(asgraph(sys), variable_dependencies(sys))

    # algebraic states which do no have cyclic dependencies between them can be reduced
    cycles = simplecycles(graph)
    reducable = pairwise_cycle_free(algebraic_idx, cycles)
    rules = [eq.lhs => eq.rhs for eq in eqs[reducable]]
    # subsitute all the equations, remove substituted
    for (i, eq) in enumerate(eqs)
        eqs[i] = eq.lhs ~ substitute(eq.rhs, rules)
    end
    deleteat!(eqs, reducable)
    return eqs
end

function pairwise_cycle_free(idx, cycles)
    pairs = collect(Iterators.product(idx, idx))
    # check for each pair whether found in cycle
    incycle(p) = p[1]==p[2] ? false : any(map(c -> p ⊆ c, cycles))
    A = incycle.(pairs)
    if !any(A)
        return idx
    else
        # find the variable with the highes number of loops
        cycle_count = reduce(+, A, dims=1)[:]
        worst_idx = findmax(cycle_count)[2]
        deleteat!(idx, worst_idx)
        return pairwise_cycle_free(idx, cycles)
    end
end
