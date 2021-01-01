export generate_io_function

"""
$(SIGNATURES)

Generate calleble functions for an `AbstractIOSystem`. An `IOSystem` will be transformed to an `IOBlock` first.

Parameters:
- `ios`: the system to build the function
- `first_states`: define states=(istates ∪ outputs) which shoud appar first
- `first_inputs`: define inputs which shoud appar first
- `simplify=true`: call simplify on the equations before building?
- `expression=Val{false}`: toggle expression and callable function output

Returns an named tuple with the fields
- `f_ip` in-place function `f(dstates, states, inputs, params, iv)`
- `f_oop` out-of-place function `f(states, inputs, params, iv)=>dstates`
- `massm` mass matrix of the system
- `states` symbols of states (in order)
- `inputs` symbols of inputs (in order)
- `params` symbols of parameters (in order)

"""
function generate_io_function(ios::AbstractIOSystem; first_states = [], first_inputs = [], simplify=true, expression = Val{false})
    if ios isa IOSystem
        ios = connect_system(ios)
    end
    # first_outputs and first_inputs may be given in namepsace version
    first_states = remove_namespace.(ios.name, value.(first_states))
    first_inputs = remove_namespace.(ios.name, value.(first_inputs))
    Set(first_states) ⊆ (Set(ios.outputs) ∪ Set(ios.istates)) || throw(ArgumentError("first_states !⊆ (outputs ∪ istates)"))
    Set(first_inputs) ⊆ Set(ios.inputs) || throw(ArgumentError("first_inputs !⊆ inputs"))

    # enforce ordering of states and inputs
    states = vcat(first_states, ios.outputs, ios.istates) |> unique
    inputs = vcat(first_inputs, ios.inputs) |> unique
    params = ios.iparams

    # equations of form o = f(...) have to be transfomed to 0 = f(...) - o

    eqs = transform_algebraic_equations(ios.system.eqs)

    # simplify equations if wanted
    if simplify
        eqs = ModelingToolkit.simplify.(eqs)
    end

    # reorder the equations to get du in the right order
    eqs = reorder_by_states(eqs, states)

    # create massmatrix, we don't use the method provided by ODESystem because of reordering
    mass_matrix = generate_massmatrix(eqs)

    # substitute x(t) by x for all terms
    states′ = makesym.(states, states=[])
    inputs′ = makesym.(inputs, states=[])
    params′ = makesym.(params, states=[])

    sub = merge(Dict(states .=> states′),
                Dict(inputs .=> inputs′),
                Dict(params .=> params′))
    formulas = [substitute(eq.rhs, sub) for eq in eqs]

    # generate functions
    f_oop, f_ip = build_function(formulas, states′, inputs′, params′, ios.system.iv; expression = expression)

    return (f_oop=f_oop, f_ip=f_ip, massm=mass_matrix, states=states′, inputs=inputs′, params=params′)
end

function transform_algebraic_equations(eqs::AbstractVector{Equation})
    eqs = deepcopy(eqs)
    for (i, eq) in enumerate(eqs)
        if isequal(eq.lhs, 0) || (eq.lhs isa Term && eq.lhs.op isa Differential)
            continue
        else
            eqs[i] = 0 ~ eq.rhs - eq.lhs
        end
    end
    return eqs
end

function reorder_by_states(eqs::AbstractVector{Equation}, states)
    @assert length(eqs) == length(states) "Numbers of eqs should be equal to states!"

    eq_idx = [findfirst(x->s ∈ Set(ModelingToolkit.vars(x.lhs)), eqs) for s in states]

    unused_idx = reverse(setdiff(1:length(eqs), eq_idx))
    for i in 1:length(eq_idx)
        if eq_idx[i] === nothing
            eq_idx[i] = pop!(unused_idx)
        end
    end
    @assert sort(unique(eq_idx)) == 1:length(eqs) "eq_idx shoud contain all idx!"
    return eqs[eq_idx]
end

function generate_massmatrix(eqs::AbstractVector{Equation})
    V = Vector{Int}(undef, length(eqs))
    for i in 1:length(eqs)
        if eqs[i].lhs isa Term && eqs[i].lhs.op isa Differential
            V[i] = 1
        elseif isequal(eqs[i].lhs, 0)
            V[i] = 0
        else
            error("Cant build mass matrix entry for $(eqs[i])")
        end
    end
    M = Diagonal(V)
    return M==I ? I : M
end
