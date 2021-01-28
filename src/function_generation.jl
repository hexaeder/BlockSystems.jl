export generate_io_function

"""
$(SIGNATURES)

Generate callable functions for an `AbstractIOSystem`. An `IOSystem` will be transformed to an `IOBlock` first.
At this level there is no more distinction between internal states and outputs:
`states=(istates ∪ outputs)`.

Parameters:
- `ios`: the system to build the function
optional:
- `type=:auto`: `:ode` or `:static`, determines the output of the function
- `first_states`: define states=(istates ∪ outputs) which should appear first
- `first_inputs`: define inputs which should appear first
- `first_params`: define parameters which should appear first
- `expression=Val{false}`: toggle expression and callable function output

Returns an named tuple with the fields
- for `type=:ode`:
  - `f_ip` in-place function `f(dstates, states, inputs, params, iv)`
  - `f_oop` out-of-place function `f(states, inputs, params, iv) => dstates`
  - `massm` mass matrix of the system
- for `type=:static`:
  - `f_ip` in-place function `f(states, inputs, params, iv)`
  - `f_oop` out-of-place function `f(inputs, params, iv) => states`
- always:
- `states` symbols of states (in order)
- `inputs` symbols of inputs (in order)
- `params` symbols of parameters (in order)

"""
function generate_io_function(ios::AbstractIOSystem; first_states = [], first_inputs = [],
                              first_params = [], first_removed = [],
                              expression = Val{false}, verbose=false, type=:auto)
    if ios isa IOSystem
        @info "Transform given system $(ios.name) to block"
        ios = connect_system(ios, verbose=verbose)
    end
    # first_outputs, first_inputs and first_params may be given in namepsace version
    first_states = remove_namespace.(ios.name, value.(first_states))
    first_inputs = remove_namespace.(ios.name, value.(first_inputs))
    first_params = remove_namespace.(ios.name, value.(first_params))

    Set(first_states) ⊆ (Set(ios.outputs) ∪ Set(ios.istates)) || throw(ArgumentError("first_states !⊆ (outputs ∪ istates)"))
    Set(first_inputs) ⊆ Set(ios.inputs) || throw(ArgumentError("first_inputs !⊆ inputs"))
    Set(first_params) ⊆ Set(ios.iparams) || throw(ArgumentError("first_params !⊆ iparams"))

    # enforce ordering of states, inputs and params
    states = vcat(first_states, ios.outputs, ios.istates) |> unique
    inputs = vcat(first_inputs, ios.inputs) |> unique
    params = vcat(first_params, ios.iparams) |> unique

    # warning
    if length(first_states) != length(states) || length(first_inputs) != length(inputs) || length(first_params) != length(params)
        @warn "The ordering of the states/inputs/params might change from run to run. Therefore it is highly recommend to provide all variables in the first_* arguments" states inputs params
    end

    # reorder the equations to get du in the right order
    eqs = reorder_by_states(ios.system.eqs, states)
    verbose && @info "Reordered eqs" eqs states

    if type == :auto
        type = all_static(eqs) ? :static : :ode
        verbose && @info "auto-equation type: $type"
    end

    local mass_matrix
    if type == :ode
        # equations of form o = f(...) have to be transformed to 0 = f(...) - o
        eqs = transform_algebraic_equations(eqs)
        verbose && @info "Transformed algebraic eqs" eqs

        # create massmatrix, we don't use the method provided by ODESystem because of reordering
        mass_matrix = generate_massmatrix(eqs)
        verbose && @info "Reordered by states and generated mass matrix" mass_matrix
    elseif type == :static
        all_static(eqs) || throw(ArgumentError("Equations of system are not static!"))
    else
        throw(ArgumentError("Unknown type $type"))
    end

    # substitute x(t) by x for all terms
    state_syms = makesym.(states, states=[])
    input_syms = makesym.(inputs, states=[])
    param_syms = makesym.(params, states=[])

    sub = merge(Dict(states .=> state_syms),
                Dict(inputs .=> input_syms),
                Dict(params .=> param_syms))
    formulas = [substitute(eq.rhs, sub) for eq in eqs]

    # generate functions
    if type == :ode
        f_oop, f_ip = build_function(formulas, state_syms, input_syms, param_syms, ios.system.iv; expression = expression)
        return (f_oop=f_oop, f_ip=f_ip,
                massm=mass_matrix,
                states=state_syms,
                inputs=input_syms,
                params=param_syms)
    elseif type == :static
        f_oop, f_ip = build_function(formulas, input_syms, param_syms, ios.system.iv; expression = expression)
        return (f_oop=f_oop, f_ip=f_ip,
                states=state_syms,
                inputs=input_syms,
                params=param_syms)
    end
end

function transform_algebraic_equations(eqs::AbstractVector{Equation})
    eqs = deepcopy(eqs)
    for (i, eq) in enumerate(eqs)
        if eq.lhs isa Term && eq.lhs.op isa Differential
            continue
        end
        eqs[i] = 0 ~ eq.rhs - eq.lhs
    end
    return eqs
end

function reorder_by_states(eqs::AbstractVector{Equation}, states)
    @assert length(eqs) == length(states) "Numbers of eqs should be equal to states!"
    # for each state, collect the eq_idx which corresponds some states (implicit
    # agebraic) don't have special equations attached to them those are the "undused_idx"
    eq_idx::Vector{Union{Int, Nothing}} = [findfirst(x->isequal(s, lhs_var(x)), eqs) for s in states]
    unused_idx = reverse(setdiff(1:length(eqs), eq_idx))
    for i in 1:length(eq_idx)
        if eq_idx[i] === nothing
            eq_idx[i] = pop!(unused_idx)
        end
    end
    @assert sort(unique(eq_idx)) == 1:length(eqs) "eq_idx should contain all idx!"
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

function all_static(eqs::AbstractVector{Equation})
    tupels = eq_type.(eqs)
    if all(first.(tupels) .== :explicit_algebraic)
        lhs = last.(tupels)
        rhs = vcat([get_variables(eq.rhs) for eq in eqs]...)
        if isempty(Set(lhs) ∩ Set(rhs))
            return true
        end
    end
    return false
end
