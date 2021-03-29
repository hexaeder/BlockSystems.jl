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
- `f_states`: define states=(istates ∪ outputs) which should appear first
- `f_inputs`: define inputs which should appear first
- `f_params`: define parameters which should appear first
- `f_rem_states`: define removed states algebraic state order
- `expression=Val{false}`: toggle expression and callable function output
- `warn=true`: toggle warnings for missing `f_*` parameters

Returns an named tuple with the fields
- for `type=:ode`:
  - `f_ip` in-place function `f(dstates, states, inputs, params, iv)`
  - `f_oop` out-of-place function `f(states, inputs, params, iv) => dstates`
- for `type=:static`:
  - `f_ip` in-place function `f(states, inputs, params, iv)`
  - `f_oop` out-of-place function `f(inputs, params, iv) => states`
- always:
- `massm` mass matrix of the system (`nothing` if :static)
- `states` symbols of states (in order)
- `inputs` symbols of inputs (in order)
- `params` symbols of parameters (in order)
- `rem_states` symbols of removed states (in order)
- `g_ip`, `g_oop` functions `g((opt. out), states, inputs, params, iv)` to calculate the
  removed states (substituted expl. algebraic equations). `nothing` if empty.
"""
function generate_io_function(ios::AbstractIOSystem; f_states = [], f_inputs = [],
                              f_params = [], f_rem_states = [],
                              expression = Val{false}, verbose=false, type=:auto, warn=true)
    if ios isa IOSystem
        @info "Transform given system $(ios.name) to block"
        ios = connect_system(ios, verbose=verbose)
    end

    # f_* may be given in namepsace version or as symbols
    f_states = prepare_f_vector(ios, f_states)
    f_inputs = prepare_f_vector(ios, f_inputs)
    f_params = prepare_f_vector(ios, f_params)
    f_rem_states = prepare_f_vector(ios, f_rem_states)

    @check Set(f_states) ⊆ (Set(ios.outputs) ∪ Set(ios.istates)) "f_states !⊆ (outputs ∪ istates)"
    @check Set(f_inputs) ⊆ Set(ios.inputs) "f_inputs !⊆ inputs"
    @check Set(f_params) ⊆ Set(ios.iparams) "f_params !⊆ iparams"
    @check Set(f_rem_states) ⊆ Set(ios.removed_states) "f_rem_states !⊆ removed_states"

    # enforce ordering of states, inputs and params
    states = vcat(f_states, ios.outputs, ios.istates) |> unique
    inputs = vcat(f_inputs, ios.inputs) |> unique
    params = vcat(f_params, ios.iparams) |> unique
    rem_states = vcat(f_rem_states, ios.removed_states) |> unique

    # warning
    if warn && (length(f_states) != length(states) ||
                length(f_inputs) != length(inputs) ||
                length(f_params) != length(params) ||
                length(f_rem_states) != length(rem_states))
        @warn "The ordering of the states/inputs/params/rem_states might change from run to run. Therefore it is highly recommend to provide all variables in the f_* arguments" states inputs params
    end

    # reorder the equations to get du in the right order
    eqs = reorder_by_states(get_eqs(ios.system), states)
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
        mass_matrix = nothing
    else
        throw(ArgumentError("Unknown type $type"))
    end

    # substitute x(t) by x for all terms
    state_syms = makesym.(states, states=[])
    input_syms = makesym.(inputs, states=[])
    param_syms = makesym.(params, states=[])
    rem_state_syms = makesym.(rem_states, states=[])

    sub = merge(Dict(states .=> state_syms),
                Dict(inputs .=> input_syms),
                Dict(params .=> param_syms))
    formulas = [substitute(eq.rhs, sub) for eq in eqs]

    # generate functions
    if type == :ode
        f_oop, f_ip = build_function(formulas, state_syms, input_syms, param_syms, get_iv(ios); expression = expression)
    elseif type == :static
        f_oop, f_ip = build_function(formulas, input_syms, param_syms, get_iv(ios); expression = expression)
    end

    # generate functions for removed states
    if isempty(rem_states)
        g_oop = nothing; g_ip = nothing
    else
        rem_eqs = reorder_by_states(ios.removed_eqs, rem_states)
        verbose && @info "Reordered removed eqs" rem_eqs rem_states

        rem_formulas = [substitute(eq.rhs, sub) for eq in rem_eqs]
        g_oop, g_ip = build_function(rem_formulas, state_syms, input_syms, param_syms, get_iv(ios); expression = expression)
    end

    return (f_oop=f_oop, f_ip=f_ip,
            massm=mass_matrix,
            states=state_syms,
            inputs=input_syms,
            params=param_syms,
            g_oop=g_oop, g_ip=g_ip,
            rem_states=rem_state_syms)
end

"""
    prepare_f_vector(iob, vector)

Prepare the user given variable lists which should appear first.
If `Symbol` type is given, it will look up the symbol as in
`iob.:sym`. The leading namespace of `iob` is removed from all variables
if it exists.
Returns new vector of `Symbolic`.
"""
function prepare_f_vector(iob::IOBlock, vector)
    newvec = Vector{Symbolic}(undef, length(vector))
    for (i, v) in enumerate(vector)
        sym = v isa Symbol ? getproperty(iob, v) : value(v)
        newvec[i] = remove_namespace(iob.name, sym)
    end
    return newvec
end

function transform_algebraic_equations(eqs::AbstractVector{Equation})
    eqs = deepcopy(eqs)
    for (i, eq) in enumerate(eqs)
        if eq.lhs isa Term && operation(eq.lhs) isa Differential
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
        if eqs[i].lhs isa Term && operation(eqs[i].lhs) isa Differential
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
