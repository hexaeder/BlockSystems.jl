export generate_io_function

"""
$(SIGNATURES)

Generate callable functions for an `AbstractIOSystem`. An `IOSystem` will be transformed to an `IOBlock` first.
At this level there is no more distinction between internal states and outputs:
`states=(istates ∪ outputs)`.

Arguments:
- `ios`: the system to build the function
optional:
- `type=:auto`: `:ode` or `:static`, determines the output of the function
- `f_states`: define states=(istates ∪ outputs) which should appear first
- `f_inputs`: define inputs which should appear first
- `f_params`: define parameters which should appear first
- `f_rem_states`: define removed states algebraic state order
- `expression=Val{false}`: toggle expression and callable function output
- `warn=WARN[]`: toggle warnings for missing `f_*` parameters
- `observed=false`: toggle creation ob "observed" function

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
- `obsf`: Only if `observed=true`. Function in the SciML "observed" style `(sym, u, p, t)`
"""
function generate_io_function(ios::AbstractIOSystem; f_states = [], f_inputs = [],
                              f_params = [], f_rem_states = [],
                              expression = Val{false}, verbose=false, type=:auto, observed=false, warn=WARN[])
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
    @check isempty(rhs_differentials(ios)) "RHS should not contain any differentials at this point."

    # enforce ordering of states, inputs and params
    states = vcat(f_states, ios.outputs, ios.istates) |> unique
    inputs = vcat(f_inputs, ios.inputs) |> unique
    params = vcat(f_params, ios.iparams) |> unique
    rem_states = vcat(f_rem_states, ios.removed_states) |> unique

    # warning
    if warn
        incomplete = false
        warnstr = ""
        if length(states)>1 && length(f_states) != length(states)
            incomplete = true
            warnstr *= "f_states, "
        end
        if length(inputs)>1 && length(f_inputs) != length(inputs)
            incomplete = true
            warnstr *= "f_inputs, "
        end
        if length(params)>1 && length(f_params) != length(params)
            incomplete = true
            warnstr *= "f_params, "
        end
        # don't warn for missing removed states ordering
        # if length(rem_states)>1 && length(f_rem_states) != length(rem_states)
        #     incomplete = true
        #     warnstr *= "f_rem_states, "
        # end
        if incomplete
            @warn "The ordering of the states might change in future versions. Therefore it is highly recommend to provide all variables in the f_* arguments. There are missing entrys in " * warnstr[begin:end-2] * "."
        end
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
    state_syms = strip_iv.(states, get_iv(ios))
    input_syms = strip_iv.(inputs, get_iv(ios))
    param_syms = strip_iv.(params, get_iv(ios))
    rem_state_syms = strip_iv.(rem_states, get_iv(ios))

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
        obsf = nothing
    else
        rem_eqs = reorder_by_states(ios.removed_eqs, rem_states)
        verbose && @info "Reordered removed eqs" rem_eqs rem_states

        rem_formulas = [substitute(eq.rhs, sub) for eq in rem_eqs]
        g_oop, g_ip = build_function(rem_formulas, state_syms, input_syms, param_syms, get_iv(ios); expression = expression)

        if observed
            @check isempty(input_syms) "Cannot create `observed` function if there are open inputs."
            obsdict = Dict{Symbol, Function}()
            for (idx, state) in pairs(rem_states)
                obs_oop, _  = build_function(rem_formulas[[idx]], state_syms, input_syms, param_syms, get_iv(ios); expression)
                obsdict[getname(state)] = obs_oop
                obsf = (sym, u, p, t) -> only(obsdict[sym](u, nothing, p, t))
            end
        else
            obsf = nothing
        end
    end

    return (f_oop=f_oop, f_ip=f_ip,
            massm=mass_matrix,
            states=state_syms,
            inputs=input_syms,
            params=param_syms,
            g_oop=g_oop, g_ip=g_ip,
            rem_states=rem_state_syms,
            obsf)
end

"""
    ODEFunction(iob::IOBlock; f_states=Symbol[], f_params=Symbol[], verbose=false)

Return an `ODEFunction` object with the corresponding mass matrix and variable names.
"""
function SciMLBase.ODEFunction(iob::IOBlock; f_states=Symbol[], f_params=Symbol[], verbose=false, warn=WARN[])
    @check isempty(iob.inputs) "all inputs must be closed"
    gen = generate_io_function(iob; f_states, f_params, verbose, warn, observed=true);

    # observed = (sym,u,p,t)->gen.g_oop(u,nothing,p,t)[findfirst(isequal(sym), Symbol.(gen.rem_states))]
    ODEFunction((du,u,p,t) -> gen.f_ip(du,u,nothing,p,t);
                mass_matrix = gen.massm,
                syms=Symbol.(gen.states),
                indepsym=get_iv(iob),
                observed=gen.obsf)
end

"""
    strip_iv(x, iv)

Strip functional dependency of the independent variable `x(iv) -> x`.
"""
function strip_iv(x::Symbolic, iv::Symbolic)
    if istree(x) && operation(x) isa Symbolic
        if (length(arguments(x)) != 1 || !isequal(arguments(x)[1], iv))
            error("Don't knowhow to handle expression $x")
        end
        return Sym{SymbolicUtils.symtype(x)}(tosymbol(operation(x)))
    else
        return x
    end
end
strip_iv(x::Num, iv::Num) = Num(strip_iv(value(x), value(iv)))

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
        if istree(eq.lhs) && operation(eq.lhs) isa Differential
            continue
        end
        eqs[i] = 0 ~ eq.rhs - eq.lhs
    end
    return eqs
end

function reorder_by_states(eqs::AbstractVector{Equation}, states)
    @assert length(eqs) == length(states) "Numbers of eqs should be equal to states!"
    # for each state, collect the eq_idx which corresponds some states (implicit
    # algebraic) don't have special equations attached to them those are the "undused_idx"
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
        if istree(eqs[i].lhs) && operation(eqs[i].lhs) isa Differential
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
