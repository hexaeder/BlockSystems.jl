export connect_system, remove_superfluous_states, substitute_algebraic_states, substitute_derivatives
export simplify_eqs, set_input, make_input, make_iparam, replace_vars

"""
$(SIGNATURES)

Recursively transform `IOSystems` to `IOBlocks`.

- substitute inputs with connected outputs
- try to eliminate equations for internal states which are not used to calculate the specified outputs of the system.
- try to eliminate explicit algebraic equations (i.e. outputs of internal blocks) by substituting each occurrence
  with their rhs. Explicit algebraic states which are marked as system outputs won't be removed.

Arguments:
- `ios`: system to connect
- `verbose=false`: toggle verbosity (show equations at different steps)
- `remove_superflous_states=true`: toggle whether the system should try to get rid of unused states
- `substitute_algebraic_states=true`: toggle whether the algorithm tries to get rid of explicit algebraic equations
- `substitute_derivatives=true`: toggle whether to expand all derivatives and try to substitute them
- `simplify_eqs=true`: toggle simplification of all equations at the end
"""
function connect_system(ios::IOSystem;
                        verbose=false,
                        simplify_eqs=true,
                        remove_superflous_states=false,
                        substitute_algebraic_states=true,
                        substitute_derivatives=true,
                        warn=WARN[])
    # recursive connect all subsystems
    for (i, subsys) in enumerate(ios.systems)
        if subsys isa IOSystem
            ios.systems[i] = connect_system(subsys;
                                            verbose,
                                            simplify_eqs,
                                            remove_superflous_states,
                                            substitute_algebraic_states,
                                            substitute_derivatives,
                                            warn)
        end
    end
    eqs = vcat([namespace_equations(iob.system) for iob in ios.systems]...)
    removed_eqs = vcat([namespace_rem_eqs(iob) for iob in ios.systems]...)

    verbose && @info "Transform IOSystem $(ios.name) to IOBlock" ios.name ios.inputs ios.outputs ios.connections eqs

    # get rid of closed inputs by substituting output states
    substitutions = reverse.(ios.connections)
    substitute_all_rhs!(eqs, substitutions)
    substitute_all_rhs!(removed_eqs, substitutions)

    # connections of type a + b => might introduce terms inside differentials, needs expansion
    diffs = rhs_differentials(eqs)
    if !isempty(diffs)
        diff_expansion_rules = Dict(diffs .=> expand_derivatives.(diffs))
        substitute_all_rhs!(eqs, diff_expansion_rules)
    end

    # keep connections around in removed states
    for con in ios.connections
        push!(removed_eqs, con.second ~ con.first)
    end

    verbose && @info "substitute inputs with outputs" eqs

    # apply the namespace transformations
    promotion_rules = ios.namespace_map
    eqs = map(eq->eqsubstitute(eq, promotion_rules), eqs)
    removed_eqs  = map(eq->eqsubstitute(eq, promotion_rules), removed_eqs)

    block = IOBlock(ios.name, eqs, ios.inputs, ios.outputs, removed_eqs; iv=get_iv(ios), warn)

    if remove_superflous_states
        block = BlockSystems.remove_superflous_states(block; verbose, warn)
    end

    if substitute_derivatives
        block = BlockSystems.substitute_derivatives(block; verbose, warn)
    end

    if substitute_algebraic_states
        block = BlockSystems.substitute_algebraic_states(block; verbose, warn)
    end

    if simplify_eqs
        block = BlockSystems.simplify_eqs(block; warn)
    end

    return block
end


"""
    remove_superfluous_states(iob::IOBlock; verbose=false, warn=WARN[])

This function removes equations from block, which are not used in order to
generate the outputs. It looks for equations which have no path to the outputs
equations in the dependency graph. Returns a new IOBlock.

The removed equations will be not available as removed equations of the new IOBlock

TODO: Maybe we should try to reduce the inputs to.
"""
function remove_superfluous_states(iob::IOBlock; verbose=false, warn=WARN[])
    iv = get_iv(iob)
    outputs = iob.outputs

    neweqs = deepcopy(equations(iob))
    sys = ODESystem(neweqs, iv; name=:tmp) # will be used for the dependency graph
    neweqs = get_eqs(sys) # the ODESystem might reorder the equations
    # generate dependency graph
    graph = eqeq_dependencies(asgraph(sys), variable_dependencies(sys))
    # find 'main' eq for each output
    output_idx = [findfirst(x->o ∈ Set(get_variables(x.lhs)), neweqs) for o in outputs]

    if any(isnothing, output_idx)
        verbose && @info "Can't remove superfluous states if outputs implicitly defined."
        return neweqs
    end

    # if there is no path from equation to output equation is not necessary
    removable = []
    for eq_node in 1:length(neweqs)
        if !any(has_path(graph, eq_node, out_node) for out_node in output_idx)
            push!(removable, eq_node)
        end
    end

    removed_eqs = neweqs[removable]
    deleteat!(neweqs, sort(removable))

    verbose && @info "Removed superflous states with equations" removed_eqs

    IOBlock(iob.name, neweqs, iob.inputs, iob.outputs, iob.removed_eqs; iv, warn)
end


"""
    substitute_algebraic_states(iob::IOBlock; verbose=false, warn=WARN[])

Reduces the number of equations by substituting explicit algebraic equations.
Returns a new IOBlock with the reduced equations. The removed eqs are stored
together with the previous `removed_eqs` in the new IOBlock.
Won't reduce algebraic states which are labeled as `output`.
"""
function substitute_algebraic_states(iob::IOBlock; verbose=false, warn=WARN[])
    reduced_eqs = deepcopy(equations(iob))

    (rules, removable) = _algebraic_substitution_rules(reduced_eqs; skip=Set(iob.outputs))

    # subsitute all the equations, remove substituted
    for (i, eq) in enumerate(reduced_eqs)
        eq = substitute_rhs(eq, rules)
        reduced_eqs[i] = _transform_implicit_algebraic(eq; trysolve=true, verbose)
    end

    # also substitute in the already removed_eqs
    removed_eqs = deepcopy(iob.removed_eqs)
    for (i, eq) in enumerate(removed_eqs)
        removed_eqs[i] = substitute_rhs(eq, rules)
    end

    verbose && @info "Substituted algebraic states:" rules

    # append the knows removed wqs with the newly removed eqs
    append!(removed_eqs, reduced_eqs[removable])
    # remove removable equations from reduced_eqs
    deleteat!(reduced_eqs, sort(removable))

    IOBlock(iob.name, reduced_eqs, iob.inputs, iob.outputs, removed_eqs; iv=get_iv(iob), warn)
end

"""
    _algebraic_substitution_rules(eqs; skip=nothing)

Extract substitution rules for explicit algebraic equations in given equations.
Makes sure that those substitutions are pairwise cycle free.

Returns Tuple of substitution dict and corresponding indices of algebraic equations
in the given equations.
"""
function _algebraic_substitution_rules(eqs; skip=nothing)
    # only consider states for reduction which are explicit algebraic and not in outputs
    condition = eq -> begin
        (type, var) = eq_type(eq)
        type == :explicit_algebraic && (skip===nothing || var ∉ skip)
    end
    algebraic_idx = findall(condition, eqs)

    rules, keep = _uncouple_algebraic_equations(eqs[algebraic_idx])
    removable = deleteat!(algebraic_idx, keep)

    return (rules, removable)
end


"""
    substitute_derivatives(iob::IOBlock; verbose=false, warn=WARN[])

Expand all derivatives in the RHS of the system. Try to substitute
in the lhs with their definition.

I.e.

    D(o) ~ 1 + D(i)   =>  D(o) ~ 2 + a
    D(i) ~ 1 + a          D(i) ~ 1 + a

Process happens in multiple steps:
- try to find explicit equation for differential
- if none found try to recursively substitute inside differential with known algebraic states
- expand derivatives and try again to substitute with known differentials

"""
function substitute_derivatives(iob::IOBlock; verbose=false, warn=WARN[])
    diffs = rhs_differentials(iob)

    # if there are none just return the old block
    if isempty(diffs)
        verbose && @info "No rhs derivatives in block!"
        return iob
    end
    verbose && @info "Substitute rhs derivatives..."

    eqs = deepcopy(equations(iob))

    algebraic_subs = Dict{Symbolic, Any}()
    function lazy_alg_subs()
        if isempty(algebraic_subs)
            verbose && println("      ...lazily initialized algebraic substitutions.")
            merge!(algebraic_subs, _algebraic_substitution_rules(eqs)[1])
        end
        return algebraic_subs
    end

    rules = Dict()
    known_differentials = Dict(eq.lhs => eq.rhs for eq in eqs if istree(eq.lhs) && operation(eq.lhs) isa Differential)
    for diff in diffs
        if haskey(known_differentials, diff)
            rules[diff] = known_differentials[diff]
            verbose && println("    ", diff, " => ", rules[diff])
        else
            substituted = substitute(diff, lazy_alg_subs())
            expanded = expand_derivatives(substituted)
            sub_known = substitute(expanded, known_differentials)

            rules[diff] = sub_known
            verbose && println("    ", diff, " => ", rules[diff])
            if verbose && !isempty(_collect_differentials(sub_known))
                println("        └ could not resolve this one!")
            end
        end
    end

    eqs = deepcopy(equations(iob))
    rem_eqs = deepcopy(iob.removed_eqs)

    for (i, eq) in enumerate(eqs)
        eq = substitute_rhs(eq, rules)
        eqs[i] = _transform_implicit_algebraic(eq; trysolve=true, verbose)
    end
    for (i, eq) in enumerate(rem_eqs)
        eq = substitute_rhs(eq, rules)
        rem_eqs[i] = _transform_implicit_algebraic(eq; trysolve=true, verbose)
    end

    newblock = IOBlock(iob.name, eqs, iob.inputs, iob.outputs, rem_eqs; iv=get_iv(iob), warn)

    return newblock
end


"""
    simplify_eqs(iob::IOBlock; verbose=false, warn=WARN[])

Simplify eqs and removed eqs and return new IOBlock.
"""
function simplify_eqs(iob::IOBlock; verbose=false, warn=WARN[], hotfix=true)
    verbose && @info "Simplify iob equations..."
    simplified_eqs = simplify.(equations(iob))
    simplified_rem_eqs = simplify.(iob.removed_eqs)

    # TODO: temporary fix until the metadata issues are solved in MetaTheory
    if hotfix
        missing_metadata1 = check_metadata(simplified_eqs)
        if !isempty(missing_metadata1)
            warn && WARN_SIMPLIFY[] && @warn "Simplification of equations of $(iob.name) lead to missing metadata of $missing_metadata1. Skip!"
            simplified_eqs = equations(iob)
        end
        missing_metadata2 = check_metadata(simplified_rem_eqs)
        if !isempty(missing_metadata2)
            warn && WARN_SIMPLIFY[] && @warn "Simplification of removed equations of $(iob.name) lead to missing metadata of $missing_metadata2. Skip!"
            simplified_rem_eqs = iob.removed_eqs
        end
    end
    IOBlock(iob.name, simplified_eqs, iob.inputs, iob.outputs, simplified_rem_eqs; iv=get_iv(iob), warn)
end

"""
    replace_vars(blk::IOBlock, p::Dict; warn=WARN[])
    replace_vars(blk::IOBlock, p::Pair; warn=WARN[])
    replace_vars(blk::IOBlock; warn=WARN[], v1=val1, v2=val2)

Replace variables, either rename them by giving a new `Symbol` or replace them by actual numerical values
(only possible for inputs and iparams). Returns new `IOBlock`.

Keys of dict can be either `Symbols` or the `Symbolic` subtypes. I.e. `blk.u => 1.0` is as valid as `:u => 1.0`.

    replace_vars(blk; π = 3.14, foo = :bar) # set blk.π to number and rename foo
    replace_vars(blk, Dict(:π => 3.14, :foo = :bar))
    replace_vars(blk, blk.π => 3.14)
"""
replace_vars(blk::IOBlock, s...; warn=WARN[]) = length(s) > 1 ? replace_vars(blk, Dict(s); warn) : replace_vars(blk, Dict(only(s)); warn)
replace_vars(blk::IOBlock; warn=WARN[], kwargs...) = replace_vars(blk, Dict(kwargs); warn)
function replace_vars(blk::IOBlock, dict::Dict; warn=WARN[])
    subs = Dict{Symbolic, Any}()
    validp  = Set(blk.iparams)
    validin = Set(blk.inputs)
    closedinputs = Set()
    for (k,v) in dict
        try
            sym = getproperty(blk, k)
        catch
            warn && @warn "Symbol :$k not present in block. Skipped."
            continue
        end
        sym = remove_namespace(blk.name, sym)

        if v isa Number
            # replace path
            if sym ∈ validin
                push!(closedinputs, sym)
            else
                @check sym ∈ validp "Symbol $sym is neither iparam nor input of block. Can not repalace with number."
            end
            subs[sym] = dict[k]
        elseif v isa Symbol || v isa Symbolic
            # rename path
            val = rename(sym, v)
            subs[sym] = val
        end
    end
    eqs     = map(eq->eqsubstitute(eq, subs), equations(blk))
    rem_eqs = map(eq->eqsubstitute(eq, subs), blk.removed_eqs)

    newinputs = setdiff(blk.inputs, closedinputs)      # kick substitutions by number out
    newinputs = map(x->substitute(x, subs), newinputs) # replace new names
    newoutputs = map(x->substitute(x, subs), blk.outputs)

    IOBlock(blk.name, eqs, newinputs, newoutputs, rem_eqs; iv=get_iv(blk), warn)
end

"""
    set_input(blk::IOBlock, p::Pair; verbose=false)

Close an input of blk. Given as an pair of (input=>substitution). The input may
be given as an symbol (i.e. :a) or symbolic (i.e. blk.a). The substitution term
can be either a numer or a term of parameters (which will become internal
parameters).
"""
function set_input(blk::IOBlock, p::Pair; verbose=false)
    sym, val = p
    input = getproperty(blk, sym)
    inputname = getname(remove_namespace(blk.name, input))
    @check input ∈ Set(namespace_inputs(blk)) "Symbol $sym is no valid input!"
    iv = get_iv(blk)
    tmp, = @variables $inputname(iv)
    tmpblk = IOBlock([tmp ~ val], [], [tmp])
    sys = IOSystem([getproperty(tmpblk, inputname) => getproperty(blk, inputname)],
                    [tmpblk, blk];
                   outputs=namespace_outputs(blk), name=blk.name)
    return connect_system(sys; verbose)
end

"""
    make_input(blk::IOBlock, sym; warn=WARN[])

Change given sym from iparam to input.
"""
function make_input(blk::IOBlock, sym; warn=WARN[])
    var_nspcd = getproperty(blk, sym)
    var = remove_namespace(blk.name, var_nspcd)
    @check var ∈ Set(blk.iparams) "$sym is not an iparam of block."

    varname = getname(var)
    iv = get_iv(blk)
    tmp, = @parameters $varname(iv)

    sub = var => tmp
    eqs     = map(eq->eqsubstitute(eq, sub), equations(blk))
    rem_eqs = map(eq->eqsubstitute(eq, sub), blk.removed_eqs)

    IOBlock(blk.name, eqs, vcat(blk.inputs,tmp), blk.outputs, rem_eqs; iv, warn)
end

"""
    make_iparam(blk::IOBlock, sym; warn=WARN[])

Change given sym from input to iparam.
"""
function make_iparam(blk::IOBlock, sym; warn=WARN[])
    var_nspcd = getproperty(blk, sym)
    var = remove_namespace(blk.name, var_nspcd)
    @check var ∈ Set(blk.inputs) "$sym is not an input of block."

    varname = getname(var)
    iv = get_iv(blk)
    tmp, = @parameters $varname
    sub = var => tmp
    eqs     = map(eq->eqsubstitute(eq, sub), equations(blk))
    rem_eqs = map(eq->eqsubstitute(eq, sub), blk.removed_eqs)

    IOBlock(blk.name, eqs, filter(!isequal(var), blk.inputs), blk.outputs, rem_eqs; iv, warn)
end
