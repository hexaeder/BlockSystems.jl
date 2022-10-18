export rename_vars, set_p

"""
    rename_vars(blk::IOBLock; warn=WARN[], kwargs...)
    rename_vars(blk::IOBlock, subs::Dict{Symbolic,Symbolic}; warn=WARN[])

Returns new IOBlock which is similar to blk but with new variable names.
Variable renaming should be provided as keyword arguments, i.e.

    rename_vars(blk; x=:newx, k=:knew)

to rename `x(t)=>newx(t)` and `k=>knew`. Substitutions can be also provided as
dict of `Symbolic` types (`Sym`s and `Term`s).
"""
function rename_vars(blk::IOBlock; warn=WARN[], kwargs...)
    substitutions = Dict{Symbolic, Symbolic}()
    for pair in kwargs
        key = remove_namespace(blk.name, getproperty(blk, pair.first))
        val = rename(key, pair.second)
        substitutions[key] = val
    end
    rename_vars(blk, substitutions; warn)
end

function rename_vars(blk::IOBlock, subs::Dict{Symbolic,Symbolic}; warn=WARN[])
    Base.depwarn("The functions `rename_vars` and `set_p` have been deprecated in favor of `replace_vars` which unifies the functionality.", :rename_vars)
    isempty(subs) && return blk
    eqs     = map(eq->eqsubstitute(eq, subs), get_eqs(blk.system))
    rem_eqs = map(eq->eqsubstitute(eq, subs), blk.removed_eqs)
    inputs  = map(x->substitute(x, subs), blk.inputs)
    outputs = map(x->substitute(x, subs), blk.outputs)
    IOBlock(blk.name, eqs, inputs, outputs, rem_eqs; iv=get_iv(blk), warn)
end


"""
    set_p(blk::IOBlock, p::Dict; warn=WARN[])
    set_p(blk::IOBlock, p::Pair; warn=WARN[])
    set_p(blk::IOBlock; warn=WARN[], p1=val1, p2=val2)

Substitutes certain parameters by actual Float values. Returns an IOBlock without those parameters. Can be used
for both iparams as well as inputs!

Keys of dict can be either `Symbols` or the `Symbolic` subtypes. I.e. `blk.u => 1.0` is as valid as `:u => 1.0`.
"""
function set_p(blk::IOBlock, p::Dict; warn=WARN[])
    Base.depwarn("The functions `rename_vars` and `set_p` have been deprecated in favor of `replace_vars` which unifies the functionality.", :set_p)
    subs = Dict{Symbolic, Float64}()
    validp  = Set(blk.iparams)
    validin = Set(blk.inputs)
    closedinputs = Set()
    for k in keys(p)
        try
            sym = getproperty(blk, k)
        catch
            warn && @warn "Symbol $k not present in block. Skipped."
            continue
        end
        sym = remove_namespace(blk.name, sym)
        if sym ∈ validin
            push!(closedinputs, sym)
        else
            @check sym ∈ validp "Symbol $sym is not iparam of block"
        end
        @check p[k] isa Number "p value has to be a number! Maybe you are looking for `rename_vars` instead?"
        subs[sym] = p[k]
    end
    eqs     = map(eq->eqsubstitute(eq, subs), equations(blk))
    rem_eqs = map(eq->eqsubstitute(eq, subs), blk.removed_eqs)
    IOBlock(blk.name, eqs, setdiff(blk.inputs,closedinputs), blk.outputs, rem_eqs; iv=get_iv(blk), warn)
end

set_p(blk::IOBlock, p...; warn=WARN[]) = length(p) > 1 ? set_p(blk, Dict(p); warn) : set_p(blk, Dict(only(p)); warn)
set_p(blk::IOBlock; warn=WARN[], kwargs...) = set_p(blk, Dict(kwargs); warn)
