"""
    @check cond msg

If `cond` evaluates false throw `ArgumentError` and print evaluation of `cond`.
"""
macro check(cond, msg)
    str = _chkmsg(cond, msg)
    quote
        if !$(esc(cond))
            throw(ArgumentError($(str)))
        end
    end
end

"""
    @checkwarn [showwarn] cond msg

If `cond` evaluates false warn and print evaluation of `cond`.
If optional argument `showwarn=false` hide warning.
"""
macro checkwarn(cond, msg)
    str = _chkmsg(cond, msg)
    quote
        if !$(esc(cond))
            @warn $(str)
        end
    end
end

macro checkwarn(showarn, cond, msg)
    str = _chkmsg(cond, msg)
    quote
        if $(esc(showarn)) && !$(esc(cond))
            @warn $(str)
        end
    end
end

"""
    _chkmsg(cond, msg)

Create an expression which evaluates as a error string which contains the original msg and
some debug information based on the condition.
"""
function _chkmsg(cond, msg)
    head = lstrip(repr(cond), ':')
    if head[begin] == '(' && head[end] == ')'
        head = head[begin+1:end-1]
    end
    head = head * " evaluated false"

    args = Expr[]

    if cond.args[1] == :(!)
        _expr_repr_list!(args, cond.args[2].args[2:end])
    elseif cond.args[1] == :(⊆)
        push!(args, :("\n   ├ set diff = " * _shortrepr(setdiff($(esc(cond.args[2])), $(esc(cond.args[3]))))))
        _expr_repr_list!(args, cond.args[2:end])
    elseif cond.args[1] == :(⊇)
        push!(args, :("\n   ├ set diff = " * _shortrepr(setdiff($(esc(cond.args[3])), $(esc(cond.args[2]))))))
        _expr_repr_list!(args, cond.args[2:end])
    elseif cond.args[1] == :(uniquenames) || cond.args[1] == :(allunique)
        push!(args, :("\n   ├ duplicates = " * _shortrepr(duplicates($(esc(cond.args[2]))))))
        _expr_repr_list!(args, cond.args[2:end])
    else
        _expr_repr_list!(args, cond.args[2:end])
    end
    return :($(esc(msg)) * "\n  " * $head * $(args...))
end

function _expr_repr_list!(list, expressions)
    for (i,a) in enumerate(expressions)
        lhs = lstrip(repr(a), ':') # remove leading : on symbols
        if lhs[begin] == '(' && lhs[end] == ')'
            lhs = lhs[begin+1:end-1]
        end
        symbol = (i == length(expressions)) ? "└ " : "├ "
        push!(list, :("\n   " * $symbol * $lhs * " = " * _shortrepr($(esc(a)))))
    end
end

_shortrepr(x) = repr(x)
function _shortrepr(@nospecialize x::Set)
    isempty(x) && return "(empty)"
    m = match(r"\[(.*)\]\)$", repr(x)::String)
    return m[1]
end
function _shortrepr(@nospecialize x::AbstractArray)
    isempty(x) && return "(empty)"
    m = match(r"\[(.*)\]$", repr(x)::String)
    return m[1]
end

"""
    remove_namespace(namespace, x)

If first namespace of `x` is `namespace` then remove it.
```
remove_namespace(A, A₊B₊x) -> B₊x
remove_namespace(B, A₊B₊x) -> A₊B₊x
````
"""
remove_namespace(namespace, name::T) where T = T(replace(String(name), Regex("^$(namespace)₊") => ""))
remove_namespace(namespace, x::Sym) = rename(x, remove_namespace(namespace, x.name))
remove_namespace(namespace, x::Term) = rename(x, remove_namespace(namespace, operation(x).name))

"""
    remove_namespace(x)

Strips `x` of its first namespace. `A₊B₊x -> B₊x`
"""
remove_namespace(name::T) where T = T(replace(String(name), Regex("^(.+?)₊") => ""))
remove_namespace(x::Sym) = rename(x, remove_namespace(x.name))
remove_namespace(x::Term) = rename(x, remove_namespace(operation(x).name))

eqsubstitute(eq::Equation, rules) = substitute(eq.lhs, rules) ~ substitute(eq.rhs, rules)
substitute_rhs(eq::Equation, rules) = eq.lhs ~ substitute(eq.rhs, rules)

uniquenames(syms) = allunique(getname.(syms))

function duplicates(X)
    X = collect(X)
    dups = Set{eltype(X)}()
    for (i, x) in enumerate(X)
        for j in i+1:lastindex(X)
            if isequal(x, X[j])
                push!(dups, x)
                break
            end
        end
    end
    return dups
end

"""
    eq_type(eq::Equation)

Checks the type of the equation. Returns:
- `(:explicit_diffeq, lhs_variable)` for explicit differential equations
- `(:implicit_diffeq, lhs_variable)` for implicit differential equations
- `(:explicit_algebraic, lhs_variable)` for explicit algebraic equations
- `(:implicit_algebraic, lhs_variable)` for implicit algebraic equations

"""
function eq_type(eq::Equation)
    if eq.lhs isa Term && operation(eq.lhs) isa Differential
        vars = get_variables(eq.lhs)
        @check length(vars) == 1 "Diff. eq $eq has more than one variable in lhs!"
        return (:explicit_diffeq, vars[1])
    elseif eq.lhs isa Sym || eq.lhs isa Term
        vars = get_variables(eq.lhs)
        @check length(vars) == 1 "Algebraic eq $eq has more than one variable in lhs!"
        diffs = _collect_differentials(eq.rhs)
        if diffs != Set{SymbolicUtils.Symbolic}()
            if operation(first(diffs.dict)[1]) isa Differential
                return (:implicit_diffeq, vars[1])
            else
                throw(ArgumentError("Unknown equation type $eq"))
            end
        end 
        if vars[1] ∈ Set(get_variables(eq.rhs))
            return (:implicit_algebraic, vars[1])
        else
            return (:explicit_algebraic, vars[1])
        end
    elseif isequal(eq.lhs, 0)
        return (:implicit_algebraic, nothing)
    else
        throw(ArgumentError("Unknown equation type $eq"))
    end
end

"""
    lhs_var(eq::Equation)

Returns the variable on the lhs of the equation for equations.
"""
lhs_var(eq::Equation) = eq_type(eq)[2]


"""
    _transform_implicit_algebraic(eq; trysolve=true, verbose=false)

Transforms implicit algebraic equations with non-nohting lhs. If `trysolve` tries to solve them
for lhs. Otherwis just transforms to `0 ~ rhs - lhs`.
"""
function _transform_implicit_algebraic(eq; trysolve=true, verbose=false)
    (type, lhs_var) = eq_type(eq)
    if type === :implicit_algebraic && !isnothing(lhs_var)
        verbose && println("    Substitution resulted in implicit equation and was transformed!")
        verbose && println("      ├ ", eq)
        eq = 0 ~ simplify(eq.rhs - eq.lhs)
        if trysolve && lhs_var ∈ Set(get_variables(eq.rhs))
            try
                eq = lhs_var ~ Symbolics.solve_for(eq, lhs_var)
            catch e
                verbose && println("      ├ could not be resolved!")
            end
        end
        verbose && println("      └ $eq")
        return eq
    else
        return eq
    end
end


"""
    recursive_substitute(term, rules::Dict)

Apply substitutions until there is no more change in `term`.
"""
function recursive_substitute(term, rules::Dict)
    new_term = substitute(term, rules)
    while !isequal(new_term, term)
        term = new_term
        new_term = substitute(term, rules)
    end
    return new_term
end

"""
    rhs_differentials(iob::IOBlock)

Return Set of all differentials which are present in the rhs of the system.
"""
function rhs_differentials(iob)
    diffs = Set{SymbolicUtils.Symbolic}()
    for eq in equations(iob)
        _collect_differentials!(diffs, eq.rhs)
    end
    return diffs
end

_collect_differentials(ex) = _collect_differentials!(Set{SymbolicUtils.Symbolic}(), ex)

function _collect_differentials!(found, ex)
    if SymbolicUtils.istree(ex)
        if operation(ex) isa Differential
            push!(found, ex)
        else
            for arg in arguments(ex)
                _collect_differentials!(found, arg)
            end
        end
    end
    return found
end
