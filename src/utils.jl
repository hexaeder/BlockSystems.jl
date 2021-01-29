"""
   check(cond, msg)

If `cond` evaluates false throw `ArgumentError` and print evaluation of `cond`.
TODO: Proper Output
"""
macro check(cond::Expr, msg)
    variables = ()
    for a in cond.args[2:end]
        if a isa Expr || a isa Symbol
            variables = (esc(a), variables...)
        end
    end
    print = :(@show($(repr(cond)),$(variables...)))
    return :(if !$(esc(cond)); $print; throw(ArgumentError(esc($msg))) end)
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
remove_namespace(namespace, x::Term) = rename(x, remove_namespace(namespace, x.op.name))

"""
    remove_namespace(x)

Strips `x` of its first namespace. `A₊B₊x -> B₊x`
"""
remove_namespace(name::T) where T = T(replace(String(name), Regex("^(.+?)₊") => ""))
remove_namespace(x::Sym) = rename(x, remove_namespace(x.name))
remove_namespace(x::Term) = rename(x, remove_namespace(x.op.name))

eqsubstitute(eq::Equation, rules) = substitute(eq.lhs, rules) ~ substitute(eq.rhs, rules)

uniquenames(syms) = allunique(getname.(syms))

"""
    eq_type(eq::Equation)

Checks the type of the equation. Returns:
- `(:diffeq, lhs_variable)` for differential equations
- `(:explicit_algebraic, lhs_variable)` for explicit algebraic equations
- `(:implicit_algebraic, lhs_variable)` for implicit algebraic equations
"""
function eq_type(eq::Equation)
    if eq.lhs isa Term && eq.lhs.op isa Differential
        vars = get_variables(eq.lhs)
        @check length(vars) == 1 "Diff. eq $eq has more than on variable in lhs!"
        return (:diffeq, vars[1])
    elseif eq.lhs isa Sym || eq.lhs isa Term
        vars = get_variables(eq.lhs)
        @check length(vars) == 1 "Algebraic eq $eq has more than on variable in lhs!"
        if vars[1] ∈ Set(get_variables(eq.rhs))
            return (:implicit_algebraic, vars[1])
        else
            return (:explicit_algebraic, vars[1])
        end
    elseif isequal(eq.lhs, 0)
        return (:implicit_algebraic, nothing)
    else
        throw(ArgumentError("Uknnown equation type $eq"))
    end
end

"""
    lhs_var(eq::Equation)

Returns the variable on the lhs of the equation for equations.
"""
lhs_var(eq::Equation) = eq_type(eq)[2]


"""
    recusive_substitute(term, rules::Dict)

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
