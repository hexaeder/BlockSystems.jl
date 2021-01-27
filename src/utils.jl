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
    return :(if !$(esc(cond)); $print; throw(ArgumentError($msg)) end)
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
    is_explicit_algebraic(eq::Equation)

True if lhs is a single symbol x and x ∉ rhs!
"""
function is_explicit_algebraic(eq::Equation)
    if eq.lhs isa Sym || eq.lhs isa Term && !(eq.lhs.op isa Differential)
        vars = get_variables(eq.lhs)
        @assert length(vars) == 1
        return vars[1] ∉ Set(get_variables(eq.rhs))
    else
        return false
    end
end
