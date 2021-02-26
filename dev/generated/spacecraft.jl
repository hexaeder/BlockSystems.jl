using BlockSystems
using ModelingToolkit
@parameters t M F(t)
@variables x(t) v(t)
D = Differential(t)

spacecraft = IOBlock([D(v) ~ F/M, D(x) ~ v], # define the equation
                     [F], # inputs of the system
                     [x], # outputs of the system
                     name = :spacecraft)

@parameters K i(t)
@variables o(t)

prop = IOBlock([o ~ K*i], [i], [o], name = :prop)
nothing # hide

@parameters p(t) m(t)
@variables Δ(t)
diff = IOBlock([Δ ~ p - m], [p, m], [Δ], name=:diff)
nothing # hide

prop_c = IOSystem([diff.Δ => prop.i], # connect output of diff to input of prop
                  [diff, prop], # subsystems
                  name=:propc)

prop_c = IOSystem([diff.Δ => prop.i], [diff, prop],
                  namespace_map = [prop.o => o,
                                   diff.p => :target,
                                   diff.m => :feedback],
                  outputs = [o],
                  name=:propc)

prop_c_block = connect_system(prop_c)
nothing #hide

@variables altitude(t)
space_controller = IOSystem([prop_c.o => spacecraft.F, spacecraft.x => prop_c.feedback],
                            [prop_c, spacecraft],
                            namespace_map = [spacecraft.x => altitude],
                            outputs = [altitude])
# we want to reduce the space_controller to a block
space_controller = connect_system(space_controller)
@info "Variables of space_controller" space_controller space_controller.system.eqs

gen = generate_io_function(space_controller)
nothing # hide

gen = generate_io_function(space_controller, f_states=[altitude, v], f_params=[K, M])
@info "Generated function" gen.massm gen.states gen.inputs gen.params
nothing # hide

using Plots
using OrdinaryDiffEq
targetfun(t) = t>1.0 ? 1.0 : 0
odefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)
p = [0.5, 1.0] # K, m
u0 = [0.0, 0.0] # altitude, v
tspan = (0.0, 30.0)
prob = ODEProblem(odefun, u0, tspan, p)
sol = solve(prob, Tsit5())
plot(t->sol(t)[1],tspan..., label="altitude", title="proportional control")
plot!(t->targetfun(t),tspan..., label="target")

spacecraft = IOBlock([D(v) ~ F/M, D(x) ~ v],
                     [F],
                     [x,v],
                     name = :spacecraft)

prop_v = IOBlock(prop, name=:prop_v)
fdiff = IOBlock(diff, name=:fdiff)

space_controller = IOSystem([spacecraft.v => prop_v.i,
                             spacecraft.x => prop_c.feedback,
                             prop_c.o => fdiff.p,
                             prop_v.o => fdiff.m,
                             fdiff.Δ => spacecraft.F],
                            [prop_v, prop_c, fdiff, spacecraft],
                            namespace_map = [spacecraft.x => altitude],
                            outputs = [altitude])

space_controller = connect_system(space_controller)
gen = generate_io_function(space_controller, f_states=[altitude, v], f_params=[prop_c.K, M, prop_v.K])

odefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)
p = [1.0, 1.0, 0.5] # propc, M, propv
u0 = [0.0, 0.0] # altitude, v
tspan = (0.0, 30.0)
prob = ODEProblem(odefun, u0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol, vars=(0,1), label="altitude", title="better control")
plot!(t->targetfun(t),tspan..., label="target")

@parameters T a1(t) a2(t)
@variables Σ(t) altitude(t)
int = IOBlock([D(o) ~ 1/T * i - o], [i], [o], name=:int)
adder = IOBlock([Σ ~ a1 + a2], [a1, a2], [Σ], name=:add)

pi_c = IOSystem([diff.Δ => prop.i,
                 prop.o => int.i,
                 prop.o => adder.a1,
                 int.o => adder.a2],
                [diff, prop, int, adder],
                namespace_map = [diff.p => :target,
                                 diff.m => :feedback,
                                 adder.Σ => o],
                outputs = [o],
                name=:pi_c)
nothing # hide

space_controller = IOSystem([pi_c.o => spacecraft.F , spacecraft.x => pi_c.feedback],
                            [pi_c, spacecraft],
                            namespace_map = [spacecraft.x => altitude],
                            outputs = [altitude])
space_controller = connect_system(space_controller, verbose=false)
@info "Variables of space_controller" space_controller space_controller.system.eqs

gen = generate_io_function(space_controller, f_states=[altitude], f_params=[K, T, M])

odefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)
p = [0.5, -1.5, 1.0] # K, T, m
u0 = [0.0, 0.0, 0.0] # altitude, int.o, v
tspan = (0.0, 50.0)
prob = ODEProblem(odefun, u0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol, vars=(0,[ 1,2 ]), label=["altitude" "integrator"], title="PT1 controller")
plot!(t->targetfun(t),tspan..., label="target")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

