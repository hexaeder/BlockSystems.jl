using IOSystems
using ModelingToolkit
@parameters t M F(t)
@variables x(t) v(t)
@derivatives D'~t

spacecraft = IOBlock([D(v) ~ F/M, D(x) ~ v], # define the equation
                     [F], # inputs of the system
                     [x], # outputs of the system
                     name = :spacecraft)
@show spacecraft.inputs
@show spacecraft.outputs
@show spacecraft.istates
nothing # hide

@parameters K i(t)
@variables o(t)

prop = IOBlock([o ~ K*i], [i], [o], name = :prop)
nothing # hide

@parameters p(t) m(t)
@variables Δ(t)
diff = IOBlock([Δ ~ p - m], [p, m], [Δ], name=:diff)
nothing # hide

prop_c = IOSystem([prop.i => diff.Δ], # connect input of prop to output of diff
                  [diff, prop], # subsystems
                  name=:propc)
@info "namespace mapping" prop_c.inputs_map prop_c.istates_map prop_c.outputs_map

@parameters target(t) feedback(t)
prop_c = IOSystem([prop.i => diff.Δ], [diff, prop],
                  outputs_map = [prop.o => o],
                  inputs_map = [diff.p => target, diff.m => feedback],
                  name=:propc)
@info "namespace mapping" prop_c.inputs_map prop_c.istates_map prop_c.outputs_map

prop_c_block = connect_system(prop_c)
nothing #hide

@variables altitude(t)
space_controller = IOSystem([spacecraft.F => prop_c.o, prop_c.feedback => spacecraft.x],
                            [prop_c, spacecraft],
                            outputs_map = [spacecraft.x => altitude])
# we want to reduce the space_controller to a block
space_controller = connect_system(space_controller)
@info "Variables of space_controller" space_controller.inputs space_controller.outputs space_controller.istates space_controller.iparams space_controller.system.eqs

gen = generate_io_function(space_controller)
nothing # hide

gen = generate_io_function(space_controller, first_states=[altitude])
@info "Generated function" gen.massm gen.states gen.inputs gen.params
nothing # hide

using Plots
using DifferentialEquations
targetfun(t) = t>1.0 ? 1.0 : 0
odefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)
p = [0.5, 1.0] # K, m
u0 = [0.0, 0.0] # altitude, v
tspan = (0.0, 30.0)
prob = ODEProblem(odefun, u0, tspan, p)
sol = solve(prob)
plot(t->sol(t)[1],tspan..., label="altitude", title="proportional control")
plot!(t->targetfun(t),tspan..., label="target")

spacecraft = IOBlock([D(v) ~ F/M, D(x) ~ v],
                     [F],
                     [x,v],
                     name = :spacecraft)

prop_v = IOBlock(prop, name=:prop_v)
fdiff = IOBlock(diff, name=:fdiff)

space_controller = IOSystem([prop_v.i => spacecraft.v,
                             prop_c.feedback => spacecraft.x,
                             fdiff.p => prop_c.o,
                             fdiff.m => prop_v.o,
                             spacecraft.F => fdiff.Δ],
                            [prop_v, prop_c, fdiff, spacecraft],
                            outputs_map = [spacecraft.x => altitude])

space_controller = connect_system(space_controller)
gen = generate_io_function(space_controller, first_states=[altitude])

odefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)
p = [1.0, 1.0, 0.5] # propc, M, propv
u0 = [0.0, 0.0] # altitude, v
tspan = (0.0, 30.0)
prob = ODEProblem(odefun, u0, tspan, p)
sol = solve(prob, dtmax=0.1)
plot(t->sol(t)[1],tspan..., label="altitude", title="better control")
plot!(t->targetfun(t),tspan..., label="target")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

