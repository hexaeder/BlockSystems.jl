using ModelingToolkit
using ModelingToolkit: value
using IOSystems
using IOSystems.PDInterface
using PowerDynamics
using Plots

# HACK: call value in order to get rid of enclosing `Num <: Type`
(t, H, P, D, Ω) = value.(@parameters t H P D Ω)
(i,) = value.(@parameters i(t))       # input
(u_r, u_i, ω) = value.(@variables u_r(t) u_i(t) ω(t)) #outputs
Dt = Differential(t)

eqs = [Dt(u_r) ~ - u_i * ω,
       Dt(u_i) ~   u_r * ω,
       Dt(ω) ~ (P - D*ω - real((u_r + im*u_i)  * conj(i))) * Ω / (2*H)]

iob = IOBlock(eqs, [i], [u_r, u_i, ω], name=:swingnode)
# iob.inputs = i(t)
# iob.outputs = u(t), ω(t)
# iob.istates = empty
# iob.iparams = D, P, Ω, H

p_dict = Dict(H=>1.0,
              P=>1.0,
              D=>1.0,
              Ω=>1.0)
ion = IONode(iob, p_dict)

construct_vertex(ion)
