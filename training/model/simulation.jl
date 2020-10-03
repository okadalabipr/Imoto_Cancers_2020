module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./observable.jl")

using .C
using .V
import .Main

using Sundials
using SteadyStateDiffEq

# Options for ODE solver
const ABSTOL = 1e-9
const RELTOL = 1e-9
const DTMIN = 1e-8

const dt = 1.0
const tspan = (0.0,7200.0)
const t2 = collect(tspan[1]:dt:tspan[end])
const t = collect(tspan[1]:1.0:tspan[end])./60.0

const conditions = [
    "MCF7_EGF", "MCF7_HRG",
    "BT474_EGF", "BT474_HRG",
    "MDAMB231_EGF", "MDAMB231_HRG"
]

simulations = Array{Float64,3}(
    undef, length(observables), length(t), length(conditions)
)

function solveode(
        f::Function,
        u0::Vector{Float64},
        t::Vector{Float64},
        p::Vector{Float64})
    local sol::ODESolution{}, is_successful::Bool
    try
        prob = ODEProblem(f,u0,(t[1],t[end]),p)
        sol = solve(
            prob,CVODE_BDF(),
            abstol=ABSTOL,reltol=RELTOL,dtmin=DTMIN,saveat=dt,verbose=false
        )
        is_successful = ifelse(sol.retcode === :Success, true, false)
    catch
        is_successful = false
    finally
        if !is_successful
            GC.gc()
        end
    end
    return is_successful ? sol : nothing
end


function get_steady_state(
        f::Function,
        u0::Vector{Float64},
        p::Vector{Float64})::Vector{Float64}
    local sol::SteadyStateSolution{}, is_successful::Bool
    try
        prob = ODEProblem(diffeq,u0,(0.0,Inf),p)
        prob = SteadyStateProblem(prob)
        sol = solve(
            prob,
            DynamicSS(
                CVODE_BDF();abstol=ABSTOL,reltol=RELTOL
            ),
            dt=dt,dtmin=DTMIN,verbose=false
        )
        is_successful = ifelse(sol.retcode === :Success, true, false)
    catch
        is_successful = false
    finally
        if !is_successful
            GC.gc()
        end
    end
    return is_successful ? sol.u : []
end


function simulate!(p::Vector{Float64}, u0::Vector{Float64})
    p_in_mcf7::Vector{Float64} = copy(p)
    u0_in_mcf7::Vector{Float64} = copy(u0)
    for (i, condition) in enumerate(conditions)
        if condition == "MCF7_EGF"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.E] = 10.0
                u0[V.H] = 0.0
            end
        elseif condition == "MCF7_HRG"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.E] = 0.0
                u0[V.H] = 10.0
            end
        elseif condition == "BT474_EGF"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p,u0) = Main.mul_ratio2mcf7!(p,u0,"BT-474")
            u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.E] = 10.0
                u0[V.H] = 0.0
            end
        elseif condition == "BT474_HRG"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p,u0) = Main.mul_ratio2mcf7!(p,u0,"BT-474")
            u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.E] = 0.0
                u0[V.H] = 10.0
            end
        elseif condition == "MDAMB231_EGF"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p,u0) = Main.mul_ratio2mcf7!(p,u0,"MDA-MB-231")
            u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.E] = 10.0
                u0[V.H] = 0.0
            end
        elseif condition == "MDAMB231_HRG"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p,u0) = Main.mul_ratio2mcf7!(p,u0,"MDA-MB-231")
            u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.E] = 0.0
                u0[V.H] = 10.0
            end
        end
        sol = solveode(diffeq,u0,t2,p)
        if sol === nothing
            return false
        else
            @inbounds @simd for j in eachindex(t)
                simulations[observables_index("Phosphorylated_Akt"),j,i] = (
                    sol.u[j][V.Aktstar]
                )
                simulations[observables_index("Phosphorylated_ERK"),j,i] = (
                    (sol.u[j][V.pERKn] + sol.u[j][V.ppERKn]) * (p[C.Vn]/p[C.Vc])
                    + sol.u[j][V.pERKc] + sol.u[j][V.ppERKc]
                )
                simulations[observables_index("Phosphorylated_cFos"),j,i] = (
                    sol.u[j][V.pcFOSn]*(p[C.Vn]/p[C.Vc]) + sol.u[j][V.pcFOSc]
                )
            end
        end
    end
end
end # module