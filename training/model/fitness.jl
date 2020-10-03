# Residual Sum of Squares
function compute_objval_rss(sim_data::Vector{Float64}, exp_data::Vector{Float64})::Float64
    error::Float64 = 0.0

    for i in eachindex(exp_data)
        @inbounds error += (sim_data[i] - exp_data[i])^2
    end

    return error
end


# Cosine similarity
function compute_objval_cos(sim_data::Vector{Float64}, exp_data::Vector{Float64})::Float64

    error::Float64 = 1.0 - dot(sim_data,exp_data)/(norm(sim_data)*norm(exp_data))

    return error
end


function conditions_index(condition_name::String)::Int

    return findfirst(isequal(condition_name),Sim.conditions)
end


function diff_sim_and_exp(sim_matrix::Matrix{Float64},exp_dict::Dict{String,Array{Float64,1}},
                            exp_timepoint::Vector{Float64},conditions::Vector{String};
                            sim_norm_max::Float64,exp_norm_max::Float64)
    sim_result::Vector{Float64} = []
    exp_result::Vector{Float64} = []

    for (i,condition) in enumerate(conditions)
        append!(sim_result,sim_matrix[Int.(exp_timepoint.+1),i])
        append!(exp_result,exp_dict[condition])
    end

    return (sim_result./sim_norm_max, exp_result./exp_norm_max)
end


# Define an objective function to be minimized.
function objective(indiv_gene::Vector{Float64})::Float64

    indiv::Vector{Float64} = decode_gene2val(indiv_gene)

    (p,u0) = update_param(indiv)

    if Sim.simulate!(p,u0) isa Nothing
        error::Vector{Float64} = zeros(18)

        sim_norm_max = maximum(
            Sim.simulations[observables_index("Phosphorylated_Akt"),:,[conditions_index("MCF7_EGF"),conditions_index("MCF7_HRG")]]
        )
        exp_t = Exp.get_timepoint("Phosphorylated_Akt")
        error[1] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_Akt"), Int.(exp_t.+1), conditions_index("MCF7_EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_Akt")]["MCF7_EGF"]
        )
        error[2] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_Akt"), Int.(exp_t.+1), conditions_index("MCF7_HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_Akt")]["MCF7_HRG"]
        )

        sim_norm_max = maximum(
            Sim.simulations[observables_index("Phosphorylated_ERK"),:,[conditions_index("MCF7_EGF"),conditions_index("MCF7_HRG")]]
        )
        exp_t = Exp.get_timepoint("Phosphorylated_ERK")
        error[3] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_ERK"), Int.(exp_t.+1), conditions_index("MCF7_EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_ERK")]["MCF7_EGF"]
        )
        error[4] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_ERK"), Int.(exp_t.+1), conditions_index("MCF7_HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_ERK")]["MCF7_HRG"]
        )

        sim_norm_max = maximum(
            Sim.simulations[observables_index("Phosphorylated_cFos"),:,[conditions_index("MCF7_EGF"),conditions_index("MCF7_HRG")]]
        )
        exp_t = Exp.get_timepoint("Phosphorylated_cFos")
        error[5] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_cFos"), Int.(exp_t.+1), conditions_index("MCF7_EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_cFos")]["MCF7_EGF"]
        )
        error[6] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_cFos"), Int.(exp_t.+1), conditions_index("MCF7_HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_cFos")]["MCF7_HRG"]
        )

        sim_norm_max = maximum(
            Sim.simulations[observables_index("Phosphorylated_Akt"),:,[conditions_index("BT474_EGF"),conditions_index("BT474_HRG")]]
        )
        exp_t = Exp.get_timepoint("Phosphorylated_Akt")
        error[7] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_Akt"), Int.(exp_t.+1), conditions_index("BT474_EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_Akt")]["BT474_EGF"]
        )
        error[8] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_Akt"), Int.(exp_t.+1), conditions_index("BT474_HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_Akt")]["BT474_HRG"]
        )

        sim_norm_max = maximum(
            Sim.simulations[observables_index("Phosphorylated_ERK"),:,[conditions_index("BT474_EGF"),conditions_index("BT474_HRG")]]
        )
        exp_t = Exp.get_timepoint("Phosphorylated_ERK")
        error[9] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_ERK"), Int.(exp_t.+1), conditions_index("BT474_EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_ERK")]["BT474_EGF"]
        )
        error[10] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_ERK"), Int.(exp_t.+1), conditions_index("BT474_HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_ERK")]["BT474_HRG"]
        )

        sim_norm_max = maximum(
            Sim.simulations[observables_index("Phosphorylated_cFos"),:,[conditions_index("BT474_EGF"),conditions_index("BT474_HRG")]]
        )
        exp_t = Exp.get_timepoint("Phosphorylated_cFos")
        error[11] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_cFos"), Int.(exp_t.+1), conditions_index("BT474_EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_cFos")]["BT474_EGF"]
        )
        error[12] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_cFos"), Int.(exp_t.+1), conditions_index("BT474_HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_cFos")]["BT474_HRG"]
        )

        sim_norm_max = maximum(
            Sim.simulations[observables_index("Phosphorylated_Akt"),:,[conditions_index("MDAMB231_EGF"),conditions_index("MDAMB231_HRG")]]
        )
        exp_t = Exp.get_timepoint("Phosphorylated_Akt")
        error[13] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_Akt"), Int.(exp_t.+1), conditions_index("MDAMB231_EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_Akt")]["MDAMB231_EGF"]
        )
        error[14] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_Akt"), Int.(exp_t.+1), conditions_index("MDAMB231_HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_Akt")]["MDAMB231_HRG"]
        )

        sim_norm_max = maximum(
            Sim.simulations[observables_index("Phosphorylated_ERK"),:,[conditions_index("MDAMB231_EGF"),conditions_index("MDAMB231_HRG")]]
        )
        exp_t = Exp.get_timepoint("Phosphorylated_ERK")
        error[15] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_ERK"), Int.(exp_t.+1), conditions_index("MDAMB231_EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_ERK")]["MDAMB231_EGF"]
        )
        error[16] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_ERK"), Int.(exp_t.+1), conditions_index("MDAMB231_HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_ERK")]["MDAMB231_HRG"]
        )

        sim_norm_max = maximum(
            Sim.simulations[observables_index("Phosphorylated_cFos"),:,[conditions_index("MDAMB231_EGF"),conditions_index("MDAMB231_HRG")]]
        )
        exp_t = Exp.get_timepoint("Phosphorylated_cFos")
        error[17] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_cFos"), Int.(exp_t.+1), conditions_index("MDAMB231_EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_cFos")]["MDAMB231_EGF"]
        )
        error[18] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_cFos"), Int.(exp_t.+1), conditions_index("MDAMB231_HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_cFos")]["MDAMB231_HRG"]
        )

        return sum(error)
    else
        return NaN
    end
end