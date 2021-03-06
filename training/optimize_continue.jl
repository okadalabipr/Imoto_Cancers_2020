include("ParamEstim.jl")
using .ParamEstim

function optimize_continue(nth_param_set::Int64)
    search_rgn::Matrix{Float64} = get_search_region()

    max_generation::Int64 = 3000
    n_population::Int64 = 3*size(search_rgn, 2)
    n_children::Int64 = 50
    n_gene::Int64 = size(search_rgn, 2)
    allowable_error::Float64 = 0.0

    p0_bounds::Vector{Float64} = [0.1, 10.0]  # [lower_bound, upper_bound]

    if !isdir("./fitparam/$nth_param_set")
        mkdir("./fitparam/$nth_param_set")
        
        (best_indiv, best_fitness) = ga_v2(
            nth_param_set,
            max_generation,
            n_population,
            n_children,
            n_gene,
            allowable_error
        )
    else
        (best_indiv, best_fitness) = ga_v2_continue(
            nth_param_set,
            max_generation,
            n_population,
            n_children,
            n_gene,
            allowable_error,
            p0_bounds
        )
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    optimize_continue(parse(Int64,ARGS[1]))
end