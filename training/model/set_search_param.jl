# Specify model parameters and/or initial values to optimize
function get_search_index()::Tuple{Array{Int64,1},Array{Int64,1}}
    # parameters
    search_idx_params::Vector{Int} = [
        C.VmaxPY
        C.KmPY
        C.kdeg
        C.kf47
        C.Vmaxr47
        C.Kmf47
        C.Kmr47
        C.kf48
        C.Kmf48
        C.Kmr48
        #C.PTEN
        C.kf49
        C.kr49
        C.Kmf49
        C.Kmr49
        C.Kmr49b
        C.kr49b
        C.kf51
        C.Vmaxr51
        C.Kmf51
        C.Kmr51
        C.Kmrb51
        C.kf52
        C.Vmaxr52
        C.Kmf52
        C.Kmr52
        C.kf54
        C.Vmaxr54
        C.Kmf54
        C.Kmr54
        C.kf55
        C.Vmaxr55
        C.Kmf55
        C.Kmr55
        C.kf38
        C.kf39
        C.kf50
        C.a98
        C.b98
        C.koff46
        C.EGF_off
        C.HRGoff_3
        C.HRGoff_4
        C.koff4
        C.koff5
        C.koff6
        C.koff7
        C.koff8
        C.koff9
        C.koff61
        C.koff62
        C.koff16
        C.koff17
        C.koff18
        C.koff19
        C.koff20
        C.koff21
        C.koff22
        C.koff23
        C.koff24
        C.koff25
        C.koff26
        C.koff27
        C.koff28
        C.koff29
        C.koff30
        C.koff31
        C.koff32
        C.koff33
        C.koff34
        C.koff35
        C.koff36
        C.koff37
        C.koff65
        C.koff66
        C.koff67
        C.koff68
        C.koff69
        C.koff70
        C.koff71
        C.koff72
        C.koff40
        C.koff41
        C.koff42
        C.koff43
        C.koff44
        C.koff45
        C.koff57
        C.koff58
        C.koff59
        C.koff60
        C.kPTP10
        C.kPTP11
        C.kPTP12
        C.kPTP13
        C.kPTP14
        C.kPTP15
        C.kPTP63
        C.kPTP64
        C.koff73
        C.koff74
        C.koff75
        C.koff76
        C.koff77
        C.koff78
        C.koff79
        C.koff80
        C.kPTP38
        C.kPTP39
        C.koff88
        C.kPTP50
        C.kf81
        C.Vmaxr81
        C.Kmf81
        C.Kmr81
        C.kf82
        C.Vmaxr82
        C.Kmf82
        C.Kmr82
        C.kf83
        C.Vmaxr83
        C.Kmf83
        C.Kmr83
        C.kf84
        C.Vmaxr84
        C.Kmf84
        C.Kmr84
        C.kf85
        C.Vmaxr85
        C.Kmf85
        C.Kmr85
        C.kcon49
        C.kon1
        C.kon86
        C.kon2
        C.kon3
        C.kon87
        C.kon4
        C.kon5
        C.kon6
        C.kon7
        C.kon8
        C.kon9
        C.kon61
        C.kon62
        C.kf10
        C.kf11
        C.kf12
        C.kf13
        C.kf14
        C.kf15
        C.kf63
        C.kf64
        C.kon16
        C.kon17
        C.kon18
        C.kon73
        C.kon19
        C.kon20
        C.kon21
        C.kon74
        C.kon22
        C.kon23
        C.kon24
        C.kon25
        C.kon75
        C.kon26
        C.kon27
        C.kon28
        C.kon29
        C.kon76
        C.kon30
        C.kon31
        C.kon32
        C.kon33
        C.kon77
        C.kon34
        C.kon35
        C.kon36
        C.kon37
        C.kon78
        C.kon65
        C.kon66
        C.kon67
        C.kon68
        C.kon79
        C.kon69
        C.kon70
        C.kon71
        C.kon72
        C.kon80
        C.kon40
        C.kon41
        C.kon42
        C.kon43
        C.kon44
        C.kon45
        C.kon88
        C.kon46
        C.kon57
        C.kon58
        C.kon59
        C.kon60
        #
        C.V1
        C.Km1
        C.V5
        C.Km5
        C.V10
        C.Km10
        #C.n10
        C.p11
        C.p12
        C.p13
        C.V14
        C.Km14
        C.V15
        C.Km15
        C.KimDUSP
        C.KexDUSP
        C.V20
        C.Km20
        C.V21
        C.Km21
        C.V24
        C.Km24
        C.V25
        C.Km25
        C.KimRSK
        C.KexRSK
        C.V27
        C.Km27
        C.V28
        C.Km28
        C.V29
        C.Km29
        C.V30
        C.Km30
        C.V31
        C.Km31
        #C.n31
        C.p32
        C.p33
        C.p34
        C.V35
        C.Km35
        C.V36
        C.Km36
        C.V37
        C.Km37
        C.KimFOS
        C.KexFOS
        C.V42
        C.Km42
        C.V43
        C.Km43
        C.V44
        C.Km44
        C.p47
        C.m47
        C.p48
        C.p49
        C.m49
        C.p50
        C.p51
        C.m51
        C.V57
        C.Km57
        #C.n57
        C.p58
        C.p59
        C.p60
        C.p61
        C.KimF
        C.KexF
        C.p63
        C.KF31
        #C.nF31
        #
        C.scale_E1
        C.scale_E2
        C.scale_E3
        C.scale_E4
        C.scale_G
        #C.scale_S
        C.a_SHC1
        C.a_SHC2
        C.a_SHC3
        C.a_SHC4
        #C.scale_I
        C.a_PIK3CA
        C.a_PIK3CB
        C.a_PIK3CD
        C.a_PIK3CG
        #ver.5
        C.scale_PTEN
        #C.scale_R
        C.a_RASA1
        C.a_RASA2
        C.a_RASA3
        #C.scale_O
        C.a_SOS1
        C.a_SOS2
        C.scale_A
        #C.scale_Akt
        C.a_AKT1
        C.a_AKT2
        #C.scale_RsD
        C.a_HRAS
        C.a_KRAS
        C.a_NRAS
        #C.scale_Raf
        C.a_ARAF
        C.a_BRAF
        C.a_RAF1
        #C.scale_MEK
        C.a_MAP2K1
        C.a_MAP2K2
        C.scale_T
        C.scale_CREBn
        #C.scale_ERKc
        C.a_MAPK1
        C.a_MAPK3
        C.scale_Elk1n
        #C.scale_RSKc
        C.a_RPS6KA1
        C.a_RPS6KA2
        C.a_RPS6KA3
        C.a_RPS6KA6
    ]

    # initial values
    search_idx_initials::Vector{Int} = [
        V.P2,
    ]

    return search_idx_params, search_idx_initials
end


function get_search_region()::Matrix{Float64}
    p::Vector{Float64} = param_values()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = get_search_index()
    search_param::Vector{Float64} = init_search_param(search_idx, p, u0)

    search_rgn::Matrix{Float64} = zeros(2, length(p)+length(u0))

    # Default: 0.1 ~ 10x
    for (i,j) in enumerate(search_idx[1])
        search_rgn[1,j] = search_param[i]*0.1  # lower bound
        search_rgn[2,j] = search_param[i]*10.0  # upper bound
    end

    # Default: 0.5 ~ 2x
    for (i,j) in enumerate(search_idx[2])
        search_rgn[1,j+length(p)] = search_param[i+length(search_idx[1])]*0.5  # lower bound
        search_rgn[2,j+length(p)] = search_param[i+length(search_idx[1])]*2.0  # upper bound
    end


    # search_rgn[:, C.param_name] = [lower_bound, upper_bound]
    # search_rgn[:, V.var_name+length(p)] = [lower_bound, upper_bound]

    #search_rgn[:,C.PTEN] = [1.0,1000.0]
    search_rgn[:,V.P2+length(p)] = [1.0,1000.0]

    search_rgn[:,C.scale_E1] = [0.1, 100.0]
    search_rgn[:,C.scale_E2] = [0.1, 100.0]
    search_rgn[:,C.scale_E3] = [0.1, 100.0]
    search_rgn[:,C.scale_E4] = [0.1, 100.0]
    search_rgn[:,C.scale_G] = [0.1, 100.0]
    #search_rgn[:,C.scale_S] = [0.1, 100.0]
    search_rgn[:,C.a_SHC1] = [0.1, 100.0]
    search_rgn[:,C.a_SHC2] = [0.1, 100.0]
    search_rgn[:,C.a_SHC3] = [0.1, 100.0]
    search_rgn[:,C.a_SHC4] = [0.1, 100.0]
    #search_rgn[:,C.scale_I] = [0.1, 100.0]
    search_rgn[:,C.a_PIK3CA] = [0.1, 100.0]
    search_rgn[:,C.a_PIK3CB] = [0.1, 100.0]
    search_rgn[:,C.a_PIK3CD] = [0.1, 100.0]
    search_rgn[:,C.a_PIK3CG] = [0.1, 100.0]
    # ver.5
    search_rgn[:,C.scale_PTEN] = [0.1, 100.0]
    #search_rgn[:,C.scale_R] = [0.1, 100.0]
    search_rgn[:,C.a_RASA1] = [0.1, 100.0]
    search_rgn[:,C.a_RASA2] = [0.1, 100.0]
    search_rgn[:,C.a_RASA3] = [0.1, 100.0]
    #search_rgn[:,C.scale_O] = [0.1, 100.0]
    search_rgn[:,C.a_SOS1] = [0.1, 100.0]
    search_rgn[:,C.a_SOS2] = [0.1, 100.0]
    search_rgn[:,C.scale_A] = [0.1, 100.0]
    #search_rgn[:,C.scale_Akt] = [0.1, 100.0]
    search_rgn[:,C.a_AKT1] = [0.1, 100.0]
    search_rgn[:,C.a_AKT2] = [0.1, 100.0]
    #search_rgn[:,C.scale_RsD] = [0.1, 100.0]
    search_rgn[:,C.a_HRAS] = [0.1, 100.0]
    search_rgn[:,C.a_KRAS] = [0.1, 100.0]
    search_rgn[:,C.a_NRAS] = [0.1, 100.0]
    #search_rgn[:,C.scale_Raf] = [0.1, 100.0]
    search_rgn[:,C.a_ARAF] = [0.1, 100.0]
    search_rgn[:,C.a_BRAF] = [0.1, 100.0]
    search_rgn[:,C.a_RAF1] = [0.1, 100.0]
    #search_rgn[:,C.scale_MEK] = [0.1, 100.0]
    search_rgn[:,C.a_MAP2K1] = [0.1, 100.0]
    search_rgn[:,C.a_MAP2K2] = [0.1, 100.0]
    search_rgn[:,C.scale_T] = [0.1, 100.0]
    search_rgn[:,C.scale_CREBn] = [0.1, 100.0]
    #search_rgn[:,C.scale_ERKc] = [0.1, 100.0]
    search_rgn[:,C.a_MAPK1] = [0.1, 100.0]
    search_rgn[:,C.a_MAPK3] = [0.1, 100.0]
    search_rgn[:,C.scale_Elk1n] = [0.1, 100.0]
    #search_rgn[:,C.scale_RSKc] = [0.1, 100.0]
    search_rgn[:,C.a_RPS6KA1] = [0.1, 100.0]
    search_rgn[:,C.a_RPS6KA2] = [0.1, 100.0]
    search_rgn[:,C.a_RPS6KA3] = [0.1, 100.0]
    search_rgn[:,C.a_RPS6KA6] = [0.1, 100.0]

    search_rgn = conv_lin2log!(search_rgn, search_idx)

    return search_rgn
end


function update_param(indiv::Vector{Float64})::Tuple{Array{Float64,1},Array{Float64,1}}
    p::Vector{Float64} = param_values()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = get_search_index()

    for (i,j) in enumerate(search_idx[1])
        @inbounds p[j] = indiv[i]
    end
    for (i,j) in enumerate(search_idx[2])
        @inbounds u0[j] = indiv[i+length(search_idx[1])]
    end

    # scaling
    (p, u0) = tpm2ival!(p, u0, "MCF7")
    #
    # constraints --------------------------------------------------------------
    p[C.V6] = p[C.V5]
    p[C.Km6] = p[C.Km5]
    p[C.KimpDUSP] = p[C.KimDUSP]
    p[C.KexpDUSP] = p[C.KexDUSP]
    p[C.KimpcFOS] = p[C.KimFOS]
    p[C.KexpcFOS] = p[C.KexFOS]
    p[C.p52] = p[C.p47]
    p[C.m52] = p[C.m47]
    p[C.p53] = p[C.p48]
    p[C.p54] = p[C.p49]
    p[C.m54] = p[C.m49]
    p[C.p55] = p[C.p50]
    p[C.p56] = p[C.p51]
    p[C.m56] = p[C.m51]
    # --------------------------------------------------------------------------

    return p, u0
end


function decode_gene2val(indiv_gene::Vector{Float64})::Vector{Float64}
    search_rgn::Matrix{Float64} = get_search_region()
    indiv::Vector{Float64} = zeros(length(indiv_gene))

    for i in eachindex(indiv_gene)
        indiv[i] = 10^(
            indiv_gene[i] * (
                search_rgn[2,i] - search_rgn[1,i]
            ) + search_rgn[1,i]
        )
    end

    return round.(indiv,sigdigits=7)
end


function encode_val2gene(indiv::Vector{Float64})
    search_rgn::Matrix{Float64} = get_search_region()
    indiv_gene::Vector{Float64} = zeros(length(indiv))

    for i in eachindex(indiv)
        indiv_gene[i] = (
            log10(indiv[i]) - search_rgn[1,i]
        ) / (
            search_rgn[2,i] - search_rgn[1,i]
        )
    end

    return indiv_gene
end


function encode_bestIndivVal2randGene(j::Int64, best_indiv::Vector{Float64},
                                        p0_bounds::Vector{Float64})::Float64
    search_rgn::Matrix{Float64} = get_search_region()
    rand_gene::Float64 = (
        log10(
            best_indiv[j]*10^(
                rand() * log10(p0_bounds[2]/p0_bounds[1]) + log10(p0_bounds[1])
            )
        ) - search_rgn[1,j]
    ) / (
        search_rgn[2,j] - search_rgn[1,j]
    )

    return rand_gene
end


function init_search_param(search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
                            p::Vector{Float64}, u0::Vector{Float64})::Vector{Float64}
    search_param = zeros(
        length(search_idx[1]) + length(search_idx[2])
    )
    for (i,j) in enumerate(search_idx[1])
        @inbounds search_param[i] = p[j]
    end
    for (i,j) in enumerate(search_idx[2])
        @inbounds search_param[i+length(search_idx[1])] = u0[j]
    end

    if any(x -> x == 0.0, search_param)
        message::String = "search_param must not contain zero."
        for (_, idx) in enumerate(search_idx[1])
            if p[idx] == 0.0
                error(
                    @sprintf(
                        "`C.%s` in search_idx_params: ", C.parameters[idx]
                    ) * message
                )
            end
        end
        for (_, idx) in enumerate(search_idx[2])
            if u0[idx] == 0.0
                error(
                    @sprintf(
                        "`V.%s` in search_idx_initials: ", V.species[idx]
                    ) * message
                )
            end
        end
    end

    return search_param
end


function conv_lin2log!(search_rgn::Matrix{Float64},
                        search_idx::Tuple{Array{Int64,1},Array{Int64,1}}
                        )::Matrix{Float64}
    for i=1:size(search_rgn,2)
        if minimum(search_rgn[:,i]) < 0.0
            message = "search_rgn[lb,ub] must be positive.\n"
            if i <= C.n_parameters
                error(
                    @sprintf(
                        "`C.%s` ", C.parameters[i]
                    ) * message
                )
            else
                error(
                    @sprintf(
                        "`V.%s` ", V.species[i-C.n_parameters]
                    ) * message
                )
            end
        elseif minimum(search_rgn[:,i]) == 0.0 && maximum(search_rgn[:,i]) != 0.0
            message = "lower_bound must be larger than 0.\n"
            if i <= C.n_parameters
                error(
                    @sprintf(
                        "`C.%s` ", C.parameters[i]
                    ) * message
                )
            else
                error(
                    @sprintf(
                        "`V.%s` ", V.species[i-C.n_parameters]
                    ) * message
                )
            end
        elseif search_rgn[2,i] - search_rgn[1,i] < 0.0
            message = "lower_bound < upper_bound\n"
            if i <= C.n_parameters
                error(
                    @sprintf(
                        "`C.%s` ", C.parameters[i]
                    ) * message
                )
            else
                error(
                    @sprintf(
                        "`V.%s` ", V.species[i-C.n_parameters]
                    ) * message
                )
            end
        end
    end

    nonzero_idx::Vector{Int} = []
    for i=1:size(search_rgn,2)
        if search_rgn[:,i] != [0.0,0.0]
            push!(nonzero_idx,i)
        end
    end
    difference::Vector{Int} = collect(
        symdiff(
            Set(nonzero_idx),
            Set(
                append!(
                    search_idx[1], C.n_parameters .+ search_idx[2]
                )
            )
        )
    )
    if length(difference) > 0
        for (i,j) in enumerate(difference)
            if j <= C.n_parameters
                print(
                    @sprintf(
                        "`C.%s`\n", C.parameters[Int(j)]
                    )
                )
            else
                print(
                    @sprintf(
                        "`V.%s`\n", V.species[Int(j)-C.n_parameters]
                    )
                )
            end
        end
        error(
            "Set these search_params in both search_idx and search_rgn."
        )
    end

    search_rgn = search_rgn[:,nonzero_idx]

    return log10.(search_rgn)
end