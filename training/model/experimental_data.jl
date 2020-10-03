module Exp
using StatsBase
using Statistics
include("./observable.jl")

const t = [i*60.0 for i in [0, 5, 15, 30, 45, 60, 90, 120]]  # sec.


function norm01(;egf1::Vector{Int}, hrg1::Vector{Int}, egf2::Vector{Int},
                hrg2::Vector{Int}, egf3::Vector{Int}, hrg3::Vector{Int})
    data1::Matrix{Float64} = hcat(egf1, hrg1)
    data2::Matrix{Float64} = hcat(egf2, hrg2)
    data3::Matrix{Float64} = hcat(egf3, hrg3)

    data1 .= data1 ./ maximum(data1)
    data2 .= data2 ./ maximum(data2)
    data3 .= data3 ./ maximum(data3)

    egf_ave::Vector{Float64} = zeros(length(t))
    hrg_ave::Vector{Float64} = zeros(length(t))
    egf_sem::Vector{Float64} = zeros(length(t))
    hrg_sem::Vector{Float64} = zeros(length(t))

    for i in eachindex(t)
        egf_ave[i] = mean([data1[i,1], data2[i,1], data3[i,1]])
        hrg_ave[i] = mean([data1[i,2], data2[i,2], data3[i,2]])
    end

    ave_vec::Vector{Float64} = vcat(egf_ave, hrg_ave)
    ave_min::Float64 = minimum(ave_vec)
    ave_max::Float64 = maximum(ave_vec)

    data1 .= (data1 .- ave_min) ./ (ave_max .- ave_min)
    data2 .= (data2 .- ave_min) ./ (ave_max .- ave_min)
    data3 .= (data3 .- ave_min) ./ (ave_max .- ave_min)

    for i in eachindex(t)
        egf_ave[i] = mean([data1[i,1], data2[i,1], data3[i,1]])
        hrg_ave[i] = mean([data1[i,2], data2[i,2], data3[i,2]])
        egf_sem[i] = std([data1[i,2], data2[i,2], data3[i,2]]) ./ sqrt(3)
        hrg_sem[i] = std([data1[i,2], data2[i,2], data3[i,2]]) ./ sqrt(3)
    end

    return egf_ave, hrg_ave, egf_sem, hrg_sem
end


experiments = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))
error_bars = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))

mcf7_data_pAkt = norm01(
    egf1=[32101, 156970, 90301, 76709, 63640, 52536, 46414, 57329],
    hrg1=[32101, 565508, 551901, 560064, 489678, 408802, 425323, 451502],
    egf2=[11612, 96189, 43622, 43238, 41007, 29902, 19255, 35079],
    hrg2=[11612, 397931, 432609, 417622, 434519, 509919, 361041, 292523],
    egf3=[66038, 208525, 102689, 117308, 125158, 92086, 68587, 78252],
    hrg3=[66038, 563079, 573540, 521062, 447462, 383774, 434807, 409615],
)
bt474_data_pAkt = norm01(
    egf1=[405198, 356865, 321475, 383445, 346872, 328052, 299123, 316633],
    hrg1=[405198, 357121, 419948, 488508, 495214, 443710, 402765, 451831],
    egf2=[432524, 619289, 581376, 481899, 429541, 399922, 376170, 334923],
    hrg2=[432524, 410919, 413878, 390581, 405359, 408471, 373108, 515120],
    egf3=[150446, 435897, 466378, 443105, 415827, 381441, 398841, 413906],
    hrg3=[150446, 556176, 560385, 539165, 589297, 552227, 540005, 539010],
)
mdamb231_data_pAkt = norm01(
    egf1=[86491, 826975, 575400, 354446, 143728, 107326, 88082, 108892],
    hrg1=[86491, 85990, 114224, 94636, 66139, 92359, 89540, 126099],
    egf2=[56004, 816102, 296684, 165282, 92711, 83090, 55002, 50605],
    hrg2=[56004, 50640, 64646, 55814, 70627, 66657, 68994, 46922],
    egf3=[36061, 729961, 218315, 166109, 68746, 53331, 45191, 21852],
    hrg3=[36061, 75029, 174050, 89546, 88262, 99161, 78673, 58412],
)
experiments[observables_index("Phosphorylated_Akt")] = Dict(
    "MCF7_EGF" => mcf7_data_pAkt[1],
    "MCF7_HRG" => mcf7_data_pAkt[2],
    "BT474_EGF" => bt474_data_pAkt[1],
    "BT474_HRG" => bt474_data_pAkt[2],
    "MDAMB231_EGF" => mdamb231_data_pAkt[1],
    "MDAMB231_HRG" => mdamb231_data_pAkt[2],
)
error_bars[observables_index("Phosphorylated_Akt")] = Dict(
    "MCF7_EGF" => mcf7_data_pAkt[3],
    "MCF7_HRG" => mcf7_data_pAkt[4],
    "BT474_EGF" => bt474_data_pAkt[3],
    "BT474_HRG" => bt474_data_pAkt[4],
    "MDAMB231_EGF" => mdamb231_data_pAkt[3],
    "MDAMB231_HRG" => mdamb231_data_pAkt[4],
)

mcf7_data_pERK = norm01(
    egf1=[65481, 446949, 221435, 283171, 265152, 266056, 204912, 188972],
    hrg1=[65481, 698717, 766252, 710005, 693622, 691856, 522173, 334410],
    egf2=[41927, 507623, 169918, 193671, 164088, 145916, 110844, 130362],
    hrg2=[41927, 605118, 699511, 654697, 579863, 490649, 299946, 229297],
    egf3=[118995, 807929, 338665, 267160, 253820, 230200, 157620, 153112],
    hrg3=[118995, 710436, 673318, 615206, 612686, 523198, 390301, 257664],
)
bt474_data_pERK = norm01(
    egf1=[358203, 550378, 633802, 632047, 500267, 394009, 339650, 221411],
    hrg1=[358203, 531893, 703437, 663640, 629213, 612525, 613871, 643056],
    egf2=[355065, 1421536, 1307969, 1101679, 939944, 689539, 507787, 468836],
    hrg2=[355065, 929915, 897601, 924274, 865529, 820386, 673456, 788623],
    egf3=[61593, 631017, 754722, 652440, 575812, 432406, 315961, 259708],
    hrg3=[61593, 480140, 487770, 463604, 438917, 452289, 470624, 531531],
)
mdamb231_data_pERK = norm01(
    egf1=[314472, 504819, 607786, 618492, 475195, 376035, 293988, 324600],
    hrg1=[314472, 156705, 183456, 277862, 141450, 199719, 253331, 407923],
    egf2=[458693, 1001334, 875594, 834259, 782639, 815888, 629576, 539187],
    hrg2=[458693, 322542, 403985, 331734, 263578, 262142, 276371, 313541],
    egf3=[365691, 747932, 937413, 945635, 870059, 706306, 510590, 451927],
    hrg3=[365691, 340876, 428510, 303543, 269653, 195660, 215972, 388050],
)
experiments[observables_index("Phosphorylated_ERK")] = Dict(
    "MCF7_EGF" => mcf7_data_pERK[1],
    "MCF7_HRG" => mcf7_data_pERK[2],
    "BT474_EGF" => bt474_data_pERK[1],
    "BT474_HRG" => bt474_data_pERK[2],
    "MDAMB231_EGF" => mdamb231_data_pERK[1],
    "MDAMB231_HRG" => mdamb231_data_pERK[2],
)
error_bars[observables_index("Phosphorylated_ERK")] = Dict(
    "MCF7_EGF" => mcf7_data_pERK[3],
    "MCF7_HRG" => mcf7_data_pERK[4],
    "BT474_EGF" => bt474_data_pERK[3],
    "BT474_HRG" => bt474_data_pERK[4],
    "MDAMB231_EGF" => mdamb231_data_pERK[3],
    "MDAMB231_HRG" => mdamb231_data_pERK[4],
)

mcf7_data_pcFos = norm01(
    egf1=[43200, 101848, 134050, 187681, 274701, 188891, 186912, 147868],
    hrg1=[43200, 243299, 340259, 537344, 583257, 527613, 551327, 630883],
    egf2=[36344, 99849, 173325, 179897, 207943, 155466, 154118, 138196],
    hrg2=[36344, 139813, 245333, 389460, 402734, 556006, 513591, 432916],
    egf3=[43604, 83374, 108733, 116103, 113879, 95504, 94969, 94662],
    hrg3=[43604, 111136, 343365, 464180, 453578, 440094, 404483, 589354],
)
bt474_data_pcFos = norm01(
    egf1=[30474, 28409, 41129, 149499, 202734, 209115, 228033, 191176],
    hrg1=[30474, 37651, 86813, 196100, 343675, 393495, 381987, 424195],
    egf2=[29037, 67984, 102222, 193780, 284045, 268596, 257846, 173705],
    hrg2=[29037, 54659, 170148, 312924, 332799, 416213, 509015, 479438],
    egf3=[41306, 106705, 223130, 255095, 311870, 301535, 296093, 274372],
    hrg3=[41306, 81420, 176053, 264250, 400465, 424547, 522829, 609155],
)
mdamb231_data_pcFos = norm01(
    egf1=[124104, 175391, 230831, 210449, 337825, 317246, 375420, 311936],
    hrg1=[124104, 124400, 102863, 125643, 124228, 144544, 156551, 163124],
    egf2=[60197, 45929, 63355, 77438, 115070, 122509, 166049, 129227],
    hrg2=[60197, 56254, 57765, 67667, 60445, 66214, 33948, 32766],
    egf3=[47449, 63064, 67735, 183991, 275103, 183211, 552227, 362380],
    hrg3=[47449, 88082, 93539, 73279, 44637, 55030, 56284, 35506],
)
experiments[observables_index("Phosphorylated_cFos")] = Dict(
    "MCF7_EGF" => mcf7_data_pcFos[1],
    "MCF7_HRG" => mcf7_data_pcFos[2],
    "BT474_EGF" => bt474_data_pcFos[1],
    "BT474_HRG" => bt474_data_pcFos[2],
    "MDAMB231_EGF" => mdamb231_data_pcFos[1],
    "MDAMB231_HRG" => mdamb231_data_pcFos[2],
)
error_bars[observables_index("Phosphorylated_cFos")] = Dict(
    "MCF7_EGF" => mcf7_data_pcFos[3],
    "MCF7_HRG" => mcf7_data_pcFos[4],
    "BT474_EGF" => bt474_data_pcFos[3],
    "BT474_HRG" => bt474_data_pcFos[4],
    "MDAMB231_EGF" => mdamb231_data_pcFos[3],
    "MDAMB231_HRG" => mdamb231_data_pcFos[4],
)


function get_timepoint(obs_name::String)::Vector{Float64}
    if obs_name in ["Phosphorylated_ErbB1", "Phosphorylated_ErbB2",
                    "Phosphorylated_ErbB3", "Phosphorylated_ErbB4",
                    "Phosphorylated_Akt", "Phosphorylated_ERK",
                    "Phosphorylated_cFos"]
        return t
    end
end
end # module