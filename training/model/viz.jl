module Viz
include("./observable.jl")

using PyPlot

const cm = PyPlot.cm.get_cmap("tab20")
options = [
    Dict(
        "divided_by" => 1.0,
        "xlim" => (),
        "xticks" => nothing,
        "xlabel" => "Time",
        "ylim" => (),
        "yticks" => nothing,
        "ylabel" => replace(replace(observables[i], "__" => "\n"), "_" => " "),
        "exp_data" => true,
        "legend_loc" => nothing,
        "cmap" => [cm.colors[j] for j in 1:20],
        "shape" => ["o", "v", "^", "<", ">", "8", "s",
                    "p", "*", "h", "H", "D", "d", "P", "X"],
        "dont_show" => [],
    ) for i in 1:length(observables)]
# ---
for i in eachindex(observables)
    options[i]["divided_by"] = 60  # sec. -> min.
    options[i]["xlim"] = (-5, 125)
    options[i]["xticks"] = [0, 30, 60, 90, 120]
    options[i]["xlabel"] = "Time (min)"
    options[i]["ylim"] = (-0.1, 1.3)
    options[i]["yticks"] = [0.0, 0.3, 0.6, 0.9, 1.2]
    options[i]["cmap"] = ["darkslateblue", "orangered"]
    options[i]["shape"] = ["D", "s"]
    options[i]["dont_show"] = []
end


function set_rcParams()
    rc("figure",figsize = (4,3))
    rc("font",family = "Arial")
    rc("font",size = 20)
    rc("axes",linewidth = 1.5)
    rc("xtick.major",width = 1.5)
    rc("ytick.major",width = 1.5)
    rc("lines",linewidth = 1.8)
    rc("lines",markersize = 12)
end

end # module