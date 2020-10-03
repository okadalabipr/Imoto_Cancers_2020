using DataFrames, CSV


const CCLE = CSV.read("../gene_expression/tpm_values.csv", DataFrame)


function get_tpm(gene::String, cell_line::String)::Float64
    return CCLE[CCLE[!,Symbol("Column1")].==gene,Symbol(cell_line)][1]
end


function init_E1(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_E1] * get_tpm("EGFR", cell_line)
end


function init_E2(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_E2] * get_tpm("ERBB2", cell_line)
end


function init_E3(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_E3] * get_tpm("ERBB3", cell_line)
end


function init_E4(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_E4] * get_tpm("ERBB4", cell_line)
end


function init_G(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_G] * get_tpm("GRB2", cell_line)
end


function init_S(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_SHC1] * get_tpm("SHC1", cell_line) +
            p[C.a_SHC2] * get_tpm("SHC2", cell_line) +
            p[C.a_SHC3] * get_tpm("SHC3", cell_line) +
            p[C.a_SHC4] * get_tpm("SHC4", cell_line))
end


function init_I(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_PIK3CA] * get_tpm("PIK3CA", cell_line) +
            p[C.a_PIK3CB] * get_tpm("PIK3CB", cell_line) +
            p[C.a_PIK3CD] * get_tpm("PIK3CD", cell_line) +
            p[C.a_PIK3CG] * get_tpm("PIK3CG", cell_line))
end

function init_PTEN(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_PTEN] * get_tpm("PTEN", cell_line)
end

function init_R(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_RASA1] * get_tpm("RASA1", cell_line) +
            p[C.a_RASA2] * get_tpm("RASA2", cell_line) +
            p[C.a_RASA3] * get_tpm("RASA3", cell_line))
end


function init_O(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_SOS1] * get_tpm("SOS1", cell_line) +
            p[C.a_SOS2] * get_tpm("SOS2", cell_line))
end


function init_A(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_A] * get_tpm("GAB1", cell_line)
end


function init_Akt(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_AKT1] * get_tpm("AKT1", cell_line) +
            p[C.a_AKT2] * get_tpm("AKT2", cell_line))
end


function init_RsD(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_HRAS] * get_tpm("HRAS", cell_line) +
            p[C.a_KRAS] * get_tpm("KRAS", cell_line) +
            p[C.a_NRAS] * get_tpm("NRAS", cell_line))
end


function init_Raf(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_ARAF] * get_tpm("ARAF", cell_line) +
            p[C.a_BRAF] * get_tpm("BRAF", cell_line) +
            p[C.a_RAF1] * get_tpm("RAF1", cell_line))
end


function init_MEK(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_MAP2K1] * get_tpm("MAP2K1", cell_line) +
            p[C.a_MAP2K2] * get_tpm("MAP2K2", cell_line))
end


function init_T(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_T] * get_tpm("PTPN1", cell_line)
end


function init_CREBn(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_CREBn] * get_tpm("CREB1", cell_line)
end


function init_ERKc(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_MAPK1] * get_tpm("MAPK1", cell_line) +
            p[C.a_MAPK3] * get_tpm("MAPK3", cell_line))
end


function init_Elk1n(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.scale_Elk1n] * get_tpm("ELK1", cell_line)
end


function init_RSKc(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.a_RPS6KA1] * get_tpm("RPS6KA1", cell_line) +
            p[C.a_RPS6KA2] * get_tpm("RPS6KA2", cell_line) +
            p[C.a_RPS6KA3] * get_tpm("RPS6KA3", cell_line) +
            p[C.a_RPS6KA6] * get_tpm("RPS6KA6", cell_line))
end


function tpm2ival!(p::Vector{Float64}, u0::Vector{Float64}, cell_line::String)
    p[C.PTEN] *= init_PTEN(p, cell_line)

    u0[V.E1] *= init_E1(p, cell_line)
    u0[V.E2] *= init_E2(p, cell_line)
    u0[V.E3] *= init_E3(p, cell_line)
    u0[V.E4] *= init_E4(p, cell_line)
    u0[V.G] *= init_G(p, cell_line)
    u0[V.S] *= init_S(p, cell_line)
    u0[V.I] *= init_I(p, cell_line)
    u0[V.R] *= init_R(p, cell_line)
    u0[V.O] *= init_O(p, cell_line)
    u0[V.A] *= init_A(p, cell_line)
    u0[V.Akt] *= init_Akt(p, cell_line)
    u0[V.RsD] *= init_RsD(p, cell_line)
    u0[V.Raf] *= init_Raf(p, cell_line)
    u0[V.MEK] *= init_MEK(p, cell_line)
    u0[V.T] *= init_T(p, cell_line)
    u0[V.CREBn] *= init_CREBn(p, cell_line)
    u0[V.ERKc] *= init_ERKc(p, cell_line)
    u0[V.Elk1n] *= init_Elk1n(p, cell_line)
    u0[V.RSKc] *= init_RSKc(p, cell_line)

    return p, u0
end


function mul_ratio2mcf7!(p::Vector{Float64}, u0::Vector{Float64}, cell_line_2::String)
    p[C.PTEN] *= init_PTEN(p, cell_line_2) / init_PTEN(p, "MCF7")

    u0[V.E1] *= init_E1(p, cell_line_2) / init_E1(p, "MCF7")
    u0[V.E2] *= init_E2(p, cell_line_2) / init_E2(p, "MCF7")
    u0[V.E3] *= init_E3(p, cell_line_2) / init_E3(p, "MCF7")
    u0[V.E4] *= init_E4(p, cell_line_2) / init_E4(p, "MCF7")
    u0[V.G] *= init_G(p, cell_line_2) / init_G(p, "MCF7")
    u0[V.S] *= init_S(p, cell_line_2) / init_S(p, "MCF7")
    u0[V.I] *= init_I(p, cell_line_2) / init_I(p, "MCF7")
    u0[V.R] *= init_R(p, cell_line_2) / init_R(p, "MCF7")
    u0[V.O] *= init_O(p, cell_line_2) / init_O(p, "MCF7")
    u0[V.A] *= init_A(p, cell_line_2) / init_A(p, "MCF7")
    u0[V.Akt] *= init_Akt(p, cell_line_2) / init_Akt(p, "MCF7")
    u0[V.RsD] *= init_RsD(p, cell_line_2) / init_RsD(p, "MCF7")
    u0[V.Raf] *= init_Raf(p, cell_line_2) / init_Raf(p, "MCF7")
    u0[V.MEK] *= init_MEK(p, cell_line_2) / init_MEK(p, "MCF7")
    u0[V.T] *= init_T(p, cell_line_2) / init_T(p, "MCF7")
    u0[V.CREBn] *= init_CREBn(p, cell_line_2) / init_CREBn(p, "MCF7")
    u0[V.ERKc] *= init_ERKc(p, cell_line_2) / init_ERKc(p, "MCF7")
    u0[V.Elk1n] *= init_Elk1n(p, cell_line_2) / init_Elk1n(p, "MCF7")
    u0[V.RSKc] *= init_RSKc(p, cell_line_2) / init_RSKc(p, "MCF7")

    return p, u0
end