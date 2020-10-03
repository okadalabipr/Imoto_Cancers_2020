import pandas as pd
from name2idx import C, V


CCLE = pd.read_csv('../gene_expression/tpm_values.csv', index_col=0)


def get_tpm(gene, cell_line, source):
    return CCLE.loc[gene, cell_line]


def init_E1(x, cell_line, source):
    return x[C.w_E1] * get_tpm("EGFR", cell_line, source)


def init_E2(x, cell_line, source):
    return x[C.w_E2] * get_tpm("ERBB2", cell_line, source)


def init_E3(x, cell_line, source):
    return x[C.w_E3] * get_tpm("ERBB3", cell_line, source)


def init_E4(x, cell_line, source):
    return x[C.w_E4] * get_tpm("ERBB4", cell_line, source)


def init_G(x, cell_line, source):
    return x[C.w_G] * get_tpm("GRB2", cell_line, source)


def init_S(x, cell_line, source):
    return (x[C.w_SHC1] * get_tpm("SHC1", cell_line, source) +
            x[C.w_SHC2] * get_tpm("SHC2", cell_line, source) +
            x[C.w_SHC3] * get_tpm("SHC3", cell_line, source) +
            x[C.w_SHC4] * get_tpm("SHC4", cell_line, source))


def init_I(x, cell_line, source):
    return (x[C.w_PIK3CA] * get_tpm("PIK3CA", cell_line, source) +
            x[C.w_PIK3CB] * get_tpm("PIK3CB", cell_line, source) +
            x[C.w_PIK3CD] * get_tpm("PIK3CD", cell_line, source) +
            x[C.w_PIK3CG] * get_tpm("PIK3CG", cell_line, source))


def init_PTEN(x, cell_line, source):
    return x[C.w_PTEN] * get_tpm("PTEN", cell_line, source)


def init_R(x, cell_line, source):
    return (x[C.w_RASA1] * get_tpm("RASA1", cell_line, source) +
            x[C.w_RASA2] * get_tpm("RASA2", cell_line, source) +
            x[C.w_RASA3] * get_tpm("RASA3", cell_line, source))


def init_O(x, cell_line, source):
    return (x[C.w_SOS1] * get_tpm("SOS1", cell_line, source) +
            x[C.w_SOS2] * get_tpm("SOS2", cell_line, source))


def init_A(x, cell_line, source):
    return x[C.w_A] * get_tpm("GAB1", cell_line, source)


def init_Akt(x, cell_line, source):
    return (x[C.w_AKT1] * get_tpm("AKT1", cell_line, source) +
            x[C.w_AKT2] * get_tpm("AKT2", cell_line, source))


def init_RsD(x, cell_line, source):
    return (x[C.w_HRAS] * get_tpm("HRAS", cell_line, source) +
            x[C.w_KRAS] * get_tpm("KRAS", cell_line, source) +
            x[C.w_NRAS] * get_tpm("NRAS", cell_line, source))


def init_Raf(x, cell_line, source):
    return (x[C.w_ARAF] * get_tpm("ARAF", cell_line, source) +
            x[C.w_BRAF] * get_tpm("BRAF", cell_line, source) +
            x[C.w_RAF1] * get_tpm("RAF1", cell_line, source))



def init_MEK(x, cell_line, source):
    return (x[C.w_MAP2K1] * get_tpm("MAP2K1", cell_line, source) +
            x[C.w_MAP2K2] * get_tpm("MAP2K2", cell_line, source))


def init_T(x, cell_line, source):
    return x[C.w_T] * get_tpm("PTPN1", cell_line, source)


def init_CREBn(x, cell_line, source):
    return x[C.w_CREBn] * get_tpm("CREB1", cell_line, source)


def init_ERKc(x, cell_line, source):
    return (x[C.w_MAPK1] * get_tpm("MAPK1", cell_line, source) +
            x[C.w_MAPK3] * get_tpm("MAPK3", cell_line, source))


def init_Elk1n(x, cell_line, source):
    return x[C.w_Elk1n] * get_tpm("ELK1", cell_line, source)


def init_RSKc(x, cell_line, source):
    return (x[C.w_RPS6KA1] * get_tpm("RPS6KA1", cell_line, source) +
            x[C.w_RPS6KA2] * get_tpm("RPS6KA2", cell_line, source) +
            x[C.w_RPS6KA3] * get_tpm("RPS6KA3", cell_line, source) +
            x[C.w_RPS6KA6] * get_tpm("RPS6KA6", cell_line, source))


def tpm2ival(x, y0, cell_line, source='CCLE'):
    x[C.PTEN] *= init_PTEN(x, cell_line, source)

    y0[V.E1] *= init_E1(x, cell_line, source)
    y0[V.E2] *= init_E2(x, cell_line, source)
    y0[V.E3] *= init_E3(x, cell_line, source)
    y0[V.E4] *= init_E4(x, cell_line, source)
    y0[V.G] *= init_G(x, cell_line, source)
    y0[V.S] *= init_S(x, cell_line, source)
    y0[V.I] *= init_I(x, cell_line, source)
    y0[V.R] *= init_R(x, cell_line, source)
    y0[V.O] *= init_O(x, cell_line, source)
    y0[V.A] *= init_A(x, cell_line, source)
    y0[V.Akt] *= init_Akt(x, cell_line, source)
    y0[V.RsD] *= init_RsD(x, cell_line, source)
    y0[V.Raf] *= init_Raf(x, cell_line, source)
    y0[V.MEK] *= init_MEK(x, cell_line, source)
    y0[V.T] *= init_T(x, cell_line, source)
    y0[V.CREBn] *= init_CREBn(x, cell_line, source)
    y0[V.ERKc] *= init_ERKc(x, cell_line, source)
    y0[V.Elk1n] *= init_Elk1n(x, cell_line, source)
    y0[V.RSKc] *= init_RSKc(x, cell_line, source)

    return x, y0


def mul_ratio2mcf7(x, y0, patientID, source='CCLE'):
    x[C.PTEN] *= init_PTEN(x, patientID, source) / init_PTEN(x, 'MCF7', source)

    y0[V.E1] *= init_E1(x, patientID, source) / init_E1(x, 'MCF7', source)
    y0[V.E2] *= init_E2(x, patientID, source) / init_E2(x, 'MCF7', source)
    y0[V.E3] *= init_E3(x, patientID, source) / init_E3(x, 'MCF7', source)
    y0[V.E4] *= init_E4(x, patientID, source) / init_E4(x, 'MCF7', source)
    y0[V.G] *= init_G(x, patientID, source) / init_G(x, 'MCF7', source)
    y0[V.S] *= init_S(x, patientID, source) / init_S(x, 'MCF7', source)
    y0[V.I] *= init_I(x, patientID, source) / init_I(x, 'MCF7', source)
    y0[V.R] *= init_R(x, patientID, source) / init_R(x, 'MCF7', source)
    y0[V.O] *= init_O(x, patientID, source) / init_O(x, 'MCF7', source)
    y0[V.A] *= init_A(x, patientID, source) / init_A(x, 'MCF7', source)
    y0[V.Akt] *= init_Akt(x, patientID, source) / init_Akt(x, 'MCF7', source)
    y0[V.RsD] *= init_RsD(x, patientID, source) / init_RsD(x, 'MCF7', source)
    y0[V.Raf] *= init_Raf(x, patientID, source) / init_Raf(x, 'MCF7', source)
    y0[V.MEK] *= init_MEK(x, patientID, source) / init_MEK(x, 'MCF7', source)
    y0[V.T] *= init_T(x, patientID, source) / init_T(x, 'MCF7', source)
    y0[V.CREBn] *= init_CREBn(x, patientID, source) / init_CREBn(x, 'MCF7', source)
    y0[V.ERKc] *= init_ERKc(x, patientID, source) / init_ERKc(x, 'MCF7', source)
    y0[V.Elk1n] *= init_Elk1n(x, patientID, source) / init_Elk1n(x, 'MCF7', source)
    y0[V.RSKc] *= init_RSKc(x, patientID, source) / init_RSKc(x, 'MCF7', source)

    return x, y0