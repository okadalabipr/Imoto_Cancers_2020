import numpy as np

from load_csv import tpm2ival
from .name2idx import C, V
from .set_model import param_values, initial_values


class SearchParam(object):
    """ Specify model parameters and/or initial values to optimize
    """
    # parameters
    idx_params = [
        C.VmaxPY,
        C.KmPY,
        C.kdeg,
        C.kf47,
        C.Vmaxr47,
        C.Kmf47,
        C.Kmr47,
        C.kf48,
        C.Kmf48,
        C.Kmr48,
        #C.PTEN,
        C.kf49,
        C.kr49,
        C.Kmf49,
        C.Kmr49,
        C.Kmr49b,
        C.kr49b,
        C.kf51,
        C.Vmaxr51,
        C.Kmf51,
        C.Kmr51,
        C.Kmrb51,
        C.kf52,
        C.Vmaxr52,
        C.Kmf52,
        C.Kmr52,
        C.kf54,
        C.Vmaxr54,
        C.Kmf54,
        C.Kmr54,
        C.kf55,
        C.Vmaxr55,
        C.Kmf55,
        C.Kmr55,
        C.kf38,
        C.kf39,
        C.kf50,
        C.a98,
        C.b98,
        C.koff46,
        C.EGF_off,
        C.HRGoff_3,
        C.HRGoff_4,
        C.koff4,
        C.koff5,
        C.koff6,
        C.koff7,
        C.koff8,
        C.koff9,
        C.koff61,
        C.koff62,
        C.koff16,
        C.koff17,
        C.koff18,
        C.koff19,
        C.koff20,
        C.koff21,
        C.koff22,
        C.koff23,
        C.koff24,
        C.koff25,
        C.koff26,
        C.koff27,
        C.koff28,
        C.koff29,
        C.koff30,
        C.koff31,
        C.koff32,
        C.koff33,
        C.koff34,
        C.koff35,
        C.koff36,
        C.koff37,
        C.koff65,
        C.koff66,
        C.koff67,
        C.koff68,
        C.koff69,
        C.koff70,
        C.koff71,
        C.koff72,
        C.koff40,
        C.koff41,
        C.koff42,
        C.koff43,
        C.koff44,
        C.koff45,
        C.koff57,
        C.koff58,
        C.koff59,
        C.koff60,
        C.kPTP10,
        C.kPTP11,
        C.kPTP12,
        C.kPTP13,
        C.kPTP14,
        C.kPTP15,
        C.kPTP63,
        C.kPTP64,
        C.koff73,
        C.koff74,
        C.koff75,
        C.koff76,
        C.koff77,
        C.koff78,
        C.koff79,
        C.koff80,
        C.kPTP38,
        C.kPTP39,
        C.koff88,
        C.kPTP50,
        C.kf81,
        C.Vmaxr81,
        C.Kmf81,
        C.Kmr81,
        C.kf82,
        C.Vmaxr82,
        C.Kmf82,
        C.Kmr82,
        C.kf83,
        C.Vmaxr83,
        C.Kmf83,
        C.Kmr83,
        C.kf84,
        C.Vmaxr84,
        C.Kmf84,
        C.Kmr84,
        C.kf85,
        C.Vmaxr85,
        C.Kmf85,
        C.Kmr85,
        C.kcon49,
        C.kon1,
        C.kon86,
        C.kon2,
        C.kon3,
        C.kon87,
        C.kon4,
        C.kon5,
        C.kon6,
        C.kon7,
        C.kon8,
        C.kon9,
        C.kon61,
        C.kon62,
        C.kf10,
        C.kf11,
        C.kf12,
        C.kf13,
        C.kf14,
        C.kf15,
        C.kf63,
        C.kf64,
        C.kon16,
        C.kon17,
        C.kon18,
        C.kon73,
        C.kon19,
        C.kon20,
        C.kon21,
        C.kon74,
        C.kon22,
        C.kon23,
        C.kon24,
        C.kon25,
        C.kon75,
        C.kon26,
        C.kon27,
        C.kon28,
        C.kon29,
        C.kon76,
        C.kon30,
        C.kon31,
        C.kon32,
        C.kon33,
        C.kon77,
        C.kon34,
        C.kon35,
        C.kon36,
        C.kon37,
        C.kon78,
        C.kon65,
        C.kon66,
        C.kon67,
        C.kon68,
        C.kon79,
        C.kon69,
        C.kon70,
        C.kon71,
        C.kon72,
        C.kon80,
        C.kon40,
        C.kon41,
        C.kon42,
        C.kon43,
        C.kon44,
        C.kon45,
        C.kon88,
        C.kon46,
        C.kon57,
        C.kon58,
        C.kon59,
        C.kon60,
        # Nakakuki et al., Cell (2010)
        C.V1,
        C.Km1,
        C.V5,
        C.Km5,
        C.V10,
        C.Km10,
        #C.n10,
        C.p11,
        C.p12,
        C.p13,
        C.V14,
        C.Km14,
        C.V15,
        C.Km15,
        C.KimDUSP,
        C.KexDUSP,
        C.V20,
        C.Km20,
        C.V21,
        C.Km21,
        C.V24,
        C.Km24,
        C.V25,
        C.Km25,
        C.KimRSK,
        C.KexRSK,
        C.V27,
        C.Km27,
        C.V28,
        C.Km28,
        C.V29,
        C.Km29,
        C.V30,
        C.Km30,
        C.V31,
        C.Km31,
        #C.n31,
        C.p32,
        C.p33,
        C.p34,
        C.V35,
        C.Km35,
        C.V36,
        C.Km36,
        C.V37,
        C.Km37,
        C.KimFOS,
        C.KexFOS,
        C.V42,
        C.Km42,
        C.V43,
        C.Km43,
        C.V44,
        C.Km44,
        C.p47,
        C.m47,
        C.p48,
        C.p49,
        C.m49,
        C.p50,
        C.p51,
        C.m51,
        C.V57,
        C.Km57,
        #C.n57,
        C.p58,
        C.p59,
        C.p60,
        C.p61,
        C.KimF,
        C.KexF,
        C.p63,
        C.KF31,
        #C.nF31,
        #
        C.w_E1,
        C.w_E2,
        C.w_E3,
        C.w_E4,
        C.w_G,
        #C.w_S,
        C.w_SHC1,
        C.w_SHC2,
        C.w_SHC3,
        C.w_SHC4,
        #C.w_I,
        C.w_PIK3CA,
        C.w_PIK3CB,
        C.w_PIK3CD,
        C.w_PIK3CG,
        # ver. 5
        C.w_PTEN,
        #C.w_R,
        C.w_RASA1,
        C.w_RASA2,
        C.w_RASA3,
        #C.w_O,
        C.w_SOS1,
        C.w_SOS2,
        C.w_A,
        #C.w_Akt,
        C.w_AKT1,
        C.w_AKT2,
        #C.w_RsD,
        C.w_HRAS,
        C.w_KRAS,
        C.w_NRAS,
        #C.w_Raf,
        C.w_ARAF,
        C.w_BRAF,
        C.w_RAF1,
        #C.w_MEK,
        C.w_MAP2K1,
        C.w_MAP2K2,
        C.w_T,
        C.w_CREBn,
        #C.w_ERKc,
        C.w_MAPK1,
        C.w_MAPK3,
        C.w_Elk1n,
        #C.w_RSKc,
        C.w_RPS6KA1,
        C.w_RPS6KA2,
        C.w_RPS6KA3,
        C.w_RPS6KA6,
    ]

    # initial values
    idx_initials = [
        V.P2,
    ]

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_param = self._init_search_param(x, y0)

        search_rgn = np.zeros((2, len(x)+len(y0)))
        # Default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            search_rgn[0, j] = search_param[i] * 0.1  # lower bound
            search_rgn[1, j] = search_param[i] * 10.  # upper bound
        # Default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j+len(x)] = \
                search_param[i+len(self.idx_params)] * 0.5  # lower bound
            search_rgn[1, j+len(x)] = \
                search_param[i+len(self.idx_params)] * 2.0  # upper bound

        # search_rgn[:, C.parameter] = [lower_bound,upper_bound]
        # search_rgn[:, V.specie+len(x)] = [lower_bound,upper_bound]
        search_rgn[:, V.P2+len(x)] = [1.0, 1000.0]

        search_rgn[:, C.w_E1] = [0.1, 100.0]
        search_rgn[:, C.w_E2] = [0.1, 100.0]
        search_rgn[:, C.w_E3] = [0.1, 100.0]
        search_rgn[:, C.w_E4] = [0.1, 100.0]
        search_rgn[:, C.w_G] = [0.1, 100.0]
        #search_rgn[:, C.w_S] = [0.1, 100.0]
        search_rgn[:, C.w_SHC1] = [0.1, 100.0]
        search_rgn[:, C.w_SHC2] = [0.1, 100.0]
        search_rgn[:, C.w_SHC3] = [0.1, 100.0]
        search_rgn[:, C.w_SHC4] = [0.1, 100.0]
        #search_rgn[:, C.w_I] = [0.1, 100.0]
        search_rgn[:, C.w_PIK3CA] = [0.1, 100.0]
        search_rgn[:, C.w_PIK3CB] = [0.1, 100.0]
        search_rgn[:, C.w_PIK3CD] = [0.1, 100.0]
        search_rgn[:, C.w_PIK3CG] = [0.1, 100.0]
        # ver.5
        search_rgn[:, C.w_PTEN] = [0.1, 100.0]
        #search_rgn[:, C.w_R] = [0.1, 100.0]
        search_rgn[:, C.w_RASA1] = [0.1, 100.0]
        search_rgn[:, C.w_RASA2] = [0.1, 100.0]
        search_rgn[:, C.w_RASA3] = [0.1, 100.0]
        #search_rgn[:, C.w_O] = [0.1, 100.0]
        search_rgn[:, C.w_SOS1] = [0.1, 100.0]
        search_rgn[:, C.w_SOS2] = [0.1, 100.0]
        search_rgn[:, C.w_A] = [0.1, 100.0]
        #search_rgn[:, C.w_Akt] = [0.1, 100.0]
        search_rgn[:, C.w_AKT1] = [0.1, 100.0]
        search_rgn[:, C.w_AKT2] = [0.1, 100.0]
        #search_rgn[:, C.w_RsD] = [0.1, 100.0]
        search_rgn[:, C.w_HRAS] = [0.1, 100.0]
        search_rgn[:, C.w_KRAS] = [0.1, 100.0]
        search_rgn[:, C.w_NRAS] = [0.1, 100.0]
        #search_rgn[:, C.w_Raf] = [0.1, 100.0]
        search_rgn[:, C.w_ARAF] = [0.1, 100.0]
        search_rgn[:, C.w_BRAF] = [0.1, 100.0]
        search_rgn[:, C.w_RAF1] = [0.1, 100.0]
        #search_rgn[:, C.w_MEK] = [0.1, 100.0]
        search_rgn[:, C.w_MAP2K1] = [0.1, 100.0]
        search_rgn[:, C.w_MAP2K2] = [0.1, 100.0]
        search_rgn[:, C.w_T] = [0.1, 100.0]
        search_rgn[:, C.w_CREBn] = [0.1, 100.0]
        #search_rgn[:, C.w_ERKc] = [0.1, 100.0]
        search_rgn[:, C.w_MAPK1] = [0.1, 100.0]
        search_rgn[:, C.w_MAPK3] = [0.1, 100.0]
        search_rgn[:, C.w_Elk1n] = [0.1, 100.0]
        #search_rgn[:, C.w_RSKc] = [0.1, 100.0]
        search_rgn[:, C.w_RPS6KA1] = [0.1, 100.0]
        search_rgn[:, C.w_RPS6KA2] = [0.1, 100.0]
        search_rgn[:, C.w_RPS6KA3] = [0.1, 100.0]
        search_rgn[:, C.w_RPS6KA6] = [0.1, 100.0]

        search_rgn = self._conv_lin2log(search_rgn)

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i+len(self.idx_params)]

        x, y0 = tpm2ival(x, y0, 'BT-474')
        # constraints --------------------------------------------------------------
        x[C.V6] = x[C.V5]
        x[C.Km6] = x[C.Km5]
        x[C.KimpDUSP] = x[C.KimDUSP]
        x[C.KexpDUSP] = x[C.KexDUSP]
        x[C.KimpcFOS] = x[C.KimFOS]
        x[C.KexpcFOS] = x[C.KexFOS]
        x[C.p52] = x[C.p47]
        x[C.m52] = x[C.m47]
        x[C.p53] = x[C.p48]
        x[C.p54] = x[C.p49]
        x[C.m54] = x[C.m49]
        x[C.p55] = x[C.p50]
        x[C.p56] = x[C.p51]
        x[C.m56] = x[C.m51]
        # --------------------------------------------------------------------------

        return x, y0

    def gene2val(self, indiv_gene):
        search_rgn = self.get_region()
        indiv = 10**(
            indiv_gene * (
                search_rgn[1, :] - search_rgn[0, :]
            ) + search_rgn[0, :]
        )

        return indiv

    def _init_search_param(self, x, y0):
        """Initialize search_param
        """
        if len(self.idx_params) != len(set(self.idx_params)):
            raise ValueError('Duplicate parameters.')
        elif len(self.idx_initials) != len(set(self.idx_initials)):
            raise ValueError('Duplicate species.')
        search_param = np.empty(
            len(self.idx_params) + len(self.idx_initials)
        )
        for i, j in enumerate(self.idx_params):
            search_param[i] = x[j]
        for i, j in enumerate(self.idx_initials):
            search_param[i+len(self.idx_params)] = y0[j]

        if np.any(search_param == 0.):
            message = 'search_param must not contain zero.'
            for idx in self.idx_params:
                if x[int(idx)] == 0.:
                    raise ValueError(
                        '"C.{}" in idx_params: '.format(
                            C.NAMES[int(idx)]
                        ) + message
                    )
            for idx in self.idx_initials:
                if y0[int(idx)] == 0.:
                    raise ValueError(
                        '"V.{}" in idx_initials: '.format(
                            V.NAMES[int(idx)]
                        ) + message
                    )

        return search_param

    def _conv_lin2log(self, search_rgn):
        """Convert Linear scale to Logarithmic scale
        """
        for i in range(search_rgn.shape[1]):
            if np.min(search_rgn[:, i]) < 0.0:
                message = 'search_rgn[lb,ub] must be positive.'
                if i <= C.NUM:
                    raise ValueError(
                        '"C.{}": '.format(C.NAMES[i]) + message
                    )
                else:
                    raise ValueError(
                        '"V.{}": '.format(V.NAMES[i-C.NUM]) + message
                    )
            elif np.min(search_rgn[:, i]) == 0 and np.max(search_rgn[:, i]) != 0:
                message = 'lower_bound must be larger than 0.'
                if i <= C.NUM:
                    raise ValueError(
                        '"C.{}" '.format(C.NAMES[i]) + message
                    )
                else:
                    raise ValueError(
                        '"V.{}" '.format(V.NAMES[i-C.NUM]) + message
                    )
            elif search_rgn[1, i] - search_rgn[0, i] < 0.0:
                message = 'lower_bound < upper_bound'
                if i <= C.NUM:
                    raise ValueError(
                        '"C.{}" : '.format(C.NAMES[i]) + message
                    )
                else:
                    raise ValueError(
                        '"V.{}" : '.format(V.NAMES[i-C.NUM]) + message
                    )
        difference = list(
            set(
                np.where(
                    np.any(search_rgn != 0., axis=0)
                )[0]
            ) ^ set(
                np.append(
                    self.idx_params, [C.NUM + idx for idx in self.idx_initials]
                )
            )
        )
        if len(difference) > 0:
            message = 'in both search_idx and search_rgn'
            for idx in difference:
                if idx <= C.NUM:
                    raise ValueError(
                        'Set "C.{}" '.format(C.NAMES[int(idx)]) + message
                    )
                else:
                    raise ValueError(
                        'Set "V.{}" '.format(V.NAMES[int(idx-C.NUM)]) + message
                    )
        search_rgn = search_rgn[:, np.any(search_rgn != 0., axis=0)]

        return np.log10(search_rgn)
