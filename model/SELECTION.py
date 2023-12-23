# FUNCTION FOR EVALUATING SELECTIVE PRESSURES

from sympy import *
import RELATEDNESS as REL
import FITNESS as FIT

def EVALUATE(o_fn, pars):
    PAR_p_f = 1
    PAR_p_m = 1
    PAR_p_ifa = 1
    PAR_p_ima = 1
    PAR_p_gfa = 1
    PAR_p_gfa_ast = PAR_p_gfa
    PAR_p_gma = 1
    PAR_p_gma_ast = PAR_p_gma
    PAR_p_ifb = 1
    PAR_p_imb = 1
    PAR_p_gfb = 1
    PAR_p_gfb_ast = PAR_p_gfb
    PAR_p_gmb = 1
    PAR_p_gmb_ast = PAR_p_gmb
    PAR_mu_ifa = .1
    PAR_mu_ima = .1
    PAR_mu_gfa = .1
    PAR_mu_gfa_ast = PAR_mu_gfa
    PAR_mu_gma = .1
    PAR_mu_gma_ast = PAR_mu_gma
    PAR_mu_ifb = .1
    PAR_mu_imb = .1
    PAR_mu_gfb = .1
    PAR_mu_gfb_ast = PAR_mu_gfb
    PAR_mu_gmb = .1
    PAR_mu_gmb_ast = PAR_mu_gmb
    cost = 1e-5
    with open(o_fn, 'w') as o_file:
        o_file.write(
            'df' + '\t' + 'dm' + '\t' + 'muf' + '\t' + 'mum' + '\t' + # 0 - 3
            'nfa' + '\t' + 'nfb' + '\t' + 'nma' + '\t' + 'nmb' + '\t' + # 4 - 7
            'u' + '\t' + 'm' + '\t' + # 8 - 9
            'sexes' + '\t' + 'patch_type' + '\t' + 'fec_or_mor' + '\t' + # 10 - 12
            'critical_b' + '\t' + 'c' + '\t' + 
            'r_ff' + '\t' + 'r_mf' + '\t' + 'r_fm' + '\t' + 'r_mm' + '\t' + 
            'age' + '\n')
        # D_FITNESS(df, dm, muf, mum, nfa, nfb, nma, nmb, u, m, sexes, patch_type, fec_or_mor)
        for u in range(len(pars)):
            if pars[u][10] == 'F':
                pars_FF = list(pars[u])
                pars_FF[10] = 'FF'
                der_exp_FF = FIT.D_FITNESS(*pars_FF[ : 13])
                pars_MF = list(pars[u])
                pars_MF[10] = 'MF'
                der_exp_MF = FIT.D_FITNESS(*pars_MF[ : 13])
                if pars[u][11] == 'ALPHA':
                    summarized_exp = ((pars[u][4] - 1) * der_exp_FF + pars[u][6] * der_exp_MF) / (pars[u][4] - 1 + pars[u][6])
                else:
                    summarized_exp = ((pars[u][5] - 1) * der_exp_FF + pars[u][7] * der_exp_MF) / (pars[u][5] - 1 + pars[u][7])
            else:
                pars_FM = list(pars[u])
                pars_FM[10] = 'FM'
                der_exp_FM = FIT.D_FITNESS(*pars_FM[ : 13])
                pars_MM = list(pars[u])
                pars_MM[10] = 'MM'
                der_exp_MM = FIT.D_FITNESS(*pars_MM[ : 13])
                if pars[u][11] == 'ALPHA':
                    summarized_exp = ((pars[u][4] - 1) * der_exp_MM + pars[u][6] * der_exp_FM) / (pars[u][4] - 1 + pars[u][6])
                else:
                    summarized_exp = ((pars[u][5] - 1) * der_exp_MM + pars[u][7] * der_exp_FM) / (pars[u][5] - 1 + pars[u][7])
            res_exp = summarized_exp.evalf(subs = {
                symbols('p_ifa') : PAR_p_ifa, symbols('p_ima') : PAR_p_ima,
                symbols('p_gfa') : PAR_p_gfa, symbols('p_gfa_ast') : PAR_p_gfa_ast,
                symbols('p_gma') : PAR_p_gma, symbols('p_gma_ast') : PAR_p_gma_ast,
                symbols('p_f') : PAR_p_f, symbols('p_m') : PAR_p_m,
                symbols('p_ifb') : PAR_p_ifb, symbols('p_imb') : PAR_p_imb,
                symbols('p_gfb') : PAR_p_gfb, symbols('p_gfb_ast') : PAR_p_gfb_ast,
                symbols('p_gmb') : PAR_p_gmb, symbols('p_gmb_ast') : PAR_p_gmb_ast,
                symbols('mu_ifa') : PAR_mu_ifa, symbols('mu_ima') : PAR_mu_ima,
                symbols('mu_gfa') : PAR_mu_gfa, symbols('mu_gfa_ast') : PAR_mu_gfa_ast,
                symbols('mu_gma') : PAR_mu_gma, symbols('mu_gma_ast') : PAR_mu_gma_ast,
                symbols('mu_ifb') : PAR_mu_ifb, symbols('mu_imb') : PAR_mu_imb,
                symbols('mu_gfb') : PAR_mu_gfb, symbols('mu_gfb_ast') : PAR_mu_gfb_ast,
                symbols('mu_gmb') : PAR_mu_gmb, symbols('mu_gmb_ast') : PAR_mu_gmb_ast
                })
            nf_bar = pars[u][8] * pars[u][4] + (1 - pars[u][8]) * pars[u][5]
            nm_bar = pars[u][8] * pars[u][6] + (1 - pars[u][8]) * pars[u][7]
            if pars[u][11] == 'ALPHA':
                # RELATEDNESS(sexes, lifespan, n_f, n_m, d_f_xxx, dm_xxx, mu_f, mu_m, m_xxx)
                d_fa = pars[u][0] * nf_bar / (pars[u][0] * nf_bar + (1 - pars[u][0]) * pars[u][4])
                d_ma = pars[u][1] * nf_bar / (pars[u][1] * nf_bar + (1 - pars[u][1]) * pars[u][4])
                m_alpha = pars[u][9] * pars[u][6] / (pars[u][9] * pars[u][6] + (1 - pars[u][9]) * nm_bar)
                r_ff_alpha = REL.RELATEDNESS('FF', pars[u][13], pars[u][4], pars[u][6], d_fa, d_ma, pars[u][2], pars[u][3], m_alpha)
                r_mf_alpha = REL.RELATEDNESS('MF', pars[u][13], pars[u][4], pars[u][6], d_fa, d_ma, pars[u][2], pars[u][3], m_alpha)
                r_fm_alpha = REL.RELATEDNESS('FM', pars[u][13], pars[u][4], pars[u][6], d_fa, d_ma, pars[u][2], pars[u][3], m_alpha)
                r_mm_alpha = REL.RELATEDNESS('MM', pars[u][13], pars[u][4], pars[u][6], d_fa, d_ma, pars[u][2], pars[u][3], m_alpha)
                for v in range(pars[u][13]):
                    res = res_exp.evalf(subs = {
                        symbols('r_ff_alpha_a') : r_ff_alpha[v], symbols('r_mf_alpha_a') : r_mf_alpha[v],
                        symbols('r_fm_alpha_a') : r_fm_alpha[v], symbols('r_mm_alpha_a') : r_mm_alpha[v],
                        symbols('c') : cost})
                    cri_b = solve(res, dict = False)[0]
                    o_file.write(
                    str(float(pars[u][0])) + '\t' + str(float(pars[u][1])) + '\t' + str(float(pars[u][2])) + '\t' + str(float(pars[u][3])) + '\t' + 
                    str(float(pars[u][4])) + '\t' + str(float(pars[u][5])) + '\t' + str(float(pars[u][6])) + '\t' + str(float(pars[u][7])) + '\t' + 
                    str(float(pars[u][8])) + '\t' + str(float(pars[u][9])) + '\t' + 
                    pars[u][10] + '\t' + pars[u][11] + '\t' + pars[u][12] + '\t' + 
                    str(cri_b) + '\t' + str(cost) + '\t' + 
                    str(r_ff_alpha[v]) + '\t' + str(r_mf_alpha[v]) + '\t' + str(r_fm_alpha[v]) + '\t' + str(r_mm_alpha[v]) + '\t' + 
                    str(v + 1) + '\n')
            else:
                # RELATEDNESS(sexes, lifespan, n_f, n_m, d_f_xxx, dm_xxx, mu_f, mu_m, m_xxx)
                d_fb = pars[u][0] * nf_bar / (pars[u][0] * nf_bar + (1 - pars[u][0]) * pars[u][5])
                d_mb = pars[u][1] * nf_bar / (pars[u][1] * nf_bar + (1 - pars[u][1]) * pars[u][5])
                m_beta = pars[u][9] * pars[u][7] / (pars[u][9] * pars[u][7] + (1 - pars[u][9]) * nm_bar)
                r_ff_beta = REL.RELATEDNESS('FF', pars[u][13], pars[u][5], pars[u][7], d_fb, d_mb, pars[u][2], pars[u][3], m_beta)
                r_mf_beta = REL.RELATEDNESS('MF', pars[u][13], pars[u][5], pars[u][7], d_fb, d_mb, pars[u][2], pars[u][3], m_beta)
                r_fm_beta = REL.RELATEDNESS('FM', pars[u][13], pars[u][5], pars[u][7], d_fb, d_mb, pars[u][2], pars[u][3], m_beta)
                r_mm_beta = REL.RELATEDNESS('MM', pars[u][13], pars[u][5], pars[u][7], d_fb, d_mb, pars[u][2], pars[u][3], m_beta)
                for v in range(pars[u][13]):
                    res = res_exp.evalf(subs = {
                        symbols('r_ff_beta_a') : r_ff_beta[v], symbols('r_mf_beta_a') : r_mf_beta[v],
                        symbols('r_fm_beta_a') : r_fm_beta[v], symbols('r_mm_beta_a') : r_mm_beta[v],
                        symbols('c') : cost})
                    cri_b = solve(res, dict = False)[0]
                    o_file.write(
                    str(float(pars[u][0])) + '\t' + str(float(pars[u][1])) + '\t' + str(float(pars[u][2])) + '\t' + str(float(pars[u][3])) + '\t' + 
                    str(float(pars[u][4])) + '\t' + str(float(pars[u][5])) + '\t' + str(float(pars[u][6])) + '\t' + str(float(pars[u][7])) + '\t' + 
                    str(float(pars[u][8])) + '\t' + str(float(pars[u][9])) + '\t' + 
                    pars[u][10] + '\t' + pars[u][11] + '\t' + pars[u][12] + '\t' + 
                    str(cri_b) + '\t' + str(cost) + '\t' + 
                    str(r_ff_beta[v]) + '\t' + str(r_mf_beta[v]) + '\t' + str(r_fm_beta[v]) + '\t' + str(r_mm_beta[v]) + '\t' + 
                    str(v + 1) + '\n')