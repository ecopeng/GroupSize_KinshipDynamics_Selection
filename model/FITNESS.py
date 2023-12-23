# EXPRESSIONS AND FUNCTIONS FOR EVALUATING FITNESS OUTCOMES

from sympy import *
import REPRODUCTIVE_VALUES as RV

df, dm, muf, mum, nfa, nfb, nma, nmb, u, m = symbols('df, dm, muf, mum, nfa, nfb, nma, nmb, u, m')
p_ifa, p_ima, p_gfa, p_gfa_ast, p_gma, p_gma_ast, p_f, p_m = symbols('p_ifa, p_ima, p_gfa, p_gfa_ast, p_gma, p_gma_ast, p_f, p_m')
p_ifb, p_imb, p_gfb, p_gfb_ast, p_gmb, p_gmb_ast = symbols('p_ifb, p_imb, p_gfb, p_gfb_ast, p_gmb, p_gmb_ast')
mu_ifa, mu_ima, mu_gfa, mu_gfa_ast, mu_gma, mu_gma_ast, muf, mum = symbols('mu_ifa, mu_ima, mu_gfa, mu_gfa_ast, mu_gma, mu_gma_ast, muf, mum')
mu_ifb, mu_imb, mu_gfb, mu_gfb_ast, mu_gmb, mu_gmb_ast = symbols('mu_ifb, mu_imb, mu_gfb, mu_gfb_ast, mu_gmb, mu_gmb_ast')
nfbar = u * nfa + (1 - u) * nfb
nmbar = u * nma + (1 - u) * nmb
#----- FA -----#
s_fafa = (
    (1 - mu_ifa) + 
    (1 - df) * p_ifa * (mu_ifa + (nfa - 1) * mu_gfa_ast) / (2 * (
        (1 - df) * p_ifa + 
        (1 - df) * (nfa - 1) * p_gfa_ast + 
        df * nfbar * p_f)) + 
    df * p_ifa * u * muf * nfa / (2 * (
        (1 - df) * p_f * nfa + 
        df * p_f * nfbar)))
s_mafa = (
    (1 - dm) * p_ifa * mu_gma * nma / (2 * (
        (1 - dm) * p_ifa + 
        (1 - dm) * (nfa - 1) * p_gfa_ast + 
        dm * nfbar * p_f)) + 
    dm * p_ifa * u * mum * nma / (2 * (
        (1 - dm) * p_f * nfa + 
        dm * p_f * nfbar)))
s_fbfa = (
    df * p_ifa * (1 - u) * muf * nfb / (2 * (
        (1 - df) * p_f * nfb + 
        df * p_f * nfbar)))
s_mbfa = (
    dm * p_ifa * (1 - u) * mum * nmb / (2 * (
        (1 - dm) * p_f * nfb + 
        dm * p_f * nfbar)))

#----- FB -----#
s_fafb = (
    df * p_ifb * u * muf * nfa / (2 * (
        (1 - df) * p_f * nfa + 
        df * p_f * nfbar)))
s_mafb = (
    dm * p_ifb * u * mum * nma / (2 * (
        (1 - dm) * p_f * nfa + 
        dm * p_f * nfbar)))
s_fbfb = (
    (1 - mu_ifb) + 
    (1 - df) * p_ifb * (
        mu_ifb + 
        (nfb - 1) * mu_gfb_ast) / (2 * (
            (1 - df) * p_ifb + 
            (1 - df) * (nfb - 1) * p_gfb_ast + 
            df * nfbar * p_f)) + 
        df * p_ifb * (1 - u) * muf * nfb / (2 * (
            (1 - df) * p_f * nfb + 
            df * p_f * nfbar)))
s_mbfb = (
    (1 - dm) * p_ifb * mu_gmb * nmb / (2 * (
        (1 - dm) * p_ifb + 
        (1 - dm) * (nfb - 1) * p_gfb_ast + 
        dm * nfbar * p_f)) + 
    dm * p_ifb * (1 - u) * mum * nmb / (2 * (
        (1 - dm) * p_f * nfb + 
        dm * p_f * nfbar)))

#----- MA -----#
s_fama = (
    (
        (1 - df) * p_gfa * nfa * mu_gfa * nfa / (2 * (
            (1 - df) * p_gfa * nfa + 
            df * p_f * nfbar)) + 
        df * p_gfa * nfa * u * muf * nfa / (2 * (
            (1 - df) * p_f * nfa + 
            df * p_f * nfbar))) * m * p_ima / (m * p_ima + m * (nma - 1) * p_gma_ast + (1 - m) * nmbar * p_m) + 
    (1 - m) * p_ima * u / (
        m * p_m * nma + 
        (1 - m) * p_m * nmbar) * (
        (1 - df) * p_f * nfa * muf * nfa / (2 * (
            (1 - df) * p_f * nfa + 
            df * p_f * nfbar)) + 
        df * p_f * nfa * u * muf * nfa / (2 * (
            (1 - df) * p_f * nfa + 
            df * p_f * nfbar))) + 
        (1 - m) * p_ima * (1 - u) / (
            m * p_m * nmb + 
            (1 - m) * p_m * nmbar) * df * p_f * nfb * u * muf * nfa / (2 * (
                (1 - df) * p_f * nfa + 
                df * p_f * nfbar)))
s_mama = (
    (1 - mu_ima) + 
    (
        (1 - dm) * p_gfa * nfa * (mu_ima + (nma - 1) * mu_gma_ast) / (2 * (
            (1 - dm) * p_gfa * nfa + 
            dm * p_f * nfbar)) + 
        dm * p_gfa * nfa * u * mum * nma / (2 * (
            (1 - dm) * p_f * nfa + 
            dm * p_f * nfbar))) * m * p_ima / (m * p_ima + m * (nma - 1) * p_gma_ast + (1 - m) * nmbar * p_m) + 
    (1 - m) * p_ima * u / (
        m * p_m * nma + 
        (1 - m) * p_m * nmbar) * (
        (1 - dm) * p_f * nfa * mum * nma / (2 * (
            (1 - dm) * p_f * nfa + 
            dm * p_f * nfbar)) + 
        dm * p_f * nfa * u * mum * nma / (2 * (
            (1 - dm) * p_f * nfa + 
            dm * p_f * nfbar))) + 
        (1 - m) * p_ima * (1 - u) / (m * p_m * nmb + (1 - m) * p_m * nmbar) * dm * p_f * nfb * u * mum * nma / (2 * (
            (1 - dm) * p_f * nfa + 
            dm * p_f * nfbar)))
s_fbma = (
    df * p_gfa * nfa * (1 - u) * muf * nfb / (2 * (
        (1 - df) * p_f * nfb + 
        df * p_f * nfbar)) * m * p_ima / (m * p_ima + m * (nma - 1) * p_gma_ast + (1 - m) * nmbar * p_f) + 
    (1 - m) * p_ima * u / (m * p_m * nma + (1 - m) * p_m * nmbar) * df * p_f * nfa * (1 - u) * muf * nfb / (2 * (
        (1 - df) * p_f * nfb + 
        df * p_f * nfbar)) + 
    (1 - m) * p_ima * (1 - u) / (m * p_m * nmb + (1 - m) * p_m * nmbar) * (
        (1 - df) * p_f * nfb * muf * nfb / (2 * (
            (1 - df) * p_f * nfb + 
            df * p_f * nfbar)) + 
        df * p_f * nfb * (1 - u) * muf * nfb / (2 * (
            (1 - df) * p_f * nfb + 
            df * p_f * nfbar))))
s_mbma = (
    dm * p_gfa * nfa * (1 - u) * mum * nmb / (2 * (
        (1 - dm) * p_f * nfb + 
        dm * p_f * nfbar)) * m * p_ima / (m * p_ima + m * (nma - 1) * p_gma_ast + (1 - m) * nmbar * p_m) + 
    (1 - m) * p_ima * u / (m * p_m * nma + (1 - m) * p_m * nmbar) * (dm * p_f * nfa * (1 - u) * mum * nmb) / (2 * (
        (1 - dm) * p_f * nfb + 
        dm * p_f * nfbar)) + 
    (1 - m) * p_ima * (1 - u) / (m * p_m * nmb + (1 - m) * p_m * nmbar) * (
        (1 - dm) * p_f * nfb * mum * nmb / (2 * ((1 - dm) * p_f * nfb + dm * p_f * nfbar)) + 
        dm * p_f * nfb * (1 - u) * mum * nmb / (2 * (
            (1 - dm) * p_f * nfb + 
            dm * p_f * nfbar))))

#----- MB -----#
s_famb = (
    df * p_gfb * nfb * u * muf * nfa / (2 * (
        (1 - df) * p_f * nfa + 
        df * p_f * nfbar)) * m * p_imb / (m * p_imb + m * (nmb - 1) * p_gmb_ast + (1 - m) * nmbar * p_m) + 
    (1 - m) * p_imb * u / (m * p_m * nma + (1 - m) * p_m * nmbar) * (
        (1 - df) * p_f * nfa * muf * nfa / (2 * (
            (1 - df) * p_f * nfa + 
            df * p_f * nfbar)) + 
        df * p_f * nfa * u * muf * nfa / (2 * (
            (1 - df) * p_f * nfa + 
            df * p_f * nfbar))) + 
    (1 - m) * p_imb * (1 - u) / (m * p_m * nmb + (1 - m) * p_m * nmbar) * df * p_f * nfb * u * muf * nfa / (2 * (
        (1 - df) * p_f * nfa + 
        df * p_f * nfbar)))
s_mamb = (
    dm * p_gfb * nfb * u * mum * nma / (2 * (
        (1 - dm) * p_f * nfa + dm * p_f * nfbar)) * m * p_imb / (m * p_imb + m * (nmb - 1) * p_gmb_ast + (1 - m) * nmbar * p_m) + 
    (1 - m) * p_imb * u / (m * p_m * nma + (1 - m) * p_m * nmbar) * (
        (1 - dm) * p_f * nfa * mum * nma / (2 * (
            (1 - dm) * p_f * nfa + 
            dm * p_f * nfbar)) + 
        dm * p_f * nfa * u * mum * nma / (2 * (
            (1 - dm) * p_f * nfa + 
            dm * p_f * nfbar))) + 
    (1 - m) * p_imb * (1 - u) / (m * p_m * nmb + (1 - m) * p_m * nmbar) * dm * p_f * nfb * u * mum * nma / (2 * (
        (1 - dm) * p_f * nfa + 
        dm * p_f * nfbar)))
s_fbmb = (
    (
        (1 - df) * p_gfb * nfb * mu_gfb * nfb / (2 * (
            (1 - df) * p_gfb * nfb + 
            df * p_f * nfbar)) + 
        df * p_gfb * nfb * (1 - u) * muf * nfb / (2 * (
            (1 - df) * p_f * nfb + 
            df * p_f * nfbar))) * m * p_imb / (m * p_imb + m * (nmb - 1) * p_gmb_ast + (1 - m) * nmbar * p_m) + 
    (1 - m) * p_imb * u / (m * p_m * nma + (1 - m) * p_m * nmbar) * df * p_f * nfa * (1 - u) * muf * nfb / (2 * (
        (1 - df) * p_f * nfb + df * p_f * nfbar)) + 
    (1 - m) * p_imb * (1 - u) / (m * p_m * nmb + (1 - m) * p_m * nmbar) * (
        (1 - df) * p_f * nfb * muf * nfb / (2 * (
            (1 - df) * p_f * nfb + 
            df * p_f * nfbar)) + 
        df * p_f * nfb * (1 - u) * muf * nfb / (2 * (
            (1 - df) * p_f * nfb + 
            df * p_f * nfbar))))
s_mbmb = (
    (1 - mu_imb) + 
    (
        (1 - dm) * p_gfb * nfb * (mu_imb + (nmb - 1) * mu_gmb_ast) / (2 * (
            (1 - dm) * p_gfb * nfb + 
            dm * p_f * nfbar)) + 
        dm * p_gfb * nfb * (1 - u) * mum * nmb / (2 * (
            (1 - dm) * p_f * nfb + 
            dm * p_f * nfbar))) * m * p_imb / (m * p_imb + m * (nmb - 1) * p_gmb_ast + (1 - m) * nmbar * p_m) + 
    (1 - m) * p_imb * u / (m * p_m * nma + (1 - m) * p_m * nmbar) * dm * p_f * nfa * (1 - u) * mum * nmb / (2 * (
        (1 - dm) * p_f * nfb + dm * p_f * nfbar)) + 
    (1 - m) * p_imb * (1 - u) / (m * p_m * nmb + (1 - m) * p_m * nmbar) * (
        (1 - dm) * p_f * nfb * mum * nmb / (2 * (
            (1 - dm) * p_f * nfb + 
            dm * p_f * nfbar)) + 
        dm * p_f * nfb * (1 - u) * mum * nmb / (2 * (
            (1 - dm) * p_f * nfb + 
            dm * p_f * nfbar))))
v_fa, v_ma, v_fb, v_mb = symbols('v_fa, v_ma, v_fb, v_mb')
w_fa = v_fa * s_fafa + v_ma * s_mafa + v_fb * s_fbfa + v_mb * s_mbfa
w_fb = v_fa * s_fafb + v_ma * s_mafb + v_fb * s_fbfb + v_mb * s_mbfb
w_ma = v_fa * s_fama + v_ma * s_mama + v_fb * s_fbma + v_mb * s_mbma
w_mb = v_fa * s_famb + v_ma * s_mamb + v_fb * s_fbmb + v_mb * s_mbmb

def W_FA_MA_FB_MB(PAR_v_fa, PAR_v_ma, PAR_v_fb, PAR_v_mb,
    PAR_df, PAR_dm, PAR_muf, PAR_mum, PAR_nfa, PAR_nfb, PAR_nma, PAR_nmb, PAR_u, PAR_m):
    par = [
    (v_fa, PAR_v_fa), (v_ma, PAR_v_ma), (v_fb, PAR_v_fb), (v_mb, PAR_v_mb),
    (df, PAR_df), (dm, PAR_dm), (muf, PAR_muf), (mum, PAR_mum),
    (nfa, PAR_nfa), (nfb, PAR_nfb), (nma, PAR_nma), (nmb, PAR_nmb),
    (u, PAR_u), (m, PAR_m)]
    return [w_fa.subs(par), w_fb.subs(par), w_ma.subs(par), w_mb.subs(par)]

def D_FITNESS(df, dm, muf, mum, nfa, nfb, nma, nmb, u, m, sexes, patch_type, fec_or_mor):
    rv = RV.REPRODUCTIVE_VALUES(df, dm, muf, mum, nfa, nfb, nma, nmb, u, m)
    v_fa = rv.get('v_fa')
    v_ma = rv.get('v_ma')
    v_fb = rv.get('v_fb')
    v_mb = rv.get('v_mb')
    W = W_FA_MA_FB_MB(v_fa, v_ma, v_fb, v_mb, df, dm, muf, mum, nfa, nfb, nma, nmb, u, m)
    w_fa = W[0]
    w_fb = W[1]
    w_ma = W[2]
    w_mb = W[3]
    c, b = symbols('c, b')
    r_ff_alpha_a, r_mf_alpha_a, r_fm_alpha_a, r_mm_alpha_a = symbols('r_ff_alpha_a, r_mf_alpha_a, r_fm_alpha_a, r_mm_alpha_a')
    r_ff_beta_a, r_mf_beta_a, r_fm_beta_a, r_mm_beta_a = symbols('r_ff_beta_a, r_mf_beta_a, r_fm_beta_a, r_mm_beta_a')
    p_ifa, p_ima, p_gfa, p_gfa_ast, p_gma_ast = symbols('p_ifa, p_ima, p_gfa, p_gfa_ast, p_gma_ast')
    p_ifb, p_imb, p_gfb, p_gfb_ast, p_gmb_ast = symbols('p_ifb, p_imb, p_gfb, p_gfb_ast, p_gmb_ast')
    mu_ifa, mu_ima, mu_gfa, mu_gfa_ast, mu_gma, mu_gma_ast = symbols('mu_ifa, mu_ima, mu_gfa, mu_gfa_ast, mu_gma, mu_gma_ast')
    mu_ifb, mu_imb, mu_gfb, mu_gfb_ast, mu_gmb, mu_gmb_ast = symbols('mu_ifb, mu_imb, mu_gfb, mu_gfb_ast, mu_gmb, mu_gmb_ast')
    if sexes == 'FF':
        if patch_type == 'ALPHA':
            if fec_or_mor == 'FECUNDITY':
                F_ff_alpha_a = (
                    -c * diff(w_fa, p_ifa) + b / (nfa - 1) * diff(w_fa, p_gfa_ast) + 
                    (nfa - 1) * r_ff_alpha_a * (
                        b / (nfa - 1) * diff(w_fa, p_ifa) + 
                        1 / (nfa - 1) * ((nfa - 2) * b / (nfa - 1) - c) * diff(w_fa, p_gfa_ast)) + 
                    nma * r_mf_alpha_a * (b - c) / nfa * diff(w_ma, p_gfa))
                return F_ff_alpha_a
            else:
                M_ff_alpha_a = (
                    c * diff(w_fa, mu_ifa) + (-b) / (nfa - 1) * diff(w_fa, mu_gfa_ast) + 
                    (nfa - 1) * r_ff_alpha_a * (
                        (-b) / (nfa - 1) * diff(w_fa, mu_ifa) + 
                        1 / (nfa - 1) * ((nfa - 2) * (-b) / (nfa - 1) + c) * diff(w_fa, mu_gfa_ast)) + 
                    nma * r_mf_alpha_a * (c - b) / nfa * diff(w_ma, mu_gfa))
                return M_ff_alpha_a
        else:
            if fec_or_mor == 'FECUNDITY':
                F_ff_beta_a = (
                    -c * diff(w_fb, p_ifb) + b / (nfb - 1) * diff(w_fb, p_gfb_ast) + 
                    (nfb - 1) * r_ff_beta_a * (
                        b / (nfb - 1) * diff(w_fb, p_ifb) + 
                        1 / (nfb - 1) * ((nfb - 2) * b / (nfb - 1) - c) * diff(w_fb, p_gfb_ast)) + 
                    nmb * r_mf_beta_a * (b - c) / nfb * diff(w_mb, p_gfb))
                return F_ff_beta_a
            else:
                M_ff_beta_a = (
                    c * diff(w_fb, mu_ifb) + (-b) / (nfb - 1) * diff(w_fb, mu_gfb_ast) + 
                    (nfb - 1) * r_ff_beta_a * (
                        (-b) / (nfb - 1) * diff(w_fb, mu_ifb) + 
                        1 / (nfb - 1) * ((nfb - 2) * (-b) / (nfb - 1) + c) * diff(w_fb, mu_gfb_ast)) + 
                    nmb * r_mf_beta_a * (c - b) / nfb * diff(w_mb, mu_gfb))
                return M_ff_beta_a
    elif sexes == 'MF':
        if patch_type == 'ALPHA':
            if fec_or_mor == 'FECUNDITY':
                F_mf_alpha_a = (
                    -c * diff(w_fa, p_ifa) + 
                    (nfa - 1) * r_ff_alpha_a * (-c) / (nfa - 1) * diff(w_fa, p_gfa_ast) + 
                    nma * r_mf_alpha_a * (
                        b / nma * diff(w_ma, p_ima) + 
                        b / nma * diff(w_ma, p_gma_ast) + 
                        (-c) / nfa * diff(w_ma, p_gfa)))
                return F_mf_alpha_a
            else:
                M_mf_alpha_a = (
                    c * diff(w_fa, mu_ifa) + (-b) / nma * diff(w_fa, mu_gma) + 
                    (nfa - 1) * r_ff_alpha_a * (
                        c / (nfa - 1) * diff(w_fa, mu_gfa_ast) + 
                        (-b) / nma * diff(w_fa, mu_gma)) + 
                    nma * r_mf_alpha_a * (
                        (-b) / nma * diff(w_ma, mu_ima) + 
                        (-b) / nma * diff(w_ma, mu_gma_ast) + 
                        c / nfa * diff(w_ma, mu_gfa)))
                return M_mf_alpha_a
        else:
            if fec_or_mor == 'FECUNDITY':
                F_mf_beta_a = (
                    -c * diff(w_fb, p_ifb) + 
                    (nfb - 1) * r_ff_beta_a * (-c) / (nfb - 1) * diff(w_fb, p_gfb_ast) + 
                    nmb * r_mf_beta_a * (
                        b / nmb * diff(w_mb, p_imb) + 
                        b / nmb * diff(w_mb, p_gmb_ast) + 
                        (-c) / nfb * diff(w_mb, p_gfb)))
                return F_mf_beta_a
            else:
                M_mf_beta_a = (
                    c * diff(w_fb, mu_ifb) + (-b) / nmb * diff(w_fb, mu_gmb) + 
                    (nfb - 1) * r_ff_beta_a * (
                        c / (nfb - 1) * diff(w_fb, mu_gfb_ast) + 
                        (-b) / nmb * diff(w_fb, mu_gmb)) + 
                    nmb * r_mf_beta_a * (
                        (-b) / nmb * diff(w_mb, mu_imb) + 
                        (-b) / nmb * diff(w_mb, mu_gmb_ast) + 
                        c / nfb * diff(w_mb, mu_gfb)))
                return M_mf_beta_a
    elif sexes == 'FM':
        if patch_type == 'ALPHA':
            if fec_or_mor == 'FECUNDITY':
                F_fm_alpha_a = (
                    -c * diff(w_ma, p_ima) + b / nfa * diff(w_ma, p_gfa) + 
                    nfa * r_fm_alpha_a * (
                        b / nfa * diff(w_fa, p_ifa) + 
                        b / nfa * diff(w_fa, p_gfa_ast)) + 
                    (nma - 1) * r_mm_alpha_a * (
                        -c / (nma - 1) * diff(w_ma, p_gma_ast) + 
                        b / nfa * diff(w_ma, p_gfa)))
                return F_fm_alpha_a
            else:
                M_fm_alpha_a = (
                    c * diff(w_ma, mu_ima) + (-b) / nfa * diff(w_ma, mu_gfa) + 
                    nfa * r_fm_alpha_a * (
                        (-b) / nfa * diff(w_fa, mu_ifa) + 
                        (-b) / nfa * diff(w_fa, mu_gfa_ast) + 
                        c / nma * diff(w_fa, mu_gma)) + 
                    (nma - 1) * r_mm_alpha_a * (
                        c / (nma - 1) * diff(w_ma, mu_gma_ast) + 
                        (-b) / nfa * diff(w_ma, mu_gfa)))
                return M_fm_alpha_a
        else:
            if fec_or_mor == 'FECUNDITY':
                F_fm_beta_a = (
                    -c * diff(w_mb, p_imb) + b / nfb * diff(w_mb, p_gfb) + 
                    nfb * r_fm_beta_a * (
                        b / nfb * diff(w_fb, p_ifb) + 
                        b / nfb * diff(w_fb, p_gfb_ast)) + 
                    (nmb - 1) * r_mm_beta_a * (
                        -c / (nmb - 1) * diff(w_mb, p_gmb_ast) + 
                        b / nfb * diff(w_mb, p_gfb)))
                return F_fm_beta_a
            else:
                M_fm_beta_a = (
                    c * diff(w_mb, mu_imb) + (-b) / nfb * diff(w_mb, mu_gfb) + 
                    nfb * r_fm_beta_a * (
                        (-b) / nfb * diff(w_fb, mu_ifb) + 
                        (-b) / nfb * diff(w_fb, mu_gfb_ast) + 
                        c / nmb * diff(w_fb, mu_gmb)) + 
                    (nmb - 1) * r_mm_beta_a * (
                        c / (nmb - 1) * diff(w_mb, mu_gmb_ast) + 
                        (-b) / nfb * diff(w_mb, mu_gfb)))
                return M_fm_beta_a
    else:
        if patch_type == 'ALPHA':
            if fec_or_mor == 'FECUNDITY':
                F_mm_alpha_a = (
                    -c * diff(w_ma, p_ima) + b / (nma - 1) * diff(w_ma, p_gma_ast) + 
                    (nma - 1) * r_mm_alpha_a * (
                        b / (nma - 1) * diff(w_ma, p_ima) + 
                        1 / (nma - 1) * ((nma - 2) * b / (nma - 1) - c) * diff(w_ma, p_gma_ast)))
                return F_mm_alpha_a
            else:
                M_mm_alpha_a = (
                    c * diff(w_ma, mu_ima) + (-b) / (nma - 1) * diff(w_ma, mu_gma_ast) + 
                    nfa * r_fm_alpha_a * (c - b) / nma * diff(w_fa, mu_gma) + 
                    (nma - 1) * r_mm_alpha_a * (
                        (-b) / (nma - 1) * diff(w_ma, mu_ima) + 
                        1 / (nma - 1) * ((nma - 2) * (-b) / (nma - 1) + c) * diff(w_ma, mu_gma_ast)))
                return M_mm_alpha_a
        else:
            if fec_or_mor == 'FECUNDITY':
                F_mm_beta_a = (
                    -c * diff(w_mb, p_imb) + b / (nmb - 1) * diff(w_mb, p_gmb_ast) + 
                    (nmb - 1) * r_mm_beta_a * (
                        b / (nmb - 1) * diff(w_mb, p_imb) + 
                        1 / (nmb - 1) * ((nmb - 2) * b / (nmb - 1) - c) * diff(w_mb, p_gmb_ast)))
                return F_mm_beta_a
            else:
                M_mm_beta_a = (
                    c * diff(w_mb, mu_imb) + (-b) / (nmb - 1) * diff(w_mb, mu_gmb_ast) + 
                    nfb * r_fm_beta_a * (c - b) / nmb * diff(w_fb, mu_gmb) + 
                    (nmb - 1) * r_mm_beta_a * (
                        (-b) / (nmb - 1) * diff(w_mb, mu_imb) + 
                        1 / (nmb - 1) * ((nmb - 2) * (-b) / (nmb - 1) + c) * diff(w_mb, mu_gmb_ast)))
                return M_mm_beta_a