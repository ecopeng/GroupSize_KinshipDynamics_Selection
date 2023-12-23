# FUNCTIONS FOR CALCULATING RELATEDNESS

from sympy import *

def RELATEDNESS(sexes, age, n_f, n_m, d_f, d_m, mu_f, mu_m, m):
    G_F, G_M, G_FF, G_FM, G_MM = symbols('G_F, G_M, G_FF, G_FM, G_MM')
    eq_1 = Eq(G_F, (1 - mu_f) * G_F + mu_f * (m * .5 * (1 + G_FM) + (1 - m) * .5))
    eq_2 = Eq(G_M, (1 - mu_m) * G_M + mu_m * (m * .5 * (1 + G_FM) + (1 - m) * .5))
    eq_3 = Eq(G_FF, (
        (1 - mu_f)**2 * G_FF + 
        2 * (1 - mu_f) * mu_f * (1 - d_f) * (
            m * .5 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
            (1 - m) * .5 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF)) + 
        mu_f**2 * (1 - d_f)**2 * (
            m**2 * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + 2 * G_FM + 1 / n_m * G_M + (1 - 1 / n_m) * G_MM) + 
            2 * m * (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
            (1 - m)**2 * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF))))
    eq_4 = Eq(G_FM, (
        (1 - mu_f) * (1 - mu_m) * G_FM + 
        (1 - mu_f) * mu_m * (1 - d_m) * (
            m * .5 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
            (1 - m) * .5 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF)) + 
        mu_f * (1 - mu_m) * (1 - d_f) * (
            m * .5 * (G_FM + 1 / n_m * G_M + (1 - 1 / n_m) * G_MM) + 
            (1 - m) * .5 * G_FM) + 
        mu_f * mu_m * (1 - d_f) * (1 - d_m) * (
            m**2 * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + 2 * G_FM + 1 / n_m * G_M + (1 - 1 / n_m) * G_MM) + 
            2 * m * (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
            (1 - m)**2 * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF))))
    eq_5 = Eq(G_MM, (
        (1 - mu_m)**2 * G_MM + 
        2 * (1 - mu_m) * mu_m * (1 - d_m) * (m * .5 * (1 / n_m * G_M + (1 - 1 / n_m) * G_MM + G_FM) + (1 - m) * .5 * G_FM) + 
        mu_m**2 * (1 - d_m)**2 * (
            m**2 * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + 2 * G_FM + 1 / n_m * G_M + (1 - 1 / n_m) * G_MM) + 
            2 * m * (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
            (1 - m)**2 * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF))))
    relatedness_equ = solve([eq_1, eq_2, eq_3, eq_4, eq_5], dict = True)[0]
    G_F = relatedness_equ[symbols('G_F')]
    G_FF = relatedness_equ[symbols('G_FF')]
    G_FM = relatedness_equ[symbols('G_FM')]
    G_M = relatedness_equ[symbols('G_M')]
    G_MM = relatedness_equ[symbols('G_MM')]
    G_FF_0 = (
        (1 - mu_f) * (1 - d_f) * (
            m * .5 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
            (1 - m) * .5 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF)) + 
        (1 - d_f) * mu_f * (1 - d_f) * (
            m * (
                m * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + 2 * G_FM + 1 / n_m * G_M + (1 - 1 / n_m) * G_MM) + 
                (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM)) + 
            (1 - m) * (
                m * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
                (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF))))
    G_FM_0 = (
        (1 - mu_f) * (1 - d_m) * (
            m * .5 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
            (1 - m) * .5 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF)) + 
        (1 - d_f) * mu_f * (1 - d_m) * (
            m * (
                m * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + 2 * G_FM + 1 / n_m * G_M + (1 - 1 / n_m) * G_MM) + 
                (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM)) + 
            (1 - m) * (
                m * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
                (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF))))
    G_MF_0 = (
        (1 - mu_m) * (1 - d_f) * (
            m * .5 * (1 / n_m * G_M + (1 - 1 / n_m) * G_MM + G_FM) + 
            (1 - m) * .5 * G_FM) + 
        (1 - d_m) * mu_m * (1 - d_f) * (
            m * (
                m * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + 2 * G_FM + 1 / n_m * G_M + (1 - 1 / n_m) * G_MM) + 
                (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM)) + 
            (1 - m) * (
                m * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
                (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF))))
    G_MM_0 = (
        (1 - mu_m) * (1 - d_m) * (
            m * .5 * (1 / n_m * G_M + (1 - 1 / n_m) * G_MM + G_FM) + 
            (1 - m) * .5 * G_FM) + 
        (1 - d_m) * mu_m * (1 - d_m) * (
            m * (
                m * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + 2 * G_FM + 1 / n_m * G_M + (1 - 1 / n_m) * G_MM) + 
                (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM)) + 
            (1 - m) * (
                m * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF + G_FM) + 
                (1 - m) * .25 * (1 / n_f * G_F + (1 - 1 / n_f) * G_FF))))
    if sexes == 'FF':
        res = [G_FF_0 / G_F]
        FFa = G_FF_0
        MFa = G_MF_0
        for i in range(age):
            t_FFa = (1 - mu_f) * FFa + mu_f * (1 - d_f) * (m * .5 * MFa + .5 * (1 / n_f * G_F + (1 - 1 / n_f) * FFa))
            t_MFa = (1 - mu_m) * MFa + mu_m * (1 - d_m) * (m * .5 * MFa + .5 * (1 / n_f * G_F + (1 - 1 / n_f) * FFa))
            res.append(t_FFa / G_F)
            FFa = t_FFa
            MFa = t_MFa
        return res
    elif sexes == 'FM':
        res = [G_FM_0 / G_M]
        MMa = G_MM_0
        FMa = G_FM_0
        for i in range(age):
            t_FMa = (1 - mu_f) * FMa + mu_f * (1 - d_f) * (m * .5 * (1 / n_m * G_M + (1 - 1 / n_m) * MMa) + .5 * FMa)
            t_MMa = (1 - mu_m) * MMa + mu_m * (1 - d_m) * (m * .5 * (1 / n_m * G_M + (1 - 1 / n_m) * MMa) + .5 * FMa)
            res.append(t_FMa / G_M)
            MMa = t_MMa
            FMa = t_FMa
        return res
    elif sexes == 'MF':
        res = [G_MF_0 / G_F]
        MFa = G_MF_0
        FFa = G_FF_0
        for i in range(age):
            t_MFa = (1 - mu_m) * MFa + mu_m * (1 - d_m) * (m * .5 * MFa + .5 * (1 / n_f * G_F + (1 - 1 / n_f) * FFa))
            t_FFa = (1 - mu_f) * FFa + mu_f * (1 - d_f) * (m * .5 * MFa + .5 * (1 / n_f * G_F + (1 - 1 / n_f) * FFa))
            res.append(t_MFa / G_F)
            MFa = t_MFa
            FFa = t_FFa
        return res
    else:
        res = [G_MM_0 / G_M]
        MMa = G_MM_0
        FMa = G_FM_0
        for i in range(age):
            t_MMa = (1 - mu_m) * MMa + mu_m * (1 - d_m) * (m * .5 * (1 / n_m * G_M + (1 - 1 / n_m) * MMa) + .5 * FMa)
            t_FMa = (1 - mu_f) * FMa + mu_f * (1 - d_f) * (m * .5 * (1 / n_m * G_M + (1 - 1 / n_m) * MMa) + .5 * FMa)
            res.append(t_MMa  / G_M)
            MMa = t_MMa
            FMa = t_FMa
        return res

def CALCULATE(fn, lifespan, u, nfa, nma, nfb, nmb, d_f, d_m, m, mu_f, mu_m):
    nf_bar = u * nfa + (1 - u) * nfb
    nm_bar = u * nma + (1 - u) * nmb
    d_fa = d_f * nf_bar / (d_f * nf_bar + (1 - d_f) * nfa)
    d_ma = d_m * nf_bar / (d_m * nf_bar + (1 - d_m) * nfa)
    d_fb = d_f * nf_bar / (d_f * nf_bar + (1 - d_f) * nfb)
    d_mb = d_m * nf_bar / (d_m * nf_bar + (1 - d_m) * nfb)
    m_a = m * nma / (m * nma + (1 - m) * nm_bar)
    m_b = m * nmb / (m * nmb + (1 - m) * nm_bar)
    FF_a = RELATEDNESS('FF', lifespan, nfa, nma, d_fa, d_ma, mu_f, mu_m, m_a)
    MF_a = RELATEDNESS('MF', lifespan, nfa, nma, d_fa, d_ma, mu_f, mu_m, m_a)
    MM_a = RELATEDNESS('MM', lifespan, nfa, nma, d_fa, d_ma, mu_f, mu_m, m_a)
    FM_a = RELATEDNESS('FM', lifespan, nfa, nma, d_fa, d_ma, mu_f, mu_m, m_a)
    FF_b = RELATEDNESS('FF', lifespan, nfb, nmb, d_fb, d_mb, mu_f, mu_m, m_b)
    MF_b = RELATEDNESS('MF', lifespan, nfb, nmb, d_fb, d_mb, mu_f, mu_m, m_b)
    MM_b = RELATEDNESS('MM', lifespan, nfb, nmb, d_fb, d_mb, mu_f, mu_m, m_b)
    FM_b = RELATEDNESS('FM', lifespan, nfb, nmb, d_fb, d_mb, mu_f, mu_m, m_b)
    with open(fn, 'w') as o_file:
        o_file.write(
            'sex' + '\t' + 'lifespan' + '\t' + 'u' + '\t' + 
            'nfa' + '\t' + 'nma' + '\t' + 'nfb' + '\t' + 'nmb' + '\t' + 
            'd_f' + '\t' + 'd_m' + '\t' + 'm' + '\t' + 
            'mu_f' + '\t' + 'mu_m' + '\t' + 
            'relatedness' + '\t' + 'age' + '\t' + 'patch_type' + '\n')
        for x in range(lifespan):
            o_file.write(
            'F' + '\t' + str(lifespan) + '\t' + str(u) + '\t' + 
            str(nfa) + '\t' + str(nma) + '\t' + str(nfb) + '\t' + str(nmb) + '\t' + 
            str(d_f) + '\t' + str(d_m) + '\t' + str(m) + '\t' + 
            str(mu_f) + '\t' + str(mu_m) + '\t' + 
            str(((nfa - 1) * FF_a[x] + nma * MF_a[x]) / (nfa + nma - 1)) + '\t' + str((x + 1)) + '\t' + 'ALPHA' + '\n' + 
            #
            'M' + '\t' + str(lifespan) + '\t' + str(u) + '\t' + 
            str(nfa) + '\t' + str(nma) + '\t' + str(nfb) + '\t' + str(nmb) + '\t' + 
            str(d_f) + '\t' + str(d_m) + '\t' + str(m) + '\t' + 
            str(mu_f) + '\t' + str(mu_m) + '\t' + 
            str((nfa * FM_a[x] + (nma - 1) * MM_a[x]) / (nfa + nma - 1)) + '\t' + str((x + 1)) + '\t' + 'ALPHA' + '\n' + 
            #
            'F' + '\t' + str(lifespan) + '\t' + str(u) + '\t' + 
            str(nfa) + '\t' + str(nma) + '\t' + str(nfb) + '\t' + str(nmb) + '\t' + 
            str(d_f) + '\t' + str(d_m) + '\t' + str(m) + '\t' + 
            str(mu_f) + '\t' + str(mu_m) + '\t' + 
            str(((nfb - 1) * FF_b[x] + nmb * MF_b[x]) / (nfb + nmb - 1)) + '\t' + str((x + 1)) + '\t' + 'BETA' + '\n' + 
            #
            'M' + '\t' + str(lifespan) + '\t' + str(u) + '\t' + 
            str(nfa) + '\t' + str(nma) + '\t' + str(nfb) + '\t' + str(nmb) + '\t' + 
            str(d_f) + '\t' + str(d_m) + '\t' + str(m) + '\t' + 
            str(mu_f) + '\t' + str(mu_m) + '\t' + 
            str((nfb * FM_b[x] + (nmb - 1) * MM_b[x]) / (nfb + nmb - 1)) + '\t' + str((x + 1)) + '\t' + 'BETA' + '\n')