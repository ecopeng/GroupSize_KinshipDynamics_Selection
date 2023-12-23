# FUNCTION FOR CALCULATING REPRODUCTIVE VALUES

from sympy import *

def REPRODUCTIVE_VALUES(df, dm, muf, mum, nfa, nfb, nma, nmb, u, m):
    nfbar = u * nfa + (1 - u) * nfb
    nmbar = u * nma + (1 - u) * nmb
    ma = m * nma / (m * nma + (1 - m) * nmbar)
    mb = m * nmb / (m * nmb + (1 - m) * nmbar)
    pifa = (1 - df) * nfa / ((1 - df) * nfa + df * nfbar)
    pima = (1 - dm) * nfa / ((1 - dm) * nfa + dm * nfbar)
    pifb = (1 - df) * nfb / ((1 - df) * nfb + df * nfbar)
    pimb = (1 - dm) * nfb / ((1 - dm) * nfb + dm * nfbar)
    #----- 1st col -----#
    w_11_fafa = (
        (1 - muf) + u * muf * nfa / (2 * u * nfa) * (
            pifa + 
            (1 - pifa) * u * nfa / nfbar))
    w_21_mafa = (
        u * mum * nma / (2 * u * nfa) * (
            pima + 
            (1 - pima) * u * nfa / nfbar))
    w_31_fbfa = (1 - u) * muf * nfb / (2 * u * nfa) * (1 - pifb) * u * nfa / nfbar
    w_41_mbfa = (1 - u) * mum * nmb / (2 * u * nfa) * (1 - pimb) * u * nfa / nfbar
    #----- 2nd col -----#
    w_12_fama = (
        u * muf * nfa / (2 * u * nma) * (
            pifa * (ma + (1 - ma) * u * nma / nmbar) + 
            (1 - pifa) * (
                u * nfa / nfbar * (
                    ma + (1 - ma) * u * nma / nmbar) + 
                (1 - u * nfa / nfbar) * (1 - mb) * u * nma / nmbar)))
    w_22_mama = (
        (1 - mum) + u * mum * nma / (2 * u * nma) * (
            pima * (ma + (1 - ma) * u * nma / nmbar) + 
            (1 - pima) * (
                u * nfa / nfbar * (ma + (1 - ma) * u * nma / nmbar) + 
                (1 - u * nfa / nfbar) * (1 - mb) * u * nma / nmbar)))
    w_32_fbma = (
        (1 - u) * muf * nfb / (2 * u * nma) * (
            pifb * (1 - mb) * u * nma / nmbar + 
            (1 - pifb) * (
                u * nfa / nfbar * (
                    ma + (1 - ma) * u * nma / nmbar) + 
                (1 - u * nfa / nfbar) * (1 - mb) * u * nma / nmbar)))
    w_42_mbma = (
        (1 - u) * mum * nmb / (2 * u * nma) * (
            pimb * (1 - mb) * u * nma / nmbar + 
            (1 - pimb) * (
                u * nfa / nfbar * (
                    ma + (1 - ma) * u * nma / nmbar) + 
                (1 - u * nfa / nfbar) * (1 - mb) * u * nma / nmbar)))
    #----- 3rd col -----#
    w_13_fafb = u * muf * nfa / (2 * (1 - u) * nfb) * (1 - pifa) * (1 - u * nfa / nfbar)
    w_23_mafb = u * mum * nma / (2 * (1 - u) * nfb) * (1 - pima) * (1 - u * nfa / nfbar)
    w_33_fbfb = (
        (1 - muf) + (1 - u) * muf * nfb / (2 * (1 - u) * nfb) * (
            pifb + 
            (1 - pifb) * (1 - u * nfa / nfbar)))
    w_43_mbfb = (
        (1 - u) * mum * nmb / (2 * (1 - u) * nfb) * (
            pimb + 
            (1 - pimb) * (1 - u * nfa / nfbar)))
    #----- 4th col -----#
    w_14_famb = (
        u * muf * nfa / (2 * (1 - u) * nmb) * (
            pifa * (1 - ma) * (1 - u * nma / nmbar) + 
            (1 - pifa) * (
                u * nfa / nfbar * (1 - ma) * (1 - u * nma / nmbar) + 
                (1 - u * nfa / nfbar) * (
                    mb + 
                    (1 - mb) * (1 - u * nma / nmbar)))))
    w_24_mamb = (
        u * mum * nma / (2 * (1 - u) * nmb) * (
            pima * (1 - ma) * (1 - u * nma / nmbar) + 
            (1 - pima) * (
                u * nfa / nfbar * (1 - ma) * (1 - u * nma / nmbar) + 
                (1 - u * nfa / nfbar) * (mb + (1 - mb) * (1 - u * nma / nmbar)))))
    w_34_fbmb = (
        (1 - u) * muf * nfb / (2 * (1 - u) * nmb) * (
            pifb * (
                mb + 
                (1 - mb) * (1 - u * nma / nmbar)) + 
            (1 - pifb) * (
                u * nfa / nfbar * (1 - ma) * (1 - u * nma / nmbar) + 
                (1 - u * nfa / nfbar) * (
                    mb + 
                    (1 - mb) * (1 - u * nma / nmbar)))))
    w_44_mbmb = (
        (1 - mum) + (1 - u) * mum * nmb / (2 * (1 - u) * nmb) * (
            pimb * (
                mb + 
                (1 - mb) * (1 - u * nma / nmbar)) + 
            (1 - pimb) * (
                u * nfa / nfbar * (1 - ma) * (1 - u * nma / nmbar) + 
                (1 - u * nfa / nfbar) * (
                    mb + 
                    (1 - mb) * (1 - u * nma / nmbar)))))
    mat = Matrix([
    [w_11_fafa, w_12_fama, w_13_fafb, w_14_famb],
    [w_21_mafa, w_22_mama, w_23_mafb, w_24_mamb], 
    [w_31_fbfa, w_32_fbma, w_33_fbfb, w_34_fbmb],
    [w_41_mbfa, w_42_mbma, w_43_mbfb, w_44_mbmb]])
    abseval = list(map(abs, [ev[0] for ev in mat.left_eigenvects()]))
    evec = mat.left_eigenvects()[abseval.index(max(abseval))][2][0]
    v_fa = evec[0]
    v_ma = evec[1]
    v_fb = evec[2]
    v_mb = evec[3]
    v_bar = (
        v_fa * u * nfa / (u * (nfa + nma) + (1 - u) * (nfb + nmb)) + 
        v_ma * u * nma / (u * (nfa + nma) + (1 - u) * (nfb + nmb)) + 
        v_fb * (1 - u) * nfb / (u * (nfa + nma) + (1 - u) * (nfb + nmb)) + 
        v_mb * (1 - u) * nmb / (u * (nfa + nma) + (1 - u) * (nfb + nmb)))
    return {
    'v_fa' : v_fa / v_bar, 
    'v_ma' : v_ma / v_bar, 
    'v_fb' : v_fb / v_bar, 
    'v_mb' : v_mb / v_bar}