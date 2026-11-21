"""
Soil hydraulic functions: Brooks and Corey (1964) retention and conductivity.

These are line-for-line ports of FNCSAT.FOR, FNCPR.FOR, FNCDSDP.FOR, FNCKREL.FOR,
FNCPOR.FOR and FNCEVAP.FOR. The order of every floating-point operation is kept,
so the Python model reproduces the Fortran results to the last printed digit.

Notation (per soil material m):
    poros  porosity n
    alpha  air-entry (bubbling) pressure head, positive (length)
    expn   Brooks-Corey exponent b = 1/lambda (lambda = pore-size index)
    rksat  saturated hydraulic conductivity K_s (length/time)
    stor   specific storage (1/length)
    sres   residual saturation
    se     saturation at the air-entry head; between he = -alpha and 0 the
           retention curve is linear from se to 1, which keeps dS/dh finite
    he     -alpha

Saturation S is the fraction of the pore space filled with water (0 to 1).
Pressure heads are negative in unsaturated soil.
"""
import numpy as np
from numba import njit

# Fortran's "1.0e-11" in FNCKREL.FOR is a single-precision literal; keep its exact value.
KREL_MIN = float(np.float32(1.0e-11))


@njit(cache=True)
def sat_from_head(h, m, expn, sres, alpha, he, se):
    """Saturation from pressure head (FNCSAT)."""
    if h >= 0.0:
        return 1.0
    elif h > he[m]:
        return 1.0 - (1.0 - se[m]) * (h / he[m])
    else:
        hh = -alpha[m] / h
        bb = 1.0 / expn[m]
        return sres[m] + (1.0 - sres[m]) * hh ** bb


@njit(cache=True)
def head_from_sat(s, m, expn, sres, alpha, he, se):
    """Pressure head from saturation (FNCPR)."""
    if s > se[m]:
        return he[m] * (1.0 - s) / (1.0 - se[m])
    else:
        return -alpha[m] * ((s - sres[m]) / (1.0 - sres[m])) ** (-expn[m])


@njit(cache=True)
def dsat_dhead(h, m, alpha, sres, expn, se, he):
    """Derivative dS/dh (FNCDSDP)."""
    if h >= 0.0:
        return 0.0
    elif h > he[m]:
        return (1.0 - se[m]) / (-he[m])
    else:
        lam = 1.0 / expn[m]
        return (1.0 - sres[m]) * lam / alpha[m] * (h / (-alpha[m])) ** (-1.0 - lam)


@njit(cache=True)
def krel(s, m, expn, sres):
    """Relative hydraulic conductivity (FNCKREL), floored at KREL_MIN."""
    k = ((s - sres[m]) / (1.0 - sres[m])) ** (2.0 * expn[m] + 3.0)
    if k < KREL_MIN:
        k = KREL_MIN
    return k


@njit(cache=True)
def porosity(h, m, poros, stor, hinit):
    """Porosity adjusted for elastic storage (FNCPOR)."""
    return poros[m] + stor[m] * (h - hinit)


@njit(cache=True)
def evap_rate(s, m, swilt, sres, zevap, evmax):
    """Soil evaporation per unit depth of the evaporation zone (FNCEVAP, 24-hour
    uniform potential). Potential above S(h_wilt), falling linearly to zero at
    residual saturation."""
    evmaxd = evmax[m] / zevap
    if s < sres[m]:
        return 0.0
    elif s < swilt[m]:
        return evmaxd * ((s - sres[m]) / (swilt[m] - sres[m]))
    else:
        return evmaxd
