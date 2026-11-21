"""
Time-stepping kernel for one soil column.

This is a port of PLANT1D.FOR and its subroutines (BCCHK, TIMCHK, ITCHK, TRANSPIR,
MATRIX, THOMAS, CNVCHK, MASBAL, PONDING, DAYAVG, PPRINT, PPRINT2, PPRINT3).
The control flow (including the GOTO-based retry of failed time steps) and the
order of floating-point operations follow the Fortran, so that with one plant and
``legacy_convergence=True`` the kernel reproduces the Fortran output files exactly.

Generalizations, all of which reduce to the Fortran when there is one variety with
cover 1 and evaporation scale 1:
  * several plant varieties, each with its own root profile, plant water potential,
    resistances and potential transpiration (Guswa et al. root uptake per variety);
  * a daily canopy-cover multiplier for each variety and a daily multiplier on soil
    evaporation (e.g. the bare-soil fraction);
  * a free-drainage (unit-gradient) bottom boundary, nbcb = 3;
  * a proper floating-point convergence test (the Fortran truncates the head
    increment to an integer, so it accepts any iterate whose increments are < 1 length unit);
  * runoff bookkeeping when ponded water is discarded (ipond = 1).

Units are whatever the inputs use (the Fortran example is cm and days).
"""
import numpy as np
from numba import njit

from .soil import sat_from_head, dsat_dhead, krel, porosity, evap_rate

ERRTPOT = float(np.float32(1.0e-8))  # single-precision literal in PPRINT.FOR

STATUS_OK = 0          # reached tmax
STATUS_NTMAX = 1       # exceeded the maximum number of time steps
STATUS_NOCONV = 2      # time step fell below dtmin without convergence


@njit(cache=True)
def _grow(a, n_new):
    b = np.zeros((n_new,) + a.shape[1:], dtype=a.dtype)
    b[: a.shape[0]] = a
    return b


@njit(cache=True)
def _interface_k(i, j, dz, imatl, rksat, sat, pr, expn, sres):
    """Harmonic-mean K_s and upwinded k_rel between block i (above) and j (below).
    Returns (K_s * k_rel, half-distance dz12). Same expressions as MATRIX.FOR."""
    xksat = rksat[imatl[i]] * rksat[imatl[j]] * (dz[i] + dz[j]) \
        / (rksat[imatl[j]] * dz[i] + rksat[imatl[i]] * dz[j])
    dz12 = (dz[i] + dz[j]) / 2.0
    delh = (pr[j] - pr[i]) / dz12 - 1.0
    if delh <= -1.0:
        coefkr = 1.0
    elif delh >= 1.0:
        coefkr = 0.0
    else:
        coefkr = 0.5 - delh / 2.0
    xkrel = krel(sat[i], imatl[i], expn, sres) * coefkr + \
        krel(sat[j], imatl[j], expn, sres) * (1.0 - coefkr)
    return xksat, xkrel, dz12


@njit(cache=True)
def _uptake_variety(v, c, pr, sat, dz, imatl, expn, sres, ksu,
                    rw, etm, nroot_v, csoil_v, croot_v, rtavg_v, hwilt_v, tpot_v, out):
    """Root water uptake of variety v (TRANSPIR.FOR). Writes per-block uptake into
    out[v, :] and returns the total. Uptake in block j is

        c * w_j dz_j m_j (h_j - h_p) / ( c_s / (rho K_s k_rel(S_j)) + c_r / rho )

    with plant head h_p = h_wilt unless that would exceed c * T_pot, in which case h_p
    is found by the secant method so that total uptake equals c * T_pot."""
    nr = nroot_v[v]
    for j in range(out.shape[1]):
        out[v, j] = 0.0
    if c <= 0.0:
        return 0.0
    tpot = tpot_v[v] * c
    hp = hwilt_v[v]
    csoil = csoil_v[v]
    croot = croot_v[v]
    rtavg = rtavg_v[v]
    sumtrans = 0.0
    for j in range(nr):
        tj = rw[v, j] * dz[j] * etm[v, j] * ((pr[j] - hp) /
             ((csoil / (rtavg * krel(sat[j], imatl[j], expn, sres) * ksu[j])) + (croot / rtavg)))
        tj = c * tj
        if tj < 0.0:
            tj = 0.0
        out[v, j] = tj
        sumtrans = sumtrans + tj
    if sumtrans <= tpot:
        return sumtrans
    ttol = tpot / 20000.0
    hp_lag = 0.0
    f_lag = -sumtrans
    f = sumtrans - tpot
    for it in range(20):
        dhp = (hp_lag - hp) * (f / (f - f_lag))
        hp_lag = hp
        f_lag = f
        hp = hp + dhp
        sumtrans = 0.0
        for jj in range(nr):
            tj = rw[v, jj] * dz[jj] * etm[v, jj] * ((pr[jj] - hp) /
                 ((csoil / (rtavg * krel(sat[jj], imatl[jj], expn, sres) * ksu[jj])) + (croot / rtavg)))
            tj = c * tj
            if tj < 0.0:
                tj = 0.0
            out[v, jj] = tj
            sumtrans = sumtrans + tj
        f = sumtrans - tpot
        if abs(f) < ttol:
            return sumtrans
    return sumtrans


@njit(cache=True, nogil=True)
def run_column(
    # grid and soil
    dz, imatl, evmult, eweight, zevap, ztop,
    poros, alpha, expn, rksat, stor, sres, se, he, swilt, evmax,
    # plants
    rw, etm, nroot_v, csoil_v, croot_v, rtavg_v, hwilt_v, tpot_v, ksu,
    cover, escale, nleak,
    # boundary conditions (1-based Fortran semantics; arrays padded with a trailing zero)
    nbctop, bctop, timet, nbcbot, bcbot, timeb,
    # initial condition
    pr_init, sat_init, prinit,
    # control
    tstart, dt0, tmax, ntmax, tday, ipond, dtprint, mprint,
    itmin, itmax, dtredu, dtincr, dtmin, dtmax, errpr, errres,
    legacy_convergence, ndepths, rtol=0.0, save_prof=True, exact_rain=False,
):
    nb = dz.shape[0]
    nv = rw.shape[0]
    nnd = ndepths.shape[0]
    ncov = cover.shape[0]
    nesc = escale.shape[0]

    pr = pr_init.copy()
    sat = sat_init.copy()
    prold = pr.copy()
    satold = sat.copy()
    princr = np.zeros(nb)
    trans = np.zeros(nb)
    trans_v = np.zeros((nv, nb))
    tsum_v = np.zeros(nv)
    uptake_prof = np.zeros((nv, nb))     # total uptake by variety and block over the run
    aa = np.zeros((nb, 3))
    rhs = np.zeros(nb)

    # ---------------- output buffers ----------------
    n_print_cap = int((tmax - tstart) / dtprint * 1.02) + 20
    n_day_cap = n_print_cap // max(mprint, 1) + 20
    # print records (PPRINT): t, transum, evsum, qleak, dt, fluxin, flxout, fluxet, fluxev, ratio, ratio1, sat_n...
    prec = np.zeros((n_print_cap, 11 + nnd))
    n_prec = 0
    # day records (PPRINT3): t, dqin, dtransum, devsum, dleak, dqout, drunoff, storage, sumtot, daysat_n..., per variety (T, Tpot, Sroot)
    nd_cols = 9 + nnd + 3 * nv
    drec = np.zeros((n_day_cap, nd_cols))
    n_drec = 0
    # profile records (PPRINT2): t, then por, sat, pr per block
    n_prof_cap = n_day_cap if save_prof else 1
    prof_t = np.zeros(n_prof_cap)
    prof = np.zeros((n_prof_cap, 3, nb))
    n_prof = 0
    # ponding records (fort.2 / .anm): t, qin, dt, qbuild, qinflx, bctop(nprtop)
    pond = np.zeros((256, 6))
    n_pond = 0
    pond_begin = np.zeros(64)
    n_pond_begin = 0

    # ---------------- INITMVS ----------------
    itime = 0
    ipr = 0
    ipr2 = 0
    ipr3 = 0
    iprad = 1
    t = tstart
    dt = dt0
    dtold = dt
    tprint = tstart + dtprint
    nprtop = 1
    nprbot = 1
    nbct = nbctop[0]
    bct = bctop[0]
    nbcb = nbcbot[0]
    bcb = bcbot[0]
    ibcflag = 0
    iqflag = 0
    qbuild = 0.0
    fluxin = 0.0
    flxout = 0.0
    fluxet = 0.0
    fluxev = 0.0
    sumtot = 0.0
    sumtotr = 0.0
    dqin = 0.0
    dtransum = 0.0
    devsum = 0.0
    dleak = 0.0
    dqout = 0.0
    drunoff = 0.0
    daysat = np.zeros(nnd)
    dtv = np.zeros(nv)
    dtpv = np.zeros(nv)
    dsr = np.zeros(nv)
    nprtad = 0
    nprbad = 0
    status = STATUS_OK

    # root-zone weights for each variety's root-weighted saturation
    wroot = np.zeros((nv, nb))
    for v in range(nv):
        ws = 0.0
        for j in range(nroot_v[v]):
            wroot[v, j] = rw[v, j] * etm[v, j] * dz[j]
            ws += wroot[v, j]
        if ws > 0.0:
            for j in range(nroot_v[v]):
                wroot[v, j] /= ws

    while True:
        # ======== label 1: new time step ========
        if (t + dt) >= (tprint - dt / 10.0):
            dtold = dt
            dt = tprint - t
            t = tprint
            ipr = 1
            iprad = 1
            ipr2 = ipr2 + iprad
            ipr3 = ipr3 + iprad
            tprint = t + dtprint
        elif (t + dt) >= (timet[nprtop - 1] - dt / 10.0):
            dtold = dt
            dt = timet[nprtop - 1] - t
            if dt < dtmin:
                dt = dtold
            t = timet[nprtop - 1]
            iprad = 0
        else:
            t = t + dt
            iprad = 0
        itime = itime + 1

        # ---- BCCHK ----
        nprtad = 0
        nprbad = 0
        tcheck = abs(t - timet[nprtop - 1])
        if (tcheck < dtmin) and (ibcflag == 0):
            nprtad = 1
            nprtop = nprtop + nprtad
            ibcflag = 1
            if iqflag == 0:
                nbct = nbctop[nprtop - ibcflag - 1]
                bct = bctop[nprtop - ibcflag - 1]
            if nbctop[nprtop - 1] == 1:
                iqflag = 0
                nbct = nbctop[nprtop - ibcflag - 1]
                bct = bctop[nprtop - ibcflag - 1]
        elif (ibcflag == 1) and (tcheck > dtmin):
            nprtad = 0
            ibcflag = 0
            if iqflag == 0:
                nbct = nbctop[nprtop - 1]
                bct = bctop[nprtop - 1]
            if nbctop[nprtop - 1] == 1:
                iqflag = 0
                nbct = nbctop[nprtop - 1]
                bct = bctop[nprtop - 1]
        elif (ibcflag == 1) and (tcheck < dtmin):
            ibcflag = 1
            nprtad = 1
            nprtop = nprtop + nprtad
            if iqflag == 0:
                nbct = nbctop[nprtop - ibcflag - 1]
                bct = bctop[nprtop - ibcflag - 1]
            if nbctop[nprtop - 1] == 1:
                iqflag = 0
                nbct = nbctop[nprtop - ibcflag - 1]
                bct = bctop[nprtop - ibcflag - 1]
        if t > timeb[nprbot - 1]:
            nprbad = 0
            nprbot = nprbot + nprbad
            nbcb = nbcbot[nprbot - 1]
            bcb = bcbot[nprbot - 1]

        # rain rate actually applied during this step (the interval that ends at or after t)
        kq = nprtop - 1 - ibcflag
        qrain = bctop[kq] if nbctop[kq] == 2 else 0.0

        # ---- TIMCHK ----
        if itime > ntmax:
            status = STATUS_NTMAX
            break
        if t > tmax:
            status = STATUS_OK
            break

        # daily forcing index (cover and evaporation multipliers)
        kday = int((t - 0.5 * dt - tstart) / tday)
        kc = min(max(kday, 0), ncov - 1)
        ke = min(max(kday, 0), nesc - 1)
        esc = escale[ke]

        # ======== label 2: Picard iterations ========
        it = 0
        itest = 0
        failed = False
        while True:
            it = it + 1
            rmax = 0.0
            # ---- ITCHK ----
            if it > itmax:
                t = t - dt
                dt = dt * dtredu
                if (dt < dtmin) or (dtredu >= 1.0):
                    failed = True
                    break
                ipr = 0
                ipr2 = ipr2 - iprad
                ipr3 = ipr3 - iprad
                tprint = tprint - iprad * dtprint
                itime = itime - 1
                itest = 1
                for i in range(nb):
                    pr[i] = prold[i]
                    sat[i] = satold[i]
                nprtop = nprtop - nprtad
                nprbot = nprbot - nprbad
                nbct = nbctop[nprtop - 1]
                bct = bctop[nprtop - 1]
                nbcb = nbcbot[nprbot - 1]
                bcb = bcbot[nprbot - 1]
                break

            # ---- TRANSPIR (per variety) ----
            for v in range(nv):
                tsum_v[v] = _uptake_variety(v, cover[kc, v], pr, sat, dz, imatl, expn, sres, ksu,
                                            rw, etm, nroot_v, csoil_v, croot_v, rtavg_v, hwilt_v,
                                            tpot_v, trans_v)
            for j in range(nb):
                s_ = 0.0
                for v in range(nv):
                    s_ = s_ + trans_v[v, j]
                trans[j] = s_

            # free drainage: outflow = K(S_ne), lagged at the current iterate
            if nbcb == 3:
                bcb = rksat[imatl[nb - 1]] * krel(sat[nb - 1], imatl[nb - 1], expn, sres)

            # ---- MATRIX ----
            for i in range(1, nb - 1):
                por = porosity(pr[i], imatl[i], poros, stor, prinit[i])
                porold = porosity(prold[i], imatl[i], poros, stor, prinit[i])
                xksat1, xkrel1, dz1 = _interface_k(i - 1, i, dz, imatl, rksat, sat, pr, expn, sres)
                diff1 = xksat1 * xkrel1 / dz1
                xksat2, xkrel2, dz2 = _interface_k(i, i + 1, dz, imatl, rksat, sat, pr, expn, sres)
                diff2 = xksat2 * xkrel2 / dz2
                aa[i, 0] = -diff1
                aa[i, 1] = dz[i] / dt * (por * dsat_dhead(pr[i], imatl[i], alpha, sres, expn, se, he)
                                         + stor[imatl[i]] * sat[i]) + diff1 + diff2
                aa[i, 2] = -diff2
                rhs[i] = -dz[i] / dt * (por * sat[i] - porold * satold[i]) \
                    - diff1 * (pr[i] - pr[i - 1] - dz1) \
                    + diff2 * (pr[i + 1] - pr[i] - dz2) \
                    - trans[i] \
                    - evap_rate(sat[i], imatl[i], swilt, sres, zevap, evmax) * evmult[i] * eweight[i] * dz[i] * esc
            # top block
            por = porosity(pr[0], imatl[0], poros, stor, prinit[0])
            porold = porosity(prold[0], imatl[0], poros, stor, prinit[0])
            xksat2, xkrel2, dz2 = _interface_k(0, 1, dz, imatl, rksat, sat, pr, expn, sres)
            diff2 = xksat2 * xkrel2 / dz2
            ev0 = evap_rate(sat[0], imatl[0], swilt, sres, zevap, evmax) * evmult[0] * eweight[0] * dz[0] * esc
            if nbct == 1:
                pr0 = 2.0 * bct - pr[0]
                sat0 = sat_from_head(bct, imatl[0], expn, sres, alpha, he, se)
                xksat1 = rksat[imatl[0]]
                dz1 = dz[0]
                delh1 = (pr[0] - pr0) / dz1 - 1.0
                if delh1 <= -1.0:
                    coefkr = 1.0
                elif delh1 >= 1.0:
                    coefkr = 0.0
                else:
                    coefkr = 0.5 - delh1 / 2.0
                xkrel1 = krel(sat0, imatl[0], expn, sres) * coefkr + \
                    krel(sat[0], imatl[0], expn, sres) * (1.0 - coefkr)
                diff1 = xksat1 * xkrel1 / dz1
                aa[0, 0] = 0.0
                aa[0, 1] = dz[0] / dt * (por * dsat_dhead(pr[0], imatl[0], alpha, sres, expn, se, he)
                                         + stor[imatl[0]] * sat[0]) + 2.0 * diff1 + diff2
                aa[0, 2] = -diff2
                rhs[0] = -dz[0] / dt * (por * sat[0] - porold * satold[0]) \
                    - diff1 * (pr[0] - pr0 - dz1) \
                    + diff2 * (pr[1] - pr[0] - dz2) \
                    - trans[0] - ev0
            else:
                aa[0, 0] = 0.0
                aa[0, 1] = dz[0] / dt * (por * dsat_dhead(pr[0], imatl[0], alpha, sres, expn, se, he)
                                         + stor[imatl[0]] * sat[0]) + diff2
                aa[0, 2] = -diff2
                rhs[0] = -dz[0] / dt * (por * sat[0] - porold * satold[0]) \
                    + bct \
                    + diff2 * (pr[1] - pr[0] - dz2) \
                    - trans[0] - ev0
            # bottom block
            ne = nb - 1
            por = porosity(pr[ne], imatl[ne], poros, stor, prinit[ne])
            porold = porosity(prold[ne], imatl[ne], poros, stor, prinit[ne])
            xksat1, xkrel1, dz1 = _interface_k(ne - 1, ne, dz, imatl, rksat, sat, pr, expn, sres)
            diff1 = xksat1 * xkrel1 / dz1
            evn = evap_rate(sat[ne], imatl[ne], swilt, sres, zevap, evmax) * evmult[ne] * eweight[ne] * dz[ne] * esc
            if nbcb == 1:
                pr0 = 2.0 * bcb - pr[ne]
                sat0 = sat_from_head(bcb, imatl[ne], expn, sres, alpha, he, se)
                xksat2 = rksat[imatl[ne]]
                dz2 = dz[ne]
                delh2 = (pr0 - pr[ne]) / dz2 - 1.0
                if delh2 <= -1.0:
                    coefkr = 1.0
                elif delh2 >= 1.0:
                    coefkr = 0.0
                else:
                    coefkr = 0.5 - delh2 / 2.0
                xkrel2 = krel(sat[ne], imatl[ne], expn, sres) * coefkr + \
                    krel(sat0, imatl[ne], expn, sres) * (1.0 - coefkr)
                diff2 = xksat2 * xkrel2 / dz2
                aa[ne, 0] = -diff1
                aa[ne, 1] = dz[ne] / dt * (por * dsat_dhead(pr[ne], imatl[ne], alpha, sres, expn, se, he)
                                           + stor[imatl[ne]] * sat[ne]) + diff1 + 2.0 * diff2
                aa[ne, 2] = 0.0
                rhs[ne] = -dz[ne] / dt * (por * sat[ne] - porold * satold[ne]) \
                    - diff1 * (pr[ne] - pr[ne - 1] - dz1) \
                    + diff2 * (pr0 - pr[ne] - dz2) \
                    - trans[ne] - evn
            else:
                aa[ne, 0] = -diff1
                aa[ne, 1] = dz[ne] / dt * (por * dsat_dhead(pr[ne], imatl[ne], alpha, sres, expn, se, he)
                                           + stor[imatl[ne]] * sat[ne]) + diff1
                aa[ne, 2] = 0.0
                rhs[ne] = -dz[ne] / dt * (por * sat[ne] - porold * satold[ne]) \
                    - diff1 * (pr[ne] - pr[ne - 1] - dz1) \
                    - bcb \
                    - trans[ne] - evn

            # ---- THOMAS ----
            aa[0, 0] = 0.0
            aa[0, 2] = aa[0, 2] / aa[0, 1]
            rhs[0] = rhs[0] / aa[0, 1]
            aa[0, 1] = 1.0
            for i in range(1, nb):
                q = aa[i, 1] - aa[i, 0] * aa[i - 1, 2]
                rhs[i] = (rhs[i] - aa[i, 0] * rhs[i - 1]) / q
                aa[i, 0] = 0.0
                aa[i, 1] = 1.0
                aa[i, 2] = aa[i, 2] / q
            princr[nb - 1] = rhs[nb - 1]
            for jj in range(nb - 2, -1, -1):
                princr[jj] = rhs[jj] - aa[jj, 2] * princr[jj + 1]

            # ---- CNVCHK ----
            prmax = 0.0
            for i in range(nb):
                pr[i] = pr[i] + princr[i]
                sat[i] = sat_from_head(pr[i], imatl[i], expn, sres, alpha, he, se)
                if legacy_convergence:
                    ptest = float(np.int64(abs(princr[i])))   # integer truncation, as in CNVCHK.FOR
                else:
                    ptest = abs(princr[i]) - rtol * abs(pr[i])
                if ptest > prmax:
                    prmax = ptest
            ntest = 0
            if prmax > errpr:
                ntest += 1
            if rmax > errres:
                ntest += 1
            if ntest == 0:
                break

        if failed:
            status = STATUS_NOCONV
            break
        if itest == 1:
            continue

        # ---- MASBAL ----
        sumfull = 0.0
        sumroot = 0.0
        for i in range(nb):
            ii = imatl[i]
            snew = sat_from_head(pr[i], ii, expn, sres, alpha, he, se)
            sold = sat_from_head(prold[i], ii, expn, sres, alpha, he, se)
            por = porosity(pr[i], ii, poros, stor, prinit[i])
            porold = porosity(prold[i], ii, poros, stor, prinit[i])
            addw = por * snew - porold * sold
            sumfull = sumfull + addw * dz[i]
            sumroot = sumroot + addw * dz[i] * etm[0, i]
        sumtot = sumtot + sumfull
        sumtotr = sumtotr + sumroot
        qin = 0.0
        if nbct == 2:
            qin = bct
        elif nbct == 1:
            pr0 = 2.0 * bct - pr[0]
            sat0 = sat_from_head(bct, imatl[0], expn, sres, alpha, he, se)
            xksat = rksat[imatl[0]]
            dz1 = dz[0]
            delh1 = (pr[0] - pr0) / dz1 - 1.0
            if delh1 <= -1.0:
                coefkr = 1.0
            elif delh1 >= 1.0:
                coefkr = 0.0
            else:
                coefkr = 0.5 - delh1 / 2.0
            xkrel = krel(sat0, imatl[0], expn, sres) * coefkr + \
                krel(sat[0], imatl[0], expn, sres) * (1.0 - coefkr)
            qin = xksat * xkrel * ((pr0 - pr[0]) / dz[0] + 1.0)
        ne = nb - 1
        qout = 0.0
        if nbcb == 2 or nbcb == 3:
            qout = bcb
        elif nbcb == 1:
            pr0 = 2.0 * bcb - pr[ne]
            sat0 = sat_from_head(bcb, imatl[ne], expn, sres, alpha, he, se)
            xksat = rksat[imatl[ne]]
            dz2 = dz[ne]
            delh2 = (pr0 - pr[ne]) / dz2 - 1.0
            if delh2 <= -1.0:
                coefkr = 1.0
            elif delh2 >= 1.0:
                coefkr = 0.0
            else:
                coefkr = 0.5 - delh2 / 2.0
            xkrel = krel(sat[ne], imatl[ne], expn, sres) * coefkr + \
                krel(sat0, imatl[ne], expn, sres) * (1.0 - coefkr)
            qout = xksat * xkrel * ((pr[ne] - pr0) / dz[ne] + 1.0)
        transum = 0.0
        for i in range(nb):
            transum = transum + dt * trans[i]
        evsum = 0.0
        for i in range(nb):
            evsum = evsum + dt * evap_rate(sat[i], imatl[i], swilt, sres, zevap, evmax) \
                * evmult[i] * eweight[i] * dz[i] * esc
        # flux across the bottom of the root zone (block nleak to nleak+1, 1-based)
        nr = nleak
        xksat, xkrel, dz2 = _interface_k(nr - 1, nr, dz, imatl, rksat, sat, pr, expn, sres)
        qleak = xksat * xkrel * ((pr[nr - 1] - pr[nr]) / dz[nr - 1] + 1.0)
        qin = qin * dt
        qout = qout * dt
        qleak = qleak * dt
        fluxin = fluxin + qin
        flxout = flxout + qout
        fluxet = fluxet + transum
        fluxev = fluxev + evsum
        denom1 = fluxin - flxout - fluxet - fluxev
        if abs(denom1) >= 1.0e-8:
            ratio = sumtot / denom1
        else:
            ratio = 0.0
        denom2 = qin - qout - transum - evsum
        if abs(denom2) >= 1.0e-8:
            ratio1 = sumfull / denom2
        else:
            ratio1 = 0.0

        # per-variety bookkeeping (not in the Fortran)
        for v in range(nv):
            sv = 0.0
            for j in range(nb):
                sv = sv + dt * trans_v[v, j]
                uptake_prof[v, j] += dt * trans_v[v, j]
            dtv[v] += sv
            dtpv[v] += dt * tpot_v[v] * cover[kc, v]
            sr = 0.0
            for j in range(nroot_v[v]):
                sr += wroot[v, j] * sat[j]
            dsr[v] += sr * dt
        dqout += qout

        # ---- update solution vector ----
        for i in range(nb):
            prold[i] = pr[i]
            satold[i] = sat[i]

        # ---- PONDING ----
        if (nbct == 2 and pr[0] > 0.0 and ipond > 0) or iqflag == 1:
            if iqflag == 0:
                if n_pond_begin >= pond_begin.shape[0]:
                    pond_begin = _grow(pond_begin, 2 * pond_begin.shape[0])
                pond_begin[n_pond_begin] = t
                n_pond_begin += 1
                iqflag = 1
                nbct = 1
                bct = pr[0]
            else:
                if exact_rain:
                    qbuild = qbuild + qrain * dt - qin
                else:
                    qbuild = qbuild + bctop[nprtop - 1] * dt - qin
                qinflx = qin / dt
                if n_pond >= pond.shape[0]:
                    pond = _grow(pond, 2 * pond.shape[0])
                pond[n_pond, 0] = t
                pond[n_pond, 1] = qin
                pond[n_pond, 2] = dt
                pond[n_pond, 3] = qbuild
                pond[n_pond, 4] = qinflx
                pond[n_pond, 5] = bctop[nprtop - 1]
                n_pond += 1
                if qbuild > 0.0:
                    nbct = 1
                    if ipond == 1:
                        drunoff += qbuild
                        qbuild = 0.0
                    elif ipond == 2:
                        bct = qbuild
                else:
                    if exact_rain and ipond == 1:
                        # surface water that infiltrated in the last ponded step instead of
                        # running off: take it back out of runoff so rain = infiltration + runoff
                        drunoff += qbuild
                    nbct = nbctop[nprtop - 1]
                    bct = bctop[nprtop - 1]
                    iqflag = 0
                    qbuild = 0.0

        # ---- DAYAVG / PPRINT3 ----
        for k in range(nnd):
            nd = ndepths[k]
            s1 = 0.0
            s2 = 0.0
            for j in range(nd):
                pj = porosity(pr[j], imatl[j], poros, stor, prinit[j])
                s1 = s1 + pj * dz[j]
                s2 = s2 + pj * sat[j] * dz[j]
            daysat[k] = daysat[k] + (s2 / s1) * dt
        dqin = dqin + qin
        dtransum = dtransum + transum
        devsum = devsum + evsum
        dleak = dleak + qleak
        if ipr3 == mprint:
            ipr3 = 0
            if n_drec >= drec.shape[0]:
                drec = _grow(drec, 2 * drec.shape[0])
            storage = 0.0
            for j in range(nb):
                storage += porosity(pr[j], imatl[j], poros, stor, prinit[j]) * sat[j] * dz[j]
            r = drec[n_drec]
            r[0] = t
            r[1] = dqin
            r[2] = dtransum
            r[3] = devsum
            r[4] = dleak
            r[5] = dqout
            r[6] = drunoff
            r[7] = storage
            r[8] = sumtot
            for k in range(nnd):
                r[9 + k] = daysat[k]
            for v in range(nv):
                r[9 + nnd + v] = dtv[v]
                r[9 + nnd + nv + v] = dtpv[v]
                r[9 + nnd + 2 * nv + v] = dsr[v]
            n_drec += 1
            dqin = 0.0
            dtransum = 0.0
            devsum = 0.0
            dleak = 0.0
            dqout = 0.0
            drunoff = 0.0
            for k in range(nnd):
                daysat[k] = 0.0
            for v in range(nv):
                dtv[v] = 0.0
                dtpv[v] = 0.0
                dsr[v] = 0.0

        # ---- PPRINT2 ----
        if ipr2 == mprint:
            ipr2 = 0
            if save_prof:
                if n_prof >= prof.shape[0]:
                    prof = _grow(prof, 2 * prof.shape[0])
                    prof_t = _grow(prof_t, 2 * prof_t.shape[0])
                prof_t[n_prof] = t
                for j in range(nb):
                    prof[n_prof, 0, j] = porosity(pr[j], imatl[j], poros, stor, prinit[j])
                    prof[n_prof, 1, j] = sat[j]
                    prof[n_prof, 2, j] = pr[j]
                n_prof += 1

        # ---- PPRINT ----
        if ipr == 1:
            ipr = 0
            if n_prec >= prec.shape[0]:
                prec = _grow(prec, 2 * prec.shape[0])
            r = prec[n_prec]
            r[0] = t
            r[1] = transum
            r[2] = evsum
            r[3] = qleak
            r[4] = dt
            r[5] = fluxin
            r[6] = flxout
            r[7] = fluxet
            r[8] = fluxev
            r[9] = ratio
            r[10] = ratio1
            for k in range(nnd):
                nd = ndepths[k]
                s1 = 0.0
                s2 = 0.0
                for j in range(nd):
                    pj = porosity(pr[j], imatl[j], poros, stor, prinit[j])
                    s1 = s1 + pj * dz[j]
                    s2 = s2 + pj * sat[j] * dz[j]
                r[11 + k] = s2 / s1
            n_prec += 1
            dt = dtold

        # ---- adjust time step ----
        if it < itmin:
            dt = dt * dtincr
            if dt > dtmax:
                dt = dtmax

    return (status, t, itime, pr, sat,
            prec[:n_prec].copy(), drec[:n_drec].copy(),
            prof_t[:n_prof].copy(), prof[:n_prof].copy(),
            pond[:n_pond].copy(), pond_begin[:n_pond_begin].copy(), uptake_prof)
