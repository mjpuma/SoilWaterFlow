"""
Run a PLANT1D case from the original Fortran input files and write the original
output files (plant1d.pl1 ... plant1d.ic), so results can be compared with the
Fortran code byte for byte.

    python -m soilwaterflow.legacy path/to/case_dir [output_dir]

The case directory holds plant1d.in and the files it names (geometry.in,
control.in, material.in, initial.in, bc.in); file names are matched without regard
to case, as on the macOS file system the Fortran was developed on.
"""
import os
import sys
import time
from dataclasses import dataclass

import numpy as np

from .fortran_io import ListReader, fe, ff
from .soil import sat_from_head, head_from_sat, krel
from ._kernel import run_column, STATUS_OK, STATUS_NTMAX, STATUS_NOCONV


def _find(case_dir, name):
    name = name.strip()
    p = os.path.join(case_dir, name)
    if os.path.exists(p):
        return p
    for f in os.listdir(case_dir):
        if f.lower() == name.lower():
            return os.path.join(case_dir, f)
    raise FileNotFoundError(f"{name} not found in {case_dir}")


@dataclass
class LegacyCase:
    """Everything read from a PLANT1D input set, plus the derived quantities that
    INITMVS computes (csoil, croot)."""
    title: str
    # geometry
    nblocks: int
    ztop: float
    zevap: float
    zroot: float
    dz: np.ndarray
    etmult: np.ndarray
    evmult: np.ndarray
    rweight: np.ndarray
    eweight: np.ndarray
    nroot: int
    nevap: int
    # control
    tstart: float
    dt: float
    tmax: float
    ntmax: int
    tday: float
    etstart: float
    etfin: float
    ipond: int
    dtprint: float
    mprint: int
    iscrn: int
    imax: int
    itmin: int
    itmax: int
    dtredu: float
    dtincr: float
    dtmin: float
    dtmax: float
    errpr: float
    errres: float
    # materials (0-based arrays over materials)
    poros: np.ndarray
    alpha: np.ndarray
    expn: np.ndarray
    rksat: np.ndarray
    stor: np.ndarray
    sres: np.ndarray
    se: np.ndarray
    he: np.ndarray
    frac: float
    rtavgdens: np.ndarray
    etmax: np.ndarray
    evmax: np.ndarray
    hwilt: np.ndarray
    hstar: np.ndarray
    swilt: np.ndarray
    imatl: np.ndarray
    # initial condition
    pr: np.ndarray
    sat: np.ndarray
    prinit: np.ndarray
    # boundary conditions
    nbctop: np.ndarray
    bctop: np.ndarray
    timet: np.ndarray
    nbcbot: np.ndarray
    bcbot: np.ndarray
    timeb: np.ndarray
    # derived
    csoil: float = 0.0
    croot: float = 0.0


def read_case(case_dir):
    main = ListReader(_find(case_dir, "plant1d.in"))
    title = main.line()
    fgeo, fcon, fmat, fini, fbc = [main.line()[:12].strip() for _ in range(5)]

    # ---------------- GEOMETRY ----------------
    r = ListReader(_find(case_dir, fgeo))
    nblocks, ztop, zevap, zroot = r.read(4, [int, float, float, float])
    ngrid, nrweight, neweight = r.read(3, int)
    dz = np.zeros(nblocks)
    if ngrid == 1:
        dzconst, nnotdz = r.read(2, int)
        dz[:] = float(dzconst)
        for _ in range(nnotdz):
            a, b, v = r.read(3, [int, int, float])
            dz[a - 1:b] = v
    else:
        dz[:] = r.read(nblocks, float)
    etmult = np.zeros(nblocks)
    evmult = np.zeros(nblocks)
    zbtop = ztop
    nroot = 0
    for i in range(nblocks):
        zbbot = zbtop + dz[i]
        if zroot > zbbot:
            etmult[i] = 1.0
            nroot = i + 1
        elif zroot <= zbtop:
            etmult[i] = 0.0
        else:
            etmult[i] = (zroot - zbtop) / dz[i]
            nroot = i + 1
        zbtop = zbbot
    zbtop = ztop
    nevap = 0
    for i in range(nblocks):
        zbbot = zbtop + dz[i]
        if zevap > zbbot:
            evmult[i] = 1.0
            nevap = i + 1
        elif zevap <= zbtop:
            evmult[i] = 0.0
        else:
            evmult[i] = (zevap - zbtop) / dz[i]
            nevap = i + 1
        zbtop = zbbot
    rweight = np.zeros(nblocks)
    eweight = np.zeros(nblocks)
    if nrweight == 1:
        rweight[:nroot] = 1.0
    else:
        for i in range(nroot):
            rweight[i] = r.read(1, float)[0]
    if neweight == 1:
        eweight[:nevap] = 1.0
    else:
        for i in range(nevap):
            eweight[i] = r.read(1, float)[0]

    # ---------------- CONTROL ----------------
    r = ListReader(_find(case_dir, fcon))
    names = ["tstart", "dt", "tmax", "ntmax", "tday", "etstart", "etfin", "ipond", "dtprint",
             "mprint", "iscrn", "imax", "itmin", "itmax", "dtredu", "dtincr", "dtmin", "dtmax",
             "errpr", "errres"]
    ints = {"ntmax", "ipond", "mprint", "iscrn", "imax", "itmin", "itmax"}
    ctl = {n: r.read(1, int if n in ints else float)[0] for n in names}

    # ---------------- MATERIALS ----------------
    r = ListReader(_find(case_dir, fmat))
    nmatl = r.read(1, int)[0]
    poros, alpha, expn, rksat, stor, sres, se = (np.zeros(nmatl) for _ in range(7))
    rtavgdens, etmax, evmax, hwilt, hstar = (np.zeros(nmatl) for _ in range(5))
    frac = 0.0
    for _ in range(nmatl):
        vals = r.read(8, [int] + [float] * 7)
        m = vals[0] - 1
        poros[m], alpha[m], expn[m], rksat[m], stor[m], sres[m], se[m] = vals[1:]
        frac, rtavgdens[m], etmax[m], evmax[m], hwilt[m], hstar[m] = r.read(6, float)
        if hwilt[m] > 0:
            hwilt[m] = -hwilt[m]
        if hstar[m] > 0:
            hstar[m] = -hstar[m]
    he = -alpha.copy()
    swilt = np.array([sat_from_head(hwilt[m], m, expn, sres, alpha, he, se) for m in range(nmatl)])
    nparam = r.read(1, int)[0]
    imatl = np.zeros(nblocks, dtype=np.int64)
    if nparam == 1:
        matconst, nnotmat = r.read(2, int)
        imatl[:] = matconst - 1
        for _ in range(nnotmat):
            a, b, mt = r.read(3, int)
            imatl[a - 1:b] = mt - 1
    else:
        imatl[:] = np.array(r.read(nblocks, int)) - 1

    # ---------------- INITIAL CONDITION ----------------
    r = ListReader(_find(case_dir, fini))
    iccode, isat = r.read(2, int)
    pr = np.zeros(nblocks)
    sat = np.zeros(nblocks)
    prinit = np.zeros(nblocks)
    hs = lambda s, m: head_from_sat(s, m, expn, sres, alpha, he, se)
    sh = lambda h, m: sat_from_head(h, m, expn, sres, alpha, he, se)
    if iccode == 1:
        ndev, value = r.read(2, [int, float])
        if isat != 1:
            prinit[:] = value
            pr[:] = value
            sat[:] = [sh(value, imatl[i]) for i in range(nblocks)]
            for _ in range(ndev):
                i, v = r.read(2, [int, float])
                prinit[i - 1] = pr[i - 1] = v
                sat[i - 1] = sh(v, imatl[i - 1])
        else:
            sat[:] = value
            prinit[:] = [hs(value, imatl[i]) for i in range(nblocks)]
            pr[:] = prinit
            for _ in range(ndev):
                i, v = r.read(2, [int, float])
                sat[i - 1] = v
                prinit[i - 1] = pr[i - 1] = hs(v, imatl[i - 1])
    elif iccode == 2:
        vals = r.read(nblocks, float)
        if isat == 1:
            sat[:] = vals
            prinit[:] = [hs(sat[i], imatl[i]) for i in range(nblocks)]
            pr[:] = prinit
        else:
            prinit[:] = vals
            pr[:] = prinit
            sat[:] = [sh(pr[i], imatl[i]) for i in range(nblocks)]
    elif iccode == 3:
        zref, prref = r.read(2, float)
        zz = ztop
        for i in range(nblocks):
            zz = zz + dz[i] / 2.0
            prinit[i] = prref + (zz - zref)
            pr[i] = prinit[i]
            sat[i] = sh(pr[i], imatl[i])
            zz = zz + dz[i] / 2.0
    else:
        raise ValueError("improper input for initial conditions")

    # ---------------- BOUNDARY CONDITIONS ----------------
    r = ListReader(_find(case_dir, fbc))
    ntop, nbot = r.read(2, int)

    def read_bc(n):
        kind = np.zeros(n + 2, dtype=np.int64)       # padded: Fortran reads past the end as zeros
        val = np.zeros(n + 2)
        tim = np.full(n + 2, 1e30)      # never reached (the Fortran has zeros here, see README)
        for i in range(n):
            k, v, isat_, tt = r.read(4, [int, float, int, float])
            if isat_ == 1 and k == 1:
                # BCINPUT indexes imatl by the BC number, not a block; for a single
                # material this is the same thing.
                v = hs(v, imatl[min(i, nblocks - 1)])
            kind[i], val[i], tim[i] = k, v, tt
        tim[n - 1] = ctl["tmax"]
        return kind, val, tim

    nbctop, bctop, timet = read_bc(ntop)
    nbcbot, bcbot, timeb = read_bc(nbot)

    case = LegacyCase(title=title, nblocks=nblocks, ztop=ztop, zevap=zevap, zroot=zroot, dz=dz,
                      etmult=etmult, evmult=evmult, rweight=rweight, eweight=eweight,
                      nroot=nroot, nevap=nevap, **ctl,
                      poros=poros, alpha=alpha, expn=expn, rksat=rksat, stor=stor, sres=sres,
                      se=se, he=he, frac=frac, rtavgdens=rtavgdens, etmax=etmax, evmax=evmax,
                      hwilt=hwilt, hstar=hstar, swilt=swilt, imatl=imatl,
                      pr=pr, sat=sat, prinit=prinit,
                      nbctop=nbctop, bctop=bctop, timet=timet,
                      nbcbot=nbcbot, bcbot=bcbot, timeb=timeb)

    # ---------------- INITMVS: root and soil resistances ----------------
    m0 = imatl[0]
    tpotnoon = etmax[m0] - evmax[m0]
    kstar = krel(sat_from_head(hstar[0], m0, expn, sres, alpha, he, se), m0, expn, sres)
    case.csoil = (zroot / tpotnoon) * ((kstar * rksat[0] * rtavgdens[0]) / (1.0 - kstar)) * \
        (hstar[m0] - hwilt[m0] + hwilt[m0] * frac)
    case.croot = (zroot / tpotnoon) * (-hwilt[m0] * frac * rtavgdens[0]) - (case.csoil / rksat[0])
    return case


def run_case(case: LegacyCase, legacy_convergence=True):
    """Run the kernel with a single plant, exactly as PLANT1D does."""
    nb = case.nblocks
    m0 = case.imatl[0]
    rw = case.rweight[None, :].copy()
    etm = case.etmult[None, :].copy()
    out = run_column(
        case.dz, case.imatl, case.evmult, case.eweight, case.zevap, case.ztop,
        case.poros, case.alpha, case.expn, case.rksat, case.stor, case.sres, case.se, case.he,
        case.swilt, case.evmax,
        rw, etm, np.array([case.nroot]), np.array([case.csoil]), np.array([case.croot]),
        np.array([case.rtavgdens[0]]), np.array([case.hwilt[m0]]),
        np.array([case.etmax[m0] - case.evmax[m0]]),
        np.full(nb, case.rksat[0]),          # TRANSPIR uses rksat(1) in every block
        np.ones((1, 1)), np.ones(1), case.nroot,
        case.nbctop, case.bctop, case.timet, case.nbcbot, case.bcbot, case.timeb,
        case.pr, case.sat, case.prinit,
        case.tstart, case.dt, case.tmax, case.ntmax, case.tday, case.ipond, case.dtprint,
        case.mprint, case.itmin, case.itmax, case.dtredu, case.dtincr, case.dtmin, case.dtmax,
        case.errpr, case.errres, legacy_convergence, np.array([40, 100]),
    )
    keys = ["status", "t", "itime", "pr", "sat", "prec", "drec", "prof_t", "prof", "pond", "pond_begin", "uptake"]
    return dict(zip(keys, out))


def write_outputs(case: LegacyCase, res, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    P = lambda n: os.path.join(out_dir, n)
    m0 = case.imatl[0]
    errtpot = float(np.float32(1.0e-8))

    with open(P("plant1d.pl1"), "w") as f1, open(P("plant1d.pl3"), "w") as f3:
        for r in res["prec"]:
            t, transum, evsum, qleak, dt = r[:5]
            fluxin, flxout, fluxet, fluxev, ratio, ratio1 = r[5:11]
            s9, s21 = r[11], r[12]
            tpot = case.etmax[m0] - case.evmax[m0]
            evpot = case.evmax[m0]
            rate_tr, rate_ev, rate_l = transum / dt, evsum / dt, qleak / dt
            if evpot < errtpot:
                evpot = 0.0
            if tpot < errtpot:
                tpot, etevnorm = 0.0, 0.0
            else:
                etevnorm = (rate_tr + rate_ev) / (tpot + evpot)
            f1.write(fe(t, 20, 6) + "".join(ff(x, 10, 4) for x in
                     (etevnorm, tpot, rate_tr, evpot, rate_ev, rate_l, s9, s21)) + "\n")
            f3.write("".join(fe(x, 14, 6) for x in (t, fluxin, flxout, fluxet, fluxev, ratio, ratio1)) + "\n")

    with open(P("plant1d.pl2"), "w") as f:
        for t, pf in zip(res["prof_t"], res["prof"]):
            depth = case.ztop
            for i in range(case.nblocks):
                depth = depth + case.dz[i] / 2.0
                f.write(fe(t, 15, 6) + fe(depth, 15, 4) + ff(pf[0, i], 9, 4) + ff(pf[1, i], 9, 4) + "\n")
                depth = depth + case.dz[i] / 2.0

    with open(P("plant1d.pl4"), "w") as f:
        for r in res["drec"]:
            t, dqin, dtr, dev, dleak = r[:5]
            f.write(fe(t, 20, 6) + "".join(ff(x, 10, 6) for x in
                    (r[9], r[10], dqin, dtr, dev, dleak, case.etmax[0], case.evmax[0])) + "\n")

    with open(P("plant1d.anm"), "w") as fa, open(P("fort.2"), "w") as f2:
        for t, qin, dt, qbuild, qinflx, bct in res["pond"]:
            fa.write("".join(ff(x, 20, 10) for x in (t, qin, dt, qbuild, qinflx, bct)) + "\n")
            f2.write("\n     PONDING INFORMATION (TIME =" + fe(t, 11, 3) + ")\n"
                     "          BUILD-UP OF HEAD AT SURFACE" + fe(qbuild, 19, 3) + "\n"
                     "          CURRENT FLUX INTO DOMAIN AT SURFACE" + fe(qinflx, 11, 3) + "\n"
                     "          ORIGINAL FLUX SPECIFICATION" + fe(bct, 19, 3) + "\n\n\n")

    with open(P("plant1d.out"), "w") as f:
        for t in res["pond_begin"]:
            f.write("\n\n     PONDING HAS BEGUN AT TIME T =" + fe(t, 15, 6) + "\n")
        if res["status"] == STATUS_NOCONV:
            f.write("\n\n     ****************  MAX ITERATIONS EXCEEDED ****************\n")

    with open(P("plant1d.ic"), "w") as f:
        pr = res["pr"]
        for k in range(0, len(pr), 6):
            f.write("  " + "".join(fe(x, 11, 4) for x in pr[k:k + 6]) + "\n")

    for n in ("plant1d.mb", "plant1d.con", "plant1d.pl5", "plant1d.low", "plant1d.md",
              "plant1d.hi", "plant1d.vhi"):
        open(P(n), "w").close()


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 2
    case_dir = argv[0]
    out_dir = argv[1] if len(argv) > 1 else case_dir
    case = read_case(case_dir)
    t0 = time.time()
    res = run_case(case)
    write_outputs(case, res, out_dir)
    print(f"{case.title.strip()}: {res['itime']} time steps, t = {res['t']:.4f}, "
          f"status {res['status']}, {time.time() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
