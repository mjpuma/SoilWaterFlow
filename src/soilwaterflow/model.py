"""
Field-scale interface: daily rain and canopy cover in, daily water balance and
transpiration by variety out. Units are mm and days throughout.

The physics is the PLANT1D model (mixed-form Richards equation, block-centered
finite differences, modified Picard iteration; Celia, Guswa and Puma), generalized
so that several crop varieties share one soil column, each with its own root
profile, plant water potential and potential transpiration:

    uptake of variety v from block j =
        c_v w_vj dz_j (h_j - h_p,v) / ( c_s,v / (K_s k_rel(S_j)) + c_r,v )

where c_v is the variety's canopy cover (planting share x seasonal cover), w_vj its
root weights (mean 1 over its root zone), h_j the soil pressure head, and h_p,v the
plant head: h_wilt,v when the soil can't meet demand, otherwise the value at which
total uptake equals c_v T_max,v. The resistances c_s,v and c_r,v are set, as in
INITMVS.FOR, so that uptake equals T_max,v when the whole root zone is at h_star,v,
and a fraction `frac` of the head drop is across the roots when the soil is saturated.

Varieties compete only through the soil water they share: one that roots deeper,
or holds lower plant water potential, can keep transpiring after another has stopped.
"""
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
import math

import numpy as np

from .soil import sat_from_head, head_from_sat, krel
from ._kernel import run_column, STATUS_OK


# ==========================================================================
# Soil and plant parameters
# ==========================================================================

@dataclass
class Soil:
    """Brooks and Corey soil, uniform with depth. Lengths in mm, times in days.

    Saturation S is the fraction of pore space filled with water, the same quantity
    as relative soil moisture s in tessera and in Laio et al. (2001)."""
    porosity: float
    air_entry: float          # air-entry (bubbling) head, mm, positive
    b: float                  # Brooks-Corey exponent (1 / pore-size index)
    k_sat: float              # saturated hydraulic conductivity, mm/day
    s_res: float = 0.0        # residual saturation
    s_e: float = 0.98         # saturation at the air-entry head (curve is linear from s_e to 1 above it)
    storage: float = 1e-11    # specific storage, 1/mm
    e_soil: float = 0.5       # soil evaporation when wetter than h_evap, mm/day
    z_evap: float = 200.0     # depth of the evaporation zone, mm
    evap_shape: tuple = (2.0, 2.0)   # beta(a, b) weights of evaporation over the evaporation zone
    h_evap: float = -316000.0        # head below which evaporation falls linearly to 0 at s_res, mm

    def _arrays(self):
        return (np.array([self.porosity]), np.array([self.air_entry]), np.array([self.b]),
                np.array([self.k_sat]), np.array([self.storage]), np.array([self.s_res]),
                np.array([self.s_e]), np.array([-self.air_entry]))

    def saturation(self, h):
        """Saturation at pressure head h (mm)."""
        _, alpha, expn, _, _, sres, se, he = self._arrays()
        return np.vectorize(lambda x: sat_from_head(float(x), 0, expn, sres, alpha, he, se))(h)

    def head(self, s):
        """Pressure head (mm) at saturation s."""
        _, alpha, expn, _, _, sres, se, he = self._arrays()
        return np.vectorize(lambda x: head_from_sat(float(x), 0, expn, sres, alpha, he, se))(s)

    def conductivity(self, s):
        """Unsaturated hydraulic conductivity K_s k_rel(S), mm/day."""
        _, _, expn, ks, _, sres, _, _ = self._arrays()
        return np.vectorize(lambda x: ks[0] * krel(float(x), 0, expn, sres))(s)

    @classmethod
    def burkea_sandy_clay(cls):
        """The sandy clay of the original PLANT1D example (MATERIAL.IN), converted from cm to mm."""
        return cls(porosity=0.423, air_entry=158.0, b=9.588, k_sat=497.8, s_res=0.40, s_e=0.98,
                   storage=1e-11, e_soil=0.456, z_evap=200.0, h_evap=-316000.0)

    @classmethod
    def sandy_loam_example(cls):
        """Illustrative, uncalibrated sandy loam chosen to sit close to tessera's default
        soil: b = 4.9 gives the Laio et al. leakage exponent beta = 2b + 4 = 13.8, and
        K_s = 800 mm/day. The air-entry head and evaporation values are illustrative."""
        return cls(porosity=0.43, air_entry=150.0, b=4.9, k_sat=800.0, s_res=0.0,
                   e_soil=1.0, z_evap=150.0, h_evap=-150000.0)


@dataclass
class Variety:
    """Water-use traits of one crop variety. Heads in mm (negative), depths in mm."""
    name: str
    t_max: float              # potential transpiration at full cover, mm/day
    h_star: float             # head at which uptake from a uniform root zone just meets t_max
    h_wilt: float             # lowest plant water potential (wilting), mm
    root_depth: float         # mm
    root_shape: tuple = (1.0, 1.0)   # beta(a, b) root weights over depth / root_depth: (1, 1) uniform,
                                     # a < b concentrates roots near the surface
    frac: float = 0.8         # share of the soil-to-leaf head drop that is across the roots in wet soil

    @classmethod
    def from_saturations(cls, name, t_max, s_star, s_w, soil: Soil, root_depth, **kw):
        """Build a variety from the relative soil moisture thresholds used in tessera
        (s* and s_w), converted to heads with the soil's retention curve."""
        return cls(name=name, t_max=t_max, h_star=float(soil.head(s_star)),
                   h_wilt=float(soil.head(s_w)), root_depth=root_depth, **kw)


def _beta_weights(zc, a, b):
    """Unnormalized beta(a, b) density at relative depths zc in (0, 1)."""
    zc = np.clip(zc, 1e-9, 1 - 1e-9)
    return zc ** (a - 1.0) * (1.0 - zc) ** (b - 1.0)


def _zone(dz, zmax, shape):
    """Fraction of each block inside [0, zmax] and beta weights normalized so that
    sum(w * dz * frac) = zmax (mean weight 1), as in GEOMETRY.FOR and betaweight.m."""
    ztop = np.concatenate([[0.0], np.cumsum(dz)[:-1]])
    zbot = ztop + dz
    frac = np.clip((zmax - ztop) / dz, 0.0, 1.0)
    n = int(np.count_nonzero(frac > 0))
    zc = (ztop + 0.5 * np.minimum(dz, np.maximum(zmax - ztop, 0))) / zmax
    w = np.zeros_like(dz)
    w[:n] = _beta_weights(zc[:n], *shape)
    s = np.sum(w * dz * frac)
    w = w * zmax / s if s > 0 else w
    return w, frac, n


# ==========================================================================
# Results
# ==========================================================================

@dataclass
class Result:
    """Daily water balance of one field. Fluxes are mm per day, summed over each day."""
    names: list
    rain: np.ndarray          # (n_days,)
    infiltration: np.ndarray  # into the soil surface
    runoff: np.ndarray        # rain that could not infiltrate
    evaporation: np.ndarray
    transpiration: np.ndarray       # (n_days, n_var) actual
    transpiration_pot: np.ndarray   # (n_days, n_var) unstressed: cover x t_max
    leakage: np.ndarray       # downward flux below the deepest root zone
    drainage: np.ndarray      # out of the bottom of the column
    storage: np.ndarray       # water in the column at the end of each day, mm
    s_root: np.ndarray        # (n_days, n_var) root-weighted mean saturation, averaged over the day
    s_column: np.ndarray      # (n_days,) column-mean saturation, averaged over the day
    storage0: float
    s_top: np.ndarray = None  # (n_days,) mean saturation of the top report_depth mm, if requested
    profiles: dict = field(default_factory=dict)  # 'h', 'sat' (n_days, n_blocks) if requested
    z: np.ndarray = None      # block centers, mm
    dz: np.ndarray = None     # block sizes, mm
    uptake_by_depth: np.ndarray = None   # (n_var, n_blocks) total uptake over the run, mm
    n_steps: int = 0

    @property
    def stress(self):
        """Daily transpiration deficit 1 - T/T_pot for each variety (0 when T_pot = 0)."""
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(self.transpiration_pot > 0,
                            1.0 - self.transpiration / self.transpiration_pot, 0.0)

    def water_balance_error(self):
        """Rain - runoff - ET - drainage - change in storage, mm (should be ~0)."""
        ds = self.storage[-1] - self.storage0
        return float(self.rain.sum() - self.runoff.sum() - self.evaporation.sum()
                     - self.transpiration.sum() - self.drainage.sum() - ds)


# ==========================================================================
# Simulation
# ==========================================================================

def rain_to_flux(rain, storm_hours):
    """Piecewise-constant top boundary flux from daily rain (mm): each day's rain falls
    at a constant rate over the first `storm_hours` hours of the day. Returns the
    PLANT1D boundary arrays (type, rate, end time), padded as the kernel expects."""
    dur = storm_hours / 24.0
    kind, val, tim = [], [], []
    for d, r in enumerate(np.asarray(rain, float)):
        if r > 0 and dur < 1.0:
            kind += [2, 2]
            val += [r / dur, 0.0]
            tim += [d + dur, d + 1.0]
        else:
            if val and val[-1] == r and r == 0.0:
                tim[-1] = d + 1.0           # merge dry days into one interval
                continue
            kind.append(2)
            val.append(r)
            tim.append(d + 1.0)
    kind += [0, 0]
    val += [0.0, 0.0]
    tim += [1e30, 1e30]
    return np.array(kind, np.int64), np.array(val), np.array(tim)


def simulate(soil: Soil, varieties, rain, cover=None, evap_scale="bare",
             depth=1500.0, dz=10.0, bottom="free", s0=0.3,
             storm_hours=4.0, dt_max=0.02, tol=0.1, rtol=1e-4, save_profiles=False,
             report_depth=None):
    """Simulate one field.

    soil        : Soil
    varieties   : list of Variety sharing the field
    rain        : (n_days,) daily rain, mm
    cover       : (n_days, n_var) canopy cover of each variety (planting share x seasonal
                  cover, 0 to 1); default 1 / n_var for every variety on every day
    evap_scale  : multiplier on soil evaporation: 'bare' (1 - total cover, default),
                  a scalar, or an (n_days,) array
    depth, dz   : column depth and block size (mm); dz may be an array of block sizes
    bottom      : 'free' (unit gradient), ('head', h_mm) or ('flux', q_mm_per_day out)
    s0          : initial saturation, scalar or (n_blocks,)
    storm_hours : each day's rain falls at a constant rate over this many hours
    dt_max      : largest time step, days
    tol, rtol   : Picard convergence: |dh| <= tol + rtol |h| in every block (mm)
    report_depth: if given (mm), also report the mean saturation of the top report_depth mm
                  as `s_top` (e.g. to compare with a bucket model's root zone)
    """
    vs = list(varieties)
    nv = len(vs)
    rain = np.asarray(rain, float)
    n_days = rain.size
    if np.ndim(dz) == 0:
        nb = int(round(depth / dz))
        dzs = np.full(nb, float(dz))
    else:
        dzs = np.asarray(dz, float)
        nb = dzs.size
    zc = np.cumsum(dzs) - 0.5 * dzs
    if max(v.root_depth for v in vs) >= dzs.sum() - dzs[-1]:
        raise ValueError("the column must extend at least one block below the deepest roots")

    poros, alpha, expn, rksat, stor, sres, se, he = soil._arrays()
    imatl = np.zeros(nb, np.int64)

    # evaporation zone
    eweight, evmult, _ = _zone(dzs, soil.z_evap, soil.evap_shape)
    swilt = np.array([sat_from_head(soil.h_evap, 0, expn, sres, alpha, he, se)])
    evmax = np.array([soil.e_soil])

    # plants
    rw = np.zeros((nv, nb))
    etm = np.zeros((nv, nb))
    nroot = np.zeros(nv, np.int64)
    csoil = np.zeros(nv)
    croot = np.zeros(nv)
    for k, v in enumerate(vs):
        rw[k], etm[k], nroot[k] = _zone(dzs, v.root_depth, v.root_shape)
        kstar = krel(sat_from_head(v.h_star, 0, expn, sres, alpha, he, se), 0, expn, sres)
        a = v.root_depth / v.t_max
        csoil[k] = a * ((kstar * soil.k_sat) / (1.0 - kstar)) * (v.h_star - v.h_wilt + v.h_wilt * v.frac)
        croot[k] = a * (-v.h_wilt * v.frac) - csoil[k] / soil.k_sat
    nleak = int(nroot.max())

    if cover is None:
        cover = np.full((n_days, nv), 1.0 / nv)
    cover = np.ascontiguousarray(np.broadcast_to(np.asarray(cover, float), (n_days, nv)))
    if isinstance(evap_scale, str):
        if evap_scale != "bare":
            raise ValueError("evap_scale must be 'bare', a number or an array")
        escale = np.clip(1.0 - cover.sum(1), 0.0, 1.0)
    else:
        escale = np.ascontiguousarray(np.broadcast_to(np.asarray(evap_scale, float), (n_days,)))

    nbctop, bctop, timet = rain_to_flux(rain, storm_hours)
    timet[np.flatnonzero(nbctop)[-1]] = float(n_days)
    if bottom == "free":
        nbcb, bcb = 3, 0.0
    elif bottom[0] == "head":
        nbcb, bcb = 1, float(bottom[1])
    elif bottom[0] == "flux":
        nbcb, bcb = 2, float(bottom[1])
    else:
        raise ValueError("bottom must be 'free', ('head', h) or ('flux', q)")
    nbcbot = np.array([nbcb, 0, 0], np.int64)
    bcbot = np.array([bcb, 0.0, 0.0])
    timeb = np.array([float(n_days), 1e30, 1e30])

    ndepths = np.array([nb], np.int64)
    if report_depth is not None:
        ntop = int(np.searchsorted(np.cumsum(dzs), report_depth - 1e-9) + 1)
        ndepths = np.array([nb, min(ntop, nb)], np.int64)

    sat = np.broadcast_to(np.asarray(s0, float), (nb,)).copy()
    pr = np.array([head_from_sat(s, 0, expn, sres, alpha, he, se) for s in sat])
    storage0 = float(np.sum(poros[0] * sat * dzs))

    out = run_column(
        dzs, imatl, evmult, eweight, float(soil.z_evap), 0.0,
        poros, alpha, expn, rksat, stor, sres, se, he, swilt, evmax,
        rw, etm, nroot, csoil, croot, np.ones(nv), np.array([v.h_wilt for v in vs]),
        np.array([v.t_max for v in vs]), np.full(nb, soil.k_sat),
        cover, escale, nleak,
        nbctop, bctop, timet, nbcbot, bcbot, timeb,
        pr, sat, pr.copy(),
        0.0, 1e-4, float(n_days), 10 ** 9, 1.0, 1, 1.0, 1,
        5, 11, 0.5, 2.0, 1e-8, float(dt_max), float(tol), 0.0,
        False, ndepths, float(rtol), bool(save_profiles), True,
    )
    status, t, itime, pr_end, sat_end, prec, drec, prof_t, prof, pond, pond_begin, uptake = out
    if status != STATUS_OK:
        raise RuntimeError(f"simulation stopped at t = {t:.4f} d (status {status})")
    if drec.shape[0] != n_days:
        raise RuntimeError(f"expected {n_days} daily records, got {drec.shape[0]}")
    o = 9 + ndepths.size
    res = Result(
        names=[v.name for v in vs], rain=rain,
        infiltration=drec[:, 1], runoff=drec[:, 6], evaporation=drec[:, 3],
        transpiration=drec[:, o:o + nv], transpiration_pot=drec[:, o + nv:o + 2 * nv],
        leakage=drec[:, 4], drainage=drec[:, 5], storage=drec[:, 7],
        s_root=drec[:, o + 2 * nv:o + 3 * nv], s_column=drec[:, 9],
        storage0=storage0, z=zc, dz=dzs, uptake_by_depth=uptake, n_steps=int(itime),
    )
    if report_depth is not None:
        res.s_top = drec[:, 10]
    if save_profiles:
        res.profiles = dict(t=prof_t, sat=prof[:, 1, :], h=prof[:, 2, :])
    return res


def simulate_fields(soil, varieties, rain, cover=None, n_workers=None, **kw):
    """Simulate many fields in parallel threads (the kernel releases the GIL).

    rain  : (n_fields, n_days)
    cover : None, (n_days, n_var) shared by all fields, or (n_fields, n_days, n_var)
    Other keywords are passed to `simulate`. Returns a list of Result."""
    rain = np.asarray(rain, float)
    nf = rain.shape[0]
    if cover is not None and np.ndim(cover) == 3:
        covs = [cover[i] for i in range(nf)]
    else:
        covs = [cover] * nf
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(simulate, soil, varieties, rain[i], covs[i], **kw) for i in range(nf)]
        return [f.result() for f in futs]
