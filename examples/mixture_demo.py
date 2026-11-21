"""
Does a maize-sorghum mixture use soil water that neither crop uses alone?

Illustrative and uncalibrated: the parameters show how the model behaves, not
findings about real farming systems.

Four stands are grown on the same soil and the same 40 random seasons:
  maize alone, sorghum alone, a 50/50 mixture, and a 50/50 mixture in which
  sorghum is given maize's (shallow) root profile. The last stand removes the
  vertical niche difference, so comparing the two mixtures isolates the effect
  of rooting depth, which a single-bucket soil model cannot represent.

    python examples/mixture_demo.py        # writes examples/mixture_demo.png
"""
import os
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from soilwaterflow.model import Soil, Variety, simulate_fields

N_DAYS = 120          # one rainy season
N_SEASONS = 40
DRY_SPELL = (40, 75)  # days with much less rain (around flowering)


def seasons(n, seed=4):
    """Daily rain for n seasons: Poisson storm arrivals with exponential depths, and a
    mid-season dry spell whose severity varies from season to season."""
    rng = np.random.default_rng(seed)
    rate = np.full(N_DAYS, 0.30)
    out = np.zeros((n, N_DAYS))
    for i in range(n):
        r = rate.copy()
        r[DRY_SPELL[0]:DRY_SPELL[1]] *= rng.uniform(0.0, 0.6)
        wet = rng.random(N_DAYS) < r
        out[i] = np.where(wet, rng.exponential(11.0, N_DAYS), 0.0)
    return out


def canopy(n_days):
    """Seasonal canopy cover: 0.2 at sowing to 1 by day 45, falling to 0.6 at harvest."""
    f = np.arange(n_days) / n_days
    return np.clip(np.minimum(0.2 + 0.8 * f / 0.375, 1.0 - 0.4 * np.clip((f - 0.8) / 0.2, 0, 1)), 0, 1)


def main():
    soil = Soil.sandy_loam_example()
    maize = Variety.from_saturations("maize", t_max=4.5, s_star=0.48, s_w=0.22, soil=soil,
                                     root_depth=700, root_shape=(1.0, 3.0))
    sorghum = Variety.from_saturations("sorghum", t_max=3.8, s_star=0.40, s_w=0.17, soil=soil,
                                       root_depth=1300, root_shape=(1.3, 1.7))
    sorghum_shallow = replace(sorghum, name="sorghum (maize roots)", root_depth=maize.root_depth,
                              root_shape=maize.root_shape)

    rain = seasons(N_SEASONS)
    kc = canopy(N_DAYS)
    kw = dict(depth=2000.0, dz=20.0, dt_max=0.05, s0=0.40, storm_hours=4.0)

    stands = {
        "maize alone": ([maize], np.c_[kc]),
        "sorghum alone": ([sorghum], np.c_[kc]),
        "mixture": ([maize, sorghum], np.c_[0.5 * kc, 0.5 * kc]),
        "mixture, sorghum with maize roots": ([maize, sorghum_shallow], np.c_[0.5 * kc, 0.5 * kc]),
    }
    res = {k: simulate_fields(soil, vs, rain, cover=cov, **kw) for k, (vs, cov) in stands.items()}

    # relative transpiration T / T_pot of each variety over the season (a yield proxy)
    def rel(rs, v):
        return np.array([r.transpiration[:, v].sum() / r.transpiration_pot[:, v].sum() for r in rs])

    rt = {
        "maize": (rel(res["maize alone"], 0), rel(res["mixture"], 0), rel(res["mixture, sorghum with maize roots"], 0)),
        "sorghum": (rel(res["sorghum alone"], 0), rel(res["mixture"], 1), rel(res["mixture, sorghum with maize roots"], 1)),
    }
    et = {k: np.array([r.transpiration.sum() for r in rs]) for k, rs in res.items()}
    drain = {k: np.array([r.drainage.sum() for r in rs]) for k, rs in res.items()}
    wb = max(abs(r.water_balance_error()) for rs in res.values() for r in rs)

    print(f"{N_SEASONS} seasons, mean rain {rain.sum(1).mean():.0f} mm (range {rain.sum(1).min():.0f}-{rain.sum(1).max():.0f})")
    print(f"largest water balance error in any run: {wb:.1e} mm\n")
    print(f"{'stand':38s} {'transpiration':>14s} {'drainage':>9s}   (mm per season, mean)")
    for k in stands:
        print(f"{k:38s} {et[k].mean():14.0f} {drain[k].mean():9.0f}")
    print(f"\nrelative transpiration T/T_pot (mean; driest quarter of seasons)")
    dry = rain.sum(1) <= np.quantile(rain.sum(1), 0.25)
    for crop, (alone, mix, mix_sh) in rt.items():
        print(f"  {crop:8s} alone {alone.mean():.2f} ({alone[dry].mean():.2f})   in mixture {mix.mean():.2f} ({mix[dry].mean():.2f})"
              f"   in mixture, sorghum with maize roots {mix_sh.mean():.2f} ({mix_sh[dry].mean():.2f})")

    # ------------------------------------------------------------------ figure
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.6))
    z = res["mixture"][0].z
    dzs = res["mixture"][0].dz
    up = np.mean([r.uptake_by_depth for r in res["mixture"]], 0) / dzs * 100   # mm per 100 mm of soil
    up_sh = np.mean([r.uptake_by_depth for r in res["mixture, sorghum with maize roots"]], 0) / dzs * 100
    ax[0].plot(up[0], z, color="#1f6f8b", label="maize")
    ax[0].plot(up[1], z, color="#c0582b", label="sorghum")
    ax[0].plot(up_sh[1], z, color="#c0582b", ls=":", label="sorghum given maize roots")
    ax[0].invert_yaxis()
    ax[0].set_xlabel("uptake per season (mm per 100 mm of soil)")
    ax[0].set_ylabel("depth (mm)")
    ax[0].set_title("A. Where each crop draws water in the mixture", loc="left", fontsize=10)
    ax[0].legend(frameon=False, fontsize=8)

    i = int(np.argsort(rain.sum(1))[N_SEASONS // 5])          # a dry season
    days = np.arange(N_DAYS)
    for k, c, ls in [("maize alone", "#1f6f8b", "-"), ("sorghum alone", "#c0582b", "-"),
                     ("mixture", "k", "-"), ("mixture, sorghum with maize roots", "k", ":")]:
        r = res[k][i]
        ax[1].plot(days, np.cumsum(r.transpiration.sum(1)), color=c, ls=ls, label=k)
    ax[1].bar(days, rain[i] * 5, color="0.8", width=1.0, label="rain (x5)")
    ax[1].axvspan(*DRY_SPELL, color="0.93", zorder=0)
    ax[1].set_xlabel("day of season")
    ax[1].set_ylabel("cumulative transpiration (mm)")
    ax[1].set_title(f"B. One dry season ({rain[i].sum():.0f} mm)", loc="left", fontsize=10)
    ax[1].legend(frameon=False, fontsize=7, loc="upper left")

    x = np.arange(2)
    w = 0.26
    for j, (lab, col) in enumerate([("alone", "0.55"), ("in mixture", "k"), ("in mixture, sorghum given maize roots", "0.85")]):
        vals = [rt["maize"][j][dry].mean(), rt["sorghum"][j][dry].mean()]
        ax[2].bar(x + (j - 1) * w, vals, w, color=col, edgecolor="k", lw=0.5, label=lab)
    ax[2].set_xticks(x, ["maize", "sorghum"])
    ax[2].set_ylim(0, 1)
    ax[2].set_ylabel("T / T_pot over the season")
    ax[2].set_title("C. Water stress in the driest quarter of seasons", loc="left", fontsize=10)
    ax[2].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("Illustrative, uncalibrated parameters", fontsize=8, color="0.4", x=0.99, ha="right")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mixture_demo.png")
    fig.savefig(out, dpi=150)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
