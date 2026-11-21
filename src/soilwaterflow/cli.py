"""
Command line:

    soilwaterflow run CONFIG.json [-o daily.csv]
    soilwaterflow legacy CASE_DIR [OUT_DIR]

`run` simulates one field from a JSON configuration and writes the daily water
balance to CSV. See examples/quickstart.json. `legacy` runs a case from the
original PLANT1D input files and writes the original output files.
"""
import argparse
import csv
import json
import os
import sys
import time

import numpy as np

from .model import Soil, Variety, simulate
from . import legacy as _legacy


def _read_columns(path):
    """Read a CSV with a header row into {column: float array}."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"{path} has no data rows")
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def _soil(spec):
    if isinstance(spec, str):
        if not hasattr(Soil, spec):
            raise ValueError(f"unknown soil preset {spec!r}; use 'sandy_loam_example', "
                             "'burkea_sandy_clay' or a dictionary of Soil parameters")
        return getattr(Soil, spec)()
    if "preset" in spec:
        base = _soil(spec["preset"])
        return Soil(**{**base.__dict__, **{k: v for k, v in spec.items() if k != "preset"}})
    return Soil(**spec)


def _variety(spec, soil):
    spec = dict(spec)
    if "root_shape" in spec:
        spec["root_shape"] = tuple(spec["root_shape"])
    if "s_star" in spec:
        return Variety.from_saturations(soil=soil, **spec)
    return Variety(**spec)


def run(config_path, out_path=None):
    with open(config_path) as f:
        cfg = json.load(f)
    base = os.path.dirname(os.path.abspath(config_path))
    soil = _soil(cfg.get("soil", "sandy_loam_example"))
    vs = [_variety(v, soil) for v in cfg["varieties"]]

    rain_spec = cfg["rain"]
    if isinstance(rain_spec, str):
        cols = _read_columns(os.path.join(base, rain_spec))
        key = "rain_mm" if "rain_mm" in cols else ("rain" if "rain" in cols else list(cols)[-1])
        rain = cols[key]
    else:
        rain = np.asarray(rain_spec, float)
    n = rain.size

    cov_spec = cfg.get("cover")
    if cov_spec is None:
        cover = None
    elif isinstance(cov_spec, str):
        cols = _read_columns(os.path.join(base, cov_spec))
        cover = np.column_stack([cols[v.name] for v in vs])
    else:
        cover = np.broadcast_to(np.asarray(cov_spec, float), (n, len(vs)))

    opts = dict(cfg.get("options", {}))
    if isinstance(opts.get("bottom"), list):
        opts["bottom"] = tuple(opts["bottom"])
    t0 = time.time()
    r = simulate(soil, vs, rain, cover=cover, **opts)

    out_path = out_path or os.path.splitext(config_path)[0] + "_daily.csv"
    header = ["day", "rain", "infiltration", "runoff", "evaporation", "drainage", "leakage",
              "storage", "s_column"]
    for v in vs:
        header += [f"T_{v.name}", f"Tpot_{v.name}", f"s_root_{v.name}"]
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for d in range(n):
            row = [d + 1, r.rain[d], r.infiltration[d], r.runoff[d], r.evaporation[d],
                   r.drainage[d], r.leakage[d], r.storage[d], r.s_column[d]]
            for k in range(len(vs)):
                row += [r.transpiration[d, k], r.transpiration_pot[d, k], r.s_root[d, k]]
            w.writerow([f"{x:.6g}" if isinstance(x, float) else x for x in row])

    print(f"{n} days in {time.time() - t0:.1f} s -> {out_path}")
    print(f"  rain {rain.sum():.1f} mm, runoff {r.runoff.sum():.1f}, evaporation {r.evaporation.sum():.1f}, "
          f"drainage {r.drainage.sum():.1f}, change in storage {r.storage[-1] - r.storage0:.1f}")
    for k, v in enumerate(vs):
        tp = r.transpiration_pot[:, k].sum()
        rel = r.transpiration[:, k].sum() / tp if tp > 0 else float("nan")
        print(f"  {v.name}: transpiration {r.transpiration[:, k].sum():.1f} mm, "
              f"{100 * rel:.0f}% of potential")
    print(f"  water balance error {r.water_balance_error():.1e} mm")
    return r


def main(argv=None):
    p = argparse.ArgumentParser(prog="soilwaterflow", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("run", help="simulate one field from a JSON configuration")
    pr.add_argument("config")
    pr.add_argument("-o", "--output", help="daily CSV (default: CONFIG_daily.csv)")
    pl = sub.add_parser("legacy", help="run original PLANT1D input files")
    pl.add_argument("case_dir")
    pl.add_argument("out_dir", nargs="?")
    a = p.parse_args(argv)
    if a.cmd == "run":
        run(a.config, a.output)
        return 0
    return _legacy.main([a.case_dir] + ([a.out_dir] if a.out_dir else []))


if __name__ == "__main__":
    sys.exit(main())
