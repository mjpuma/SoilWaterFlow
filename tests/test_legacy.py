"""The Python port must reproduce the original Fortran (PLANT1D, 2005) exactly.

Reference outputs in tests/data/fortran_reference were produced by compiling the
Fortran in legacy/fortran with gfortran (and they match the outputs committed to the
original SoilWaterFlow repository in 2018). The *_fullprecision files come from the
same code with the .ic and .pl3 formats widened to 17 significant digits."""
import gzip
import os
import shutil

import numpy as np
import pytest

from soilwaterflow.legacy import read_case, run_case, write_outputs

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
REF = os.path.join(HERE, "data", "fortran_reference")
FILES = ["plant1d.pl1", "plant1d.pl2", "plant1d.pl3", "plant1d.pl4",
         "plant1d.anm", "fort.2", "plant1d.out", "plant1d.ic"]


@pytest.fixture(scope="module")
def legacy_run(tmp_path_factory):
    case_dir = tmp_path_factory.mktemp("case")
    src = os.path.join(ROOT, "legacy", "fortran")
    for f in os.listdir(src):
        if f.endswith(".IN"):
            shutil.copy(os.path.join(src, f), case_dir)
    case = read_case(str(case_dir))
    res = run_case(case)
    out = tmp_path_factory.mktemp("out")
    write_outputs(case, res, str(out))
    return case, res, out


def _ref(name):
    with gzip.open(os.path.join(REF, name + ".gz"), "rb") as f:
        return f.read()


@pytest.mark.parametrize("name", FILES)
def test_output_files_identical(legacy_run, name):
    _, _, out = legacy_run
    with open(os.path.join(out, name), "rb") as f:
        assert f.read() == _ref(name), f"{name} differs from the Fortran output"


def test_full_precision(legacy_run):
    _, res, _ = legacy_run
    ic = np.array([float(x) for x in _ref("plant1d_fullprecision.ic").split()])
    assert np.array_equal(ic, res["pr"])
    pl3 = np.array([[float(x) for x in ln.split()] for ln in _ref("plant1d_fullprecision.pl3").decode().splitlines()])
    assert np.array_equal(pl3, res["prec"][:, [0, 5, 6, 7, 8, 9, 10]])


def test_step_count(legacy_run):
    _, res, _ = legacy_run
    assert res["status"] == 0
    assert res["itime"] == 202687
