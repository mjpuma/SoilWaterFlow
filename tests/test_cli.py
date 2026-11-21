import csv
import os

from soilwaterflow.cli import main

HERE = os.path.dirname(__file__)


def test_quickstart(tmp_path):
    out = tmp_path / "daily.csv"
    assert main(["run", os.path.join(HERE, "..", "examples", "quickstart.json"), "-o", str(out)]) == 0
    rows = list(csv.DictReader(open(out)))
    assert len(rows) == 365
    assert {"T_maize", "Tpot_sorghum", "s_root_sorghum"} <= set(rows[0])
