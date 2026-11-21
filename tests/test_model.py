"""Checks of the field-scale interface and the multi-variety generalization."""
import numpy as np
import pytest

from soilwaterflow.model import Soil, Variety, simulate, simulate_fields

N = 120


@pytest.fixture(scope="module")
def setup():
    rng = np.random.default_rng(3)
    rain = np.where(rng.random(N) < 0.2, rng.exponential(12.0, N), 0.0)
    soil = Soil.sandy_loam_example()
    maize = Variety.from_saturations("maize", 4.5, 0.48, 0.22, soil, root_depth=600, root_shape=(1.0, 3.0))
    sorghum = Variety.from_saturations("sorghum", 3.8, 0.40, 0.17, soil, root_depth=1100, root_shape=(1.2, 1.8))
    return rain, soil, maize, sorghum


def test_water_balance_closes(setup):
    rain, soil, maize, sorghum = setup
    cov = np.full((N, 2), 0.4)
    r = simulate(soil, [maize, sorghum], rain, cover=cov, dz=20.0, dt_max=0.05)
    assert abs(r.water_balance_error()) < 1e-3 * rain.sum()
    # infiltration + runoff = rain
    assert np.allclose(r.infiltration + r.runoff, rain, atol=1e-6)
    # uptake by depth adds up to the daily totals
    assert np.allclose(r.uptake_by_depth.sum(1), r.transpiration.sum(0))


def test_identical_varieties_split_equals_one(setup):
    """Two copies of one variety at cover 0.5 each must behave exactly like one at cover 1."""
    rain, soil, maize, _ = setup
    one = simulate(soil, [maize], rain, cover=np.ones((N, 1)), dz=20.0, dt_max=0.05)
    two = simulate(soil, [maize, maize], rain, cover=np.full((N, 2), 0.5), dz=20.0, dt_max=0.05)
    assert np.allclose(two.transpiration.sum(1), one.transpiration[:, 0], rtol=1e-9, atol=1e-9)
    assert np.allclose(two.transpiration[:, 0], two.transpiration[:, 1], rtol=1e-12, atol=1e-12)
    assert np.allclose(two.storage, one.storage, rtol=1e-9)


def test_actual_never_exceeds_potential(setup):
    rain, soil, maize, sorghum = setup
    r = simulate(soil, [maize, sorghum], rain, cover=np.full((N, 2), 0.5), dz=20.0, dt_max=0.05)
    assert np.all(r.transpiration <= r.transpiration_pot * (1 + 1e-6) + 1e-9)
    assert np.all((r.stress >= -1e-6) & (r.stress <= 1 + 1e-9))


def test_deep_roots_take_deep_water(setup):
    rain, soil, maize, sorghum = setup
    r = simulate(soil, [maize, sorghum], np.zeros(N), cover=np.full((N, 2), 0.5),
                 s0=0.45, dz=20.0, dt_max=0.05)
    below = r.z > maize.root_depth
    assert r.uptake_by_depth[0, below].sum() == 0.0
    assert r.uptake_by_depth[1, below].sum() > 0.0


def test_intense_rain_runs_off(setup):
    _, soil, maize, _ = setup
    rain = np.zeros(10)
    rain[2] = 200.0                     # 200 mm in one hour >> K_s = 800 mm/day
    r = simulate(soil, [maize], rain, storm_hours=1.0, dz=10.0, dt_max=0.02)
    assert r.runoff.sum() > 50.0
    assert abs(r.water_balance_error()) < 0.05


def test_fields_in_parallel(setup):
    rain, soil, maize, sorghum = setup
    rains = np.stack([rain, rain[::-1]])
    rs = simulate_fields(soil, [maize, sorghum], rains, cover=np.full((N, 2), 0.5), dz=25.0, dt_max=0.05)
    single = simulate(soil, [maize, sorghum], rains[1], cover=np.full((N, 2), 0.5), dz=25.0, dt_max=0.05)
    assert np.array_equal(rs[1].transpiration, single.transpiration)
