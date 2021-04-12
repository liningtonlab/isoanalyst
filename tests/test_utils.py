import pytest

from isoanalyst import utils


def test_ppm_tolerance_default():
    mass = 300.0
    expected = (299.997, 300.003)
    actual = utils.ppm_tolerance(mass)
    assert pytest.approx(expected[0]) == actual[0]
    assert pytest.approx(expected[1]) == actual[1]


def test_ppm_tolerance_error_specified():
    mass = 300.0
    ppm_error = 100
    expected = (299.97, 300.03)
    actual = utils.ppm_tolerance(mass, error=100)
    assert pytest.approx(expected[0]) == actual[0]
    assert pytest.approx(expected[1]) == actual[1]


def test_c_isotope_plus():
    mass = 300.0
    assert pytest.approx(301.00335) == utils.c_isotope(mass)


def test_c_isotope_minus():
    mass = 300.0
    assert pytest.approx(298.99665) == utils.c_isotope(mass, sign=-1)


def test_n_isotope_plus():
    mass = 300.0
    assert pytest.approx(300.9970349) == utils.n_isotope(mass)


def test_n_isotope_minus():
    mass = 300.0
    assert pytest.approx(299.0029651) == utils.n_isotope(mass, sign=-1)
