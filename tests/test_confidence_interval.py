from snc.main import calculate_confidence_interval
import pytest

def test_confidence_interval():
    r = 0.56
    n = 30
    lower, upper = calculate_confidence_interval(n, r)
    assert pytest.approx(lower, abs=1e-3) == 0.2502
    assert pytest.approx(upper, abs=1e-3) == 0.7658
