from math import inf, pi, tan, isclose

from libdenavit.effective_length_factor import effective_length_factor


def test_k_uninhibited():
    target_K = 3.0
    GA = inf
    GB = 6 / (pi / target_K) / tan(pi / target_K)
    K = effective_length_factor("uninhibited", GA, GB)
    assert isclose(K, target_K)


def test_k_uninhibited_one_zero():
    target_K = 1.5
    GA = 0.0
    GB = -6 * tan(pi / target_K) / (pi / target_K)
    K = effective_length_factor("uninhibited", GA, GB)
    assert isclose(K, target_K)


def test_k_inhibited():
    target_K = 0.8
    GA = 10.0
    GB = 0.45469570131810483
    K = effective_length_factor('inhibited', GA, GB)
    assert isclose(K, target_K)


def test_k_inhibited_one_zero():
    target_K = 0.6
    GA = 0.0
    GB = 0.6067769722307135
    K = effective_length_factor('inhibited', GA, GB)
    assert isclose(K, target_K)


def test_k_basic_cases():
    assert effective_length_factor("uninhibited", inf, inf) == inf
    assert effective_length_factor("uninhibited", inf, 0.0) == 2.0
    assert effective_length_factor("uninhibited", 0.0, inf) == 2.0
    assert effective_length_factor("uninhibited", 0.0, 0.0) == 1.0
    assert effective_length_factor("inhibited", inf, inf) == 1.0
    assert effective_length_factor("inhibited", inf, 0.0) == 0.7
    assert effective_length_factor("inhibited", 0.0, inf) == 0.7
    assert effective_length_factor("inhibited", 0.0, 0.0) == 0.5
