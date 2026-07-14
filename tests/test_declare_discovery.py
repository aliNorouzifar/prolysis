"""Unit tests for the pure-Python MINERful replacement (measures + discovery + pruning).

Measures were validated bit-exactly against the MINERful jar on real logs (all 28
templates, thousands of constraints, and the sliding-window slider). These tests pin the
hand-verifiable cases and the discovery/pruning wiring; they need no scientific stack.
"""
import math

from prolysis.declare.measures import measure
from prolysis.declare.discovery import discover
from prolysis.declare.signals import sliding_window_signals

# Tiny log: "ab" x2, "ba" x1.  N=3, E=6.
VARIANTS = [(("a", "b"), 2), (("b", "a"), 1)]
N, E = 3, 6


def approx(x, y):
    return math.isclose(x, y, abs_tol=1e-9)


def test_response_measures():
    m = measure("Response", ["a", "b"], VARIANTS, N, E)
    assert approx(m["coverage"], 3 / 6)        # 3 a-activations / 6 events
    assert approx(m["support"], 2 / 6)         # only the two "ab" traces fulfil
    assert approx(m["confidence"], 2 / 3)
    assert approx(m["tr_coverage"], 1.0)       # every trace has an a
    assert approx(m["tr_support"], 2 / 3)


def test_existence_measures():
    absc = measure("Absence", ["a"], VARIANTS, N, E)
    assert approx(absc["coverage"], N / E) and approx(absc["support"], 0.0)
    atl = measure("AtLeast1", ["a"], VARIANTS, N, E)
    assert approx(atl["support"], N / E) and approx(atl["confidence"], 1.0)


def test_init_measures_and_universal_quirk():
    init = measure("Init", ["a"], VARIANTS, N, E)   # a is first in 2/3 traces
    assert approx(init["coverage"], N / E) and approx(init["support"], 2 / 6)
    # universally-first activity -> MINERful reports 1.0/1.0/1.0
    uni = [(("x", "a"), 1), (("x", "b"), 1)]
    m = measure("Init", ["x"], uni, 2, 4)
    assert approx(m["coverage"], 1.0) and approx(m["support"], 1.0) and approx(m["confidence"], 1.0)


def test_alternate_succession_uses_response_plus_altprecedence():
    # "a a b": Response(a) fulfils both a's; AlternatePrecedence(b) fulfils the b.
    v = [(("a", "a", "b"), 1)]
    m = measure("AlternateSuccession", ["a", "b"], v, 1, 3)
    assert approx(m["support"], 3 / 3)   # all 3 activations fulfilled


def test_discover_enumeration_and_pruning():
    acts = ["a", "b", "c"]
    full = discover(VARIANTS, acts, support=0.0, confidence=0.0, prune="none")
    # 9 unary templates x 3 activities + 19 binary templates x 6 ordered pairs
    assert len(full["constraints"]) == 9 * 3 + 19 * 6
    assert set(full["tasks"]) == set(acts)
    # hierarchy pruning must not increase the set and must drop subsumed constraints
    pruned = discover(VARIANTS, acts, support=0.0, confidence=0.0, prune="hierarchy")
    assert len(pruned["constraints"]) < len(full["constraints"])


def test_pruning_keeps_strongest_and_drops_subsumed():
    # "ab" only: ChainSuccession(a,b) holds; it subsumes the weaker ordering constraints.
    v = [(("a", "b"), 5)]
    pruned = discover(v, ["a", "b"], support=0.1, confidence=0.0, prune="hierarchy")
    present = {(c["template"], tuple(p[0] for p in c["parameters"])) for c in pruned["constraints"]}
    assert ("ChainSuccession", ("a", "b")) in present       # strongest ordering constraint
    for weaker in ("ChainResponse", "Response", "Succession", "RespondedExistence", "CoExistence"):
        assert (weaker, ("a", "b")) not in present           # all subsumed


def test_sliding_window_signals_shape():
    traces = [("a", "b")] * 10 + [("b", "a")] * 10
    sig = sliding_window_signals(traces, window_size=10, slide=5,
                                 measures=("Confidence", "Coverage"))
    # windows at starts 0,5,10 -> 3 windows
    n_windows = len(sig["Confidence"][0]) - 3   # rows are [template, act1, act2, *values]
    assert n_windows == 3
    assert set(sig) == {"Confidence", "Coverage"}
