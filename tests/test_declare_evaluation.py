"""Unit tests for the pure-Python Declare evaluator (``prolysis.declare``).

These are dependency-free (no scientific stack) — the evaluator is plain Python. The
verdicts below were cross-checked for exact parity against Janus'
``JanusMeasurementsStarter -detailsLevel event`` on the sample logs (9450/9450 labels).
"""
from prolysis.declare import (
    SATISFIED,
    VIOLATED,
    VAC_SATISFIED,
    constraint_key,
    evaluate_trace,
    evaluate_log,
)


def ev(trace, template, *params):
    return evaluate_trace(list(trace), template, list(params))


# --- unary existence / position templates (never vacuous) -------------------------

def test_absence():
    assert ev("xyz", "Absence", "a") == SATISFIED
    assert ev("xay", "Absence", "a") == VIOLATED


def test_participation():
    assert ev("xay", "Participation", "a") == SATISFIED
    assert ev("xyz", "Participation", "a") == VIOLATED


def test_init_end():
    assert ev("abc", "Init", "a") == SATISFIED
    assert ev("bac", "Init", "a") == VIOLATED
    assert ev("abc", "End", "c") == SATISFIED
    assert ev("abc", "End", "a") == VIOLATED


# --- vacuity: activation absent -> vacuously satisfied ----------------------------

def test_vacuity_when_activation_absent():
    assert ev("xyz", "Response", "a", "b") == VAC_SATISFIED       # no a
    assert ev("xyz", "Precedence", "a", "b") == VAC_SATISFIED     # no b
    assert ev("xyz", "RespondedExistence", "a", "b") == VAC_SATISFIED
    assert ev("xyz", "CoExistence", "a", "b") == VAC_SATISFIED    # neither
    assert ev("xyz", "NotCoExistence", "a", "b") == VAC_SATISFIED  # neither


# --- ordering relations -----------------------------------------------------------

def test_response():
    assert ev("axb", "Response", "a", "b") == SATISFIED
    assert ev("axbxa", "Response", "a", "b") == VIOLATED   # trailing a with no later b
    assert ev("axbxax b".replace(" ", ""), "Response", "a", "b") == SATISFIED


def test_precedence():
    assert ev("axb", "Precedence", "a", "b") == SATISFIED
    assert ev("bxa", "Precedence", "a", "b") == VIOLATED   # b before any a


def test_chain_response_precedence():
    assert ev("ab", "ChainResponse", "a", "b") == SATISFIED
    assert ev("axb", "ChainResponse", "a", "b") == VIOLATED
    assert ev("ab", "ChainPrecedence", "a", "b") == SATISFIED
    assert ev("xb", "ChainPrecedence", "a", "b") == VIOLATED


def test_alternate_response_precedence():
    assert ev("abab", "AlternateResponse", "a", "b") == SATISFIED
    assert ev("aab", "AlternateResponse", "a", "b") == VIOLATED   # two a's, no b between
    assert ev("abab", "AlternatePrecedence", "a", "b") == SATISFIED
    assert ev("abb", "AlternatePrecedence", "a", "b") == VIOLATED  # two b's, no a between


def test_coexistence_variants():
    assert ev("xayb", "CoExistence", "a", "b") == SATISFIED   # both present
    assert ev("xay", "CoExistence", "a", "b") == VIOLATED     # only a
    assert ev("xay", "NotCoExistence", "a", "b") == SATISFIED  # only a -> not both
    assert ev("xayb", "NotCoExistence", "a", "b") == VIOLATED  # both -> violated


def test_succession_family():
    assert ev("ab", "Succession", "a", "b") == SATISFIED
    assert ev("ba", "Succession", "a", "b") == VIOLATED       # b before a and a after b
    assert ev("ab", "ChainSuccession", "a", "b") == SATISFIED
    assert ev("axb", "ChainSuccession", "a", "b") == VIOLATED


def test_negative_succession():
    assert ev("axb", "NotSuccession", "a", "b") == VIOLATED   # a eventually followed by b
    assert ev("bxa", "NotSuccession", "a", "b") == SATISFIED  # a never followed by b
    assert ev("ab", "NotChainSuccession", "a", "b") == VIOLATED
    assert ev("axb", "NotChainSuccession", "a", "b") == SATISFIED


# --- keys / log-level API ---------------------------------------------------------

def test_constraint_key_format():
    assert constraint_key("Response", ["a", "b"]) == "Response(a,b)"   # no space
    assert constraint_key("Absence", ["a"]) == "Absence(a)"


def test_evaluate_log_shape_and_key_filtering():
    model = {
        "tasks": ["a", "b"],
        "constraints": [
            {"template": "Response", "parameters": [["a"], ["b"]]},
            {"template": "Absence", "parameters": [["a"]]},
            {"template": "Bogus", "parameters": [["a"], ["b"]]},  # unsupported -> skipped
        ],
    }
    out = evaluate_log(["<a,b>", "<b,a>"], model)
    assert set(out) == {"<a,b>", "<b,a>"}
    assert set(out["<a,b>"]) == {"Response(a,b)", "Absence(a)"}   # Bogus dropped
    assert out["<a,b>"]["Response(a,b)"] == SATISFIED
    assert out["<b,a>"]["Response(a,b)"] == VIOLATED
    assert out["<a,b>"]["Absence(a)"] == VIOLATED
