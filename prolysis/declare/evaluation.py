"""Per-trace evaluation of Declare constraints (pure Python).

This module reproduces the verdicts the discriminative and LLM tools used to obtain by
calling Janus' ``JanusMeasurementsStarter -detailsLevel event`` and reducing the
per-event RV-LTL codes: a constraint over a trace is

* ``"violated"``      – the constraint is *activated* and its obligation fails;
* ``"satisfied"``     – the constraint is activated and fulfilled;
* ``"v_satisfied"``   – the constraint is never activated (vacuously satisfied).

The three labels are exactly the strings the legacy ``encode()`` step produced from the
Janus output (a ``2`` anywhere in the code array -> violated; else a ``3`` -> satisfied;
else vacuous), so downstream feature encoding is unchanged.

Templates follow the Janus naming used in the model JSON (``Participation`` = AtLeast1,
etc.). Activation = *presence* of the activating activity; when the activation is absent
the verdict is vacuous. When a constraint is activated several times, a single violation
dominates (matches "``2`` present -> violated").
"""
from __future__ import annotations

from typing import Callable, Dict, Iterable, List, Sequence

SATISFIED = "satisfied"
VIOLATED = "violated"
VAC_SATISFIED = "v_satisfied"

Trace = Sequence[str]


# --- unary existence / position templates (never vacuous) -------------------------

def _absence(t: Trace, a: str) -> str:
    # "Absence" here is Janus' Absence1: a must occur zero times.
    return VIOLATED if a in t else SATISFIED


def _participation(t: Trace, a: str) -> str:
    # Participation = AtLeast1.
    return SATISFIED if a in t else VIOLATED


def _init(t: Trace, a: str) -> str:
    return SATISFIED if t and t[0] == a else VIOLATED


def _end(t: Trace, a: str) -> str:
    return SATISFIED if t and t[-1] == a else VIOLATED


# --- binary relation templates ----------------------------------------------------
# Each returns a verdict given the trace and the two activities. Vacuity is decided by
# the activation condition (presence of the activating activity/ies).

def _responded_existence(t: Trace, a: str, b: str) -> str:
    if a not in t:
        return VAC_SATISFIED
    return SATISFIED if b in t else VIOLATED


def _co_existence(t: Trace, a: str, b: str) -> str:
    if a not in t and b not in t:
        return VAC_SATISFIED
    return SATISFIED if (a in t and b in t) else VIOLATED


def _not_co_existence(t: Trace, a: str, b: str) -> str:
    if a not in t and b not in t:
        return VAC_SATISFIED
    return VIOLATED if (a in t and b in t) else SATISFIED


def _response(t: Trace, a: str, b: str) -> str:
    if a not in t:
        return VAC_SATISFIED
    # every a must be eventually followed by a b -> the last a needs a later b.
    last_a = _rindex(t, a)
    return SATISFIED if (b in t[last_a + 1:]) else VIOLATED


def _precedence(t: Trace, a: str, b: str) -> str:
    if b not in t:
        return VAC_SATISFIED
    # every b must be preceded by some a -> the first b needs an earlier a.
    first_b = t.index(b)
    return SATISFIED if (a in t[:first_b]) else VIOLATED


def _alternate_response(t: Trace, a: str, b: str) -> str:
    if a not in t:
        return VAC_SATISFIED
    pending = False  # an a awaiting its b
    for e in t:
        if e == a:
            if pending:
                return VIOLATED  # two a's without a b between
            pending = True
        elif e == b:
            pending = False
    return VIOLATED if pending else SATISFIED


def _alternate_precedence(t: Trace, a: str, b: str) -> str:
    if b not in t:
        return VAC_SATISFIED
    has_a_since_b = False
    for e in t:
        if e == a:
            has_a_since_b = True
        elif e == b:
            if not has_a_since_b:
                return VIOLATED
            has_a_since_b = False
    return SATISFIED


def _chain_response(t: Trace, a: str, b: str) -> str:
    if a not in t:
        return VAC_SATISFIED
    n = len(t)
    for i, e in enumerate(t):
        if e == a and (i + 1 >= n or t[i + 1] != b):
            return VIOLATED
    return SATISFIED


def _chain_precedence(t: Trace, a: str, b: str) -> str:
    if b not in t:
        return VAC_SATISFIED
    for i, e in enumerate(t):
        if e == b and (i == 0 or t[i - 1] != a):
            return VIOLATED
    return SATISFIED


def _succession(t: Trace, a: str, b: str) -> str:
    if a not in t and b not in t:
        return VAC_SATISFIED
    resp_ok = (a not in t) or (b in t[_rindex(t, a) + 1:])
    prec_ok = (b not in t) or (a in t[:t.index(b)])
    return SATISFIED if (resp_ok and prec_ok) else VIOLATED


def _alternate_succession(t: Trace, a: str, b: str) -> str:
    if a not in t and b not in t:
        return VAC_SATISFIED
    resp = _alternate_response(t, a, b)
    prec = _alternate_precedence(t, a, b)
    return SATISFIED if (resp != VIOLATED and prec != VIOLATED) else VIOLATED


def _chain_succession(t: Trace, a: str, b: str) -> str:
    if a not in t and b not in t:
        return VAC_SATISFIED
    resp = _chain_response(t, a, b)
    prec = _chain_precedence(t, a, b)
    return SATISFIED if (resp != VIOLATED and prec != VIOLATED) else VIOLATED


def _not_succession(t: Trace, a: str, b: str) -> str:
    # a is never eventually followed by b.
    if a not in t and b not in t:
        return VAC_SATISFIED
    first_a = t.index(a) if a in t else None
    if first_a is not None and b in t[first_a + 1:]:
        return VIOLATED
    return SATISFIED


def _not_chain_succession(t: Trace, a: str, b: str) -> str:
    # a is never immediately followed by b.
    if a not in t and b not in t:
        return VAC_SATISFIED
    for i, e in enumerate(t[:-1]):
        if e == a and t[i + 1] == b:
            return VIOLATED
    return SATISFIED


def _rindex(t: Trace, x: str) -> int:
    """Index of the last occurrence of ``x`` in ``t``."""
    for i in range(len(t) - 1, -1, -1):
        if t[i] == x:
            return i
    return -1


# Janus template name -> (evaluator, arity)
_TEMPLATES: Dict[str, tuple] = {
    "Absence": (_absence, 1),
    "Participation": (_participation, 1),
    "Init": (_init, 1),
    "End": (_end, 1),
    "RespondedExistence": (_responded_existence, 2),
    "CoExistence": (_co_existence, 2),
    "NotCoExistence": (_not_co_existence, 2),
    "Response": (_response, 2),
    "Precedence": (_precedence, 2),
    "AlternateResponse": (_alternate_response, 2),
    "AlternatePrecedence": (_alternate_precedence, 2),
    "ChainResponse": (_chain_response, 2),
    "ChainPrecedence": (_chain_precedence, 2),
    "Succession": (_succession, 2),
    "AlternateSuccession": (_alternate_succession, 2),
    "ChainSuccession": (_chain_succession, 2),
    "NotSuccession": (_not_succession, 2),
    "NotChainSuccession": (_not_chain_succession, 2),
}


def supported_templates() -> Iterable[str]:
    """Janus template names this engine can evaluate."""
    return _TEMPLATES.keys()


def constraint_key(template: str, params: Sequence[str]) -> str:
    """Reproduce the Janus constraint key, e.g. ``"Response(a,b)"`` (no space)."""
    return f"{template}({','.join(params)})"


def evaluate_trace(trace: Trace, template: str, params: Sequence[str]) -> str:
    """Return the verdict (:data:`SATISFIED` / :data:`VIOLATED` / :data:`VAC_SATISFIED`)
    for a single constraint on a single trace."""
    try:
        fn, arity = _TEMPLATES[template]
    except KeyError as exc:  # pragma: no cover - guarded by the model builder
        raise ValueError(f"Unsupported Declare template: {template!r}") from exc
    return fn(trace, *params[:arity])


def _parse_trace(key: str) -> List[str]:
    """Split a variant key ``"<a,b,c>"`` into its activity list (mirrors the legacy
    ``extract_trace_variants`` encoding, which joins activities with ``,``)."""
    inner = key[1:-1] if key.startswith("<") and key.endswith(">") else key
    return inner.split(",") if inner else []


def evaluate_log(trace_keys: Iterable[str], model: dict) -> Dict[str, Dict[str, str]]:
    """Evaluate every model constraint against every trace variant.

    Parameters
    ----------
    trace_keys:
        Variant keys in the ``"<a,b,c>"`` form produced by ``extract_trace_variants``.
    model:
        A Declare model dict ``{"tasks": [...], "constraints": [{"template", "parameters"}]}``
        using Janus template names (as written by the discriminative tool's model
        builder).

    Returns
    -------
    ``{trace_key: {constraint_key: verdict}}`` — a drop-in replacement for the reduced
    Janus ``[eventsEvaluation]`` output the legacy ``encode()`` step consumed.
    """
    constraints = []
    for c in model.get("constraints", []):
        template = c["template"]
        if template not in _TEMPLATES:
            continue
        params = [p[0] for p in c["parameters"]]
        constraints.append((template, params, constraint_key(template, params)))

    out: Dict[str, Dict[str, str]] = {}
    for key in trace_keys:
        trace = _parse_trace(key)
        verdicts: Dict[str, str] = {}
        for template, params, ckey in constraints:
            verdicts[ckey] = evaluate_trace(trace, template, params)
        out[key] = verdicts
    return out
