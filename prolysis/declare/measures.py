"""MINERful-compatible Declare measures (pure Python).

Reproduces MINERful's per-constraint measures — event-based ``support`` / ``confidence``
/ ``coverage`` and trace-based ``tr_support`` / ``tr_confidence`` / ``tr_coverage`` — so
the discovery step no longer needs the MINERful jar.

Definitions (validated bit-exactly against MINERful ``MinerFulMinerStarter`` on real logs,
all 28 templates, thousands of constraints):

* An **activation** is an event that triggers the constraint's obligation; it is
  **fulfilled** when the obligation is met. Existence/position templates
  (Absence/AtLeast/AtMost/Init/End) have exactly one activation per trace.
* ``coverage``   = (#activation events) / (#events in the log)
* ``support``    = (#fulfilled activations) / (#events in the log)
* ``confidence`` = support / coverage   (= #fulfilled / #activations)
* ``tr_coverage``   = (#traces with >=1 activation) / (#traces)
* ``tr_support``    = (#traces where every activation is fulfilled, and >=1) / (#traces)
* ``tr_confidence`` = tr_support / tr_coverage

Two MINERful quirks are reproduced exactly:

* ``AlternateSuccession`` event fulfilment = Response(a) for the a-part + AlternatePrecedence(b)
  for the b-part (not AlternateResponse+AlternatePrecedence).
* ``Init`` / ``End``: when the ordering constraint is satisfied by *every* trace it is
  reported with ``coverage = support = confidence = 1.0`` (otherwise the 1-activation-per-trace
  formula applies).
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

Trace = Sequence[str]
Variant = Tuple[Sequence[str], int]  # (trace, frequency)

# Templates whose activation is 1-per-trace (existence / position).
UNARY_TEMPLATES = {
    "AtLeast1", "AtLeast2", "AtLeast3", "AtMost1", "AtMost2", "AtMost3",
    "Absence", "Init", "End",
}
BINARY_TEMPLATES = {
    "RespondedExistence", "NotRespondedExistence", "Response", "NotResponse",
    "ChainResponse", "NotChainResponse", "Precedence", "NotPrecedence",
    "ChainPrecedence", "NotChainPrecedence", "AlternateResponse", "AlternatePrecedence",
    "CoExistence", "NotCoExistence", "Succession", "NotSuccession",
    "ChainSuccession", "NotChainSuccession", "AlternateSuccession",
}
ALL_TEMPLATES = UNARY_TEMPLATES | BINARY_TEMPLATES


def _counts(template: str, t: Trace, a: str, b: str | None) -> Tuple[int, int]:
    """(#activations, #fulfilled activations) for one trace."""
    n = len(t)
    if template == "AtLeast1":  return (1, int(t.count(a) >= 1))
    if template == "AtLeast2":  return (1, int(t.count(a) >= 2))
    if template == "AtLeast3":  return (1, int(t.count(a) >= 3))
    if template == "AtMost1":   return (1, int(t.count(a) <= 1))
    if template == "AtMost2":   return (1, int(t.count(a) <= 2))
    if template == "AtMost3":   return (1, int(t.count(a) <= 3))
    if template == "Absence":   return (1, int(t.count(a) == 0))
    if template == "Init":      return (1, int(n > 0 and t[0] == a))
    if template == "End":       return (1, int(n > 0 and t[-1] == a))

    ia = [i for i, e in enumerate(t) if e == a]
    ib = [i for i, e in enumerate(t) if e == b]
    has_a, has_b = bool(ia), bool(ib)

    if template == "RespondedExistence":
        return (len(ia), len(ia) if has_b else 0)
    if template == "NotRespondedExistence":
        return (len(ia), len(ia) if not has_b else 0)
    if template == "Response":
        return (len(ia), sum(1 for i in ia if any(j > i for j in ib)))
    if template == "NotResponse":
        return (len(ia), sum(1 for i in ia if not any(j > i for j in ib)))
    if template == "ChainResponse":
        return (len(ia), sum(1 for i in ia if i + 1 < n and t[i + 1] == b))
    if template == "NotChainResponse":
        return (len(ia), sum(1 for i in ia if not (i + 1 < n and t[i + 1] == b)))
    if template == "Precedence":
        return (len(ib), sum(1 for i in ib if any(j < i for j in ia)))
    if template == "NotPrecedence":
        return (len(ib), sum(1 for i in ib if not any(j < i for j in ia)))
    if template == "ChainPrecedence":
        return (len(ib), sum(1 for i in ib if i - 1 >= 0 and t[i - 1] == a))
    if template == "NotChainPrecedence":
        return (len(ib), sum(1 for i in ib if not (i - 1 >= 0 and t[i - 1] == a)))
    if template == "AlternateResponse":
        ful = 0
        for k, i in enumerate(ia):
            nxt = ia[k + 1] if k + 1 < len(ia) else n
            if any(i < j < nxt for j in ib): ful += 1
        return (len(ia), ful)
    if template == "AlternatePrecedence":
        ful = 0
        for k, i in enumerate(ib):
            prv = ib[k - 1] if k > 0 else -1
            if any(prv < j < i for j in ia): ful += 1
        return (len(ib), ful)
    if template == "CoExistence":
        return (len(ia) + len(ib), (len(ia) if has_b else 0) + (len(ib) if has_a else 0))
    if template == "NotCoExistence":
        return (len(ia) + len(ib), (len(ia) if not has_b else 0) + (len(ib) if not has_a else 0))
    if template == "Succession":
        fa = sum(1 for i in ia if any(j > i for j in ib))
        fb = sum(1 for i in ib if any(j < i for j in ia))
        return (len(ia) + len(ib), fa + fb)
    if template == "NotSuccession":
        fa = sum(1 for i in ia if not any(j > i for j in ib))
        fb = sum(1 for i in ib if not any(j < i for j in ia))
        return (len(ia) + len(ib), fa + fb)
    if template == "ChainSuccession":
        fa = sum(1 for i in ia if i + 1 < n and t[i + 1] == b)
        fb = sum(1 for i in ib if i - 1 >= 0 and t[i - 1] == a)
        return (len(ia) + len(ib), fa + fb)
    if template == "NotChainSuccession":
        fa = sum(1 for i in ia if not (i + 1 < n and t[i + 1] == b))
        fb = sum(1 for i in ib if not (i - 1 >= 0 and t[i - 1] == a))
        return (len(ia) + len(ib), fa + fb)
    if template == "AlternateSuccession":
        # MINERful: Response(a) for the a-part + AlternatePrecedence(b) for the b-part.
        fa = sum(1 for i in ia if any(j > i for j in ib))
        fb = 0
        for k, i in enumerate(ib):
            prv = ib[k - 1] if k > 0 else -1
            if any(prv < j < i for j in ia): fb += 1
        return (len(ia) + len(ib), fa + fb)
    raise ValueError(f"Unsupported template: {template!r}")


def measure(template: str, params: Sequence[str], variants: List[Variant],
            n_traces: int, n_events: int) -> Dict[str, float]:
    """Compute the six MINERful measures for one constraint over the log ``variants``."""
    a = params[0]
    b = params[1] if len(params) > 1 else None
    act = ful = 0            # event-level totals
    tr_cov = tr_sup = 0      # trace-level counts
    for t, f in variants:
        na, nf = _counts(template, t, a, b)
        act += na * f
        ful += nf * f
        if na > 0:
            tr_cov += f
            if nf == na:
                tr_sup += f
    coverage = act / n_events if n_events else 0.0
    support = ful / n_events if n_events else 0.0
    confidence = support / coverage if coverage else 0.0
    # MINERful Init/End quirk: universally-satisfied ordering constraint -> 1.0/1.0/1.0.
    if template in ("Init", "End") and ful == n_traces:
        coverage = support = confidence = 1.0
    tr_coverage = tr_cov / n_traces if n_traces else 0.0
    tr_support = tr_sup / n_traces if n_traces else 0.0
    tr_confidence = tr_support / tr_coverage if tr_coverage else 0.0
    return {
        "support": support, "confidence": confidence, "coverage": coverage,
        "tr_support": tr_support, "tr_confidence": tr_confidence, "tr_coverage": tr_coverage,
    }
