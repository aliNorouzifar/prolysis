"""Declare constraint discovery (pure Python) — a MINERful replacement.

Enumerates candidate Declare constraints over a log's activities, computes the six
MINERful measures (:mod:`prolysis.declare.measures`), applies threshold filtering, and
optionally removes subsumed constraints (``prune="hierarchy"``). Emits the same JSON
structure as MINERful's ``MinerFulMinerStarter -oJSON`` so existing consumers
(``rules_from_json``, the discriminative model builder) work unchanged.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from prolysis.declare.measures import (
    ALL_TEMPLATES, BINARY_TEMPLATES, UNARY_TEMPLATES, measure,
)

Variant = Tuple[Sequence[str], int]

# --- Declare subsumption hierarchy: stronger template -> weaker templates it implies ---
# (same activity parameters). Used by hierarchy pruning to drop redundant constraints.
_SUBSUMES: Dict[str, set] = {
    "ChainResponse": {"Response", "AlternateResponse", "RespondedExistence"},
    "AlternateResponse": {"Response", "RespondedExistence"},
    "Response": {"RespondedExistence"},
    "ChainPrecedence": {"Precedence", "AlternatePrecedence"},
    "AlternatePrecedence": {"Precedence"},
    "ChainSuccession": {"ChainResponse", "ChainPrecedence", "AlternateSuccession",
                        "Succession", "Response", "Precedence", "AlternateResponse",
                        "AlternatePrecedence", "CoExistence", "RespondedExistence"},
    "AlternateSuccession": {"Succession", "AlternateResponse", "AlternatePrecedence",
                            "Response", "Precedence", "CoExistence", "RespondedExistence"},
    "Succession": {"Response", "Precedence", "CoExistence", "RespondedExistence"},
    "CoExistence": {"RespondedExistence"},
    # Init(a)/End(a) guarantee a occurs, so they subsume AtLeast1(a) (MINERful prunes it).
    "Init": {"AtLeast1"},
    "End": {"AtLeast1"},
    "NotCoExistence": {"NotSuccession", "NotResponse", "NotPrecedence",
                       "NotRespondedExistence", "NotChainSuccession", "NotChainResponse",
                       "NotChainPrecedence"},
    "NotSuccession": {"NotResponse", "NotPrecedence", "NotChainSuccession",
                      "NotChainResponse", "NotChainPrecedence"},
    "NotResponse": {"NotChainResponse"},
    "NotPrecedence": {"NotChainPrecedence"},
    "NotChainSuccession": {"NotChainResponse", "NotChainPrecedence"},
    "AtLeast3": {"AtLeast2", "AtLeast1"},
    "AtLeast2": {"AtLeast1"},
    "AtMost1": {"AtMost2", "AtMost3"},
    "AtMost2": {"AtMost3"},
}


def _params_key(c: dict) -> tuple:
    return tuple(p[0] for p in c["parameters"])


def _prune_hierarchy(constraints: List[dict]) -> List[dict]:
    """Drop a constraint when a strictly stronger one over the same activities is present."""
    present = {(c["template"], _params_key(c)) for c in constraints}
    kept = []
    for c in constraints:
        key = _params_key(c)
        subsumed = any((strong, key) in present and c["template"] in _SUBSUMES.get(strong, ())
                       for strong in _SUBSUMES)
        if not subsumed:
            kept.append(c)
    return kept


def discover(variants: List[Variant], activities: Iterable[str], *,
             support: float = 0.0, confidence: float = 0.0,
             prune: str = "hierarchy", templates: Optional[Iterable[str]] = None,
             name: str = "prolysis-declare") -> dict:
    """Discover Declare constraints from ``variants`` (list of ``(trace, frequency)``).

    Returns a MINERful-compatible dict ``{"name", "tasks", "constraints"}``; each
    constraint carries the six measures. With ``prune="none"`` all candidates above the
    thresholds are returned; with ``prune="hierarchy"`` subsumed constraints are removed.
    """
    acts = sorted(activities)
    n_traces = sum(f for _, f in variants)
    n_events = sum(len(t) * f for t, f in variants)
    tmpls = set(templates) if templates is not None else set(ALL_TEMPLATES)

    out: List[dict] = []
    for tmpl in tmpls & UNARY_TEMPLATES:
        for a in acts:
            m = measure(tmpl, [a], variants, n_traces, n_events)
            if m["support"] >= support and m["confidence"] >= confidence:
                out.append({"template": tmpl, "parameters": [[a]], **m})
    for tmpl in tmpls & BINARY_TEMPLATES:
        for a in acts:
            for b in acts:
                if a == b:
                    continue
                m = measure(tmpl, [a, b], variants, n_traces, n_events)
                if m["support"] >= support and m["confidence"] >= confidence:
                    out.append({"template": tmpl, "parameters": [[a], [b]], **m})

    if prune == "hierarchy":
        out = _prune_hierarchy(out)
    elif prune not in ("none", None):
        raise ValueError(f"Unknown prune mode: {prune!r}")

    return {"name": name, "tasks": acts, "constraints": out}


def variants_from_log(log) -> Tuple[List[Variant], set]:
    """Extract ``[(trace, frequency), ...]`` and the activity set from a pm4py log."""
    from pm4py.algo.filtering.log.variants import variants_filter
    V = variants_filter.get_variants(log)
    variants = [(tuple(k), len(v)) for k, v in V.items()]
    activities = {a for t, _ in variants for a in t}
    return variants, activities


def discover_from_log(log_path: str, output_path: Optional[str] = None, *,
                      support: float = 0.0, confidence: float = 0.0,
                      prune: str = "hierarchy",
                      templates: Optional[Iterable[str]] = None) -> dict:
    """Read an XES log, discover constraints, and (optionally) write MINERful-style JSON."""
    import json
    import pm4py

    log = pm4py.read_xes(str(log_path), variant="rustxes")
    variants, activities = variants_from_log(log)
    model = discover(variants, activities, support=support, confidence=confidence,
                     prune=prune, templates=templates)
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as fp:
            json.dump(model, fp, indent=4)
    return model
