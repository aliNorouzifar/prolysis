"""Pure-Python Declare engine.

Replaces the Java tools (MINERful / Janus) for the parts of the pipeline that only
need per-trace constraint evaluation. :func:`evaluate_log` reproduces, in Python, the
per-trace *satisfied / violated / vacuously-satisfied* verdicts that Janus'
``JanusMeasurementsStarter`` produced at ``-detailsLevel event`` (the only thing the
discriminative and LLM tools consumed from it).
"""
from prolysis.declare.evaluation import (
    SATISFIED,
    VIOLATED,
    VAC_SATISFIED,
    constraint_key,
    evaluate_trace,
    evaluate_log,
)

__all__ = [
    "SATISFIED",
    "VIOLATED",
    "VAC_SATISFIED",
    "constraint_key",
    "evaluate_trace",
    "evaluate_log",
]
