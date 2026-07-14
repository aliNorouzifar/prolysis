"""Declare discovery / behavioural-signal helpers.

Historically these shelled out to the ``MINERful.jar`` (Java). They now delegate to the
pure-Python engine in :mod:`prolysis.declare`, which reproduces MINERful's measures and
hierarchy pruning bit-exactly (validated against the jar on real logs). Signatures and
output files are unchanged, so existing callers keep working — no Java required.
"""
from __future__ import annotations

import csv
import json

import pandas as pd
import pm4py

from prolysis.declare.discovery import discover_from_log
from prolysis.declare.discovery import _SUBSUMES, _params_key
from prolysis.declare.signals import write_slider_csv


def _ordered_traces(xes_path: str):
    """Read an XES log and return its traces (as tuples) in file order."""
    log = pm4py.read_xes(xes_path, variant="rustxes")
    df = log if isinstance(log, pd.DataFrame) else pm4py.convert_to_dataframe(log)
    return [tuple(g["concept:name"]) for _, g in df.groupby("case:concept:name", sort=False)]


def mine_minerful_for_declare_constraints(window_size, sliding_window_size):
    """Sliding-window behavioural signals -> ``output_files/behavioral_signals.csv``.

    Pure-Python replacement for ``MinerFulMinerSlider``; the ``Confidence`` / ``Coverage``
    columns the X-PVI pipeline consumes are bit-exact with MINERful.
    """
    traces = _ordered_traces("output_files/log_ordered.xes")
    write_slider_csv(traces, int(window_size), int(sliding_window_size),
                     "output_files/behavioral_signals.csv")


def discover_declare(input_log_path, output_log_path, support, confidence):
    """Discover Declare constraints -> MINERful-style JSON at ``output_log_path``.

    Pure-Python replacement for ``MinerFulMinerStarter`` with ``-prune hierarchy``;
    measures and the pruned constraint set are bit-exact with MINERful.
    """
    discover_from_log(input_log_path, output_log_path,
                      support=float(support), confidence=float(confidence),
                      prune="hierarchy")


def prune_constraints_minerful(output_constraint_path, output_constraint_path_pruned):
    """Mark hierarchy-redundant constraints -> CSV (``Constraint``; ``Redudant``).

    Pure-Python replacement for ``MinerFulSimplificationStarter``: a constraint is
    redundant when a strictly stronger one over the same activities is present.
    """
    with open(output_constraint_path, "r") as fh:
        model = json.load(fh)
    constraints = model.get("constraints", [])
    present = {(c["template"], _params_key(c)) for c in constraints}

    with open(output_constraint_path_pruned, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=";")
        w.writerow(["Constraint", "Redudant"])
        for c in constraints:
            key = _params_key(c)
            redundant = any((strong, key) in present
                            and c["template"] in _SUBSUMES.get(strong, ())
                            for strong in _SUBSUMES)
            name = f"{c['template']}({','.join(key)})"
            w.writerow([name, redundant])
