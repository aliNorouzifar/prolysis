"""Sliding-window Declare behavioural signals (pure Python).

Replaces MINERful's ``MinerFulMinerSlider`` (the ``behavioral_signals.csv`` producer
used by the X-PVI segmentation pipeline). For each window of ``window_size`` consecutive
traces, sliding by ``slide`` traces, it computes every candidate constraint's measures
and returns a per-constraint time series — one value per window.

The output mirrors what ``import_minerful_constraints_timeseries_data`` produced from the
CSV: a list of ``[template, act1, act2, v0, v1, ...]`` rows per measure. As in the CSV,
binary constraints carry ``act2`` with a leading space (``" b"``) and unary constraints
carry ``act2 = ""``; downstream grouping relies on that convention. Measure values are
scaled by 100 (percentages), matching the CSV.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from prolysis.declare.measures import (
    ALL_TEMPLATES, BINARY_TEMPLATES, UNARY_TEMPLATES, measure,
)

Trace = Sequence[str]

# slider measure name -> measure() key
_MEASURE_KEYS = {
    "Support": "support", "Confidence": "confidence", "Coverage": "coverage",
    "Trace support": "tr_support", "Trace confidence": "tr_confidence",
    "Trace coverage": "tr_coverage",
}


def candidate_constraints(activities: Iterable[str],
                          templates: Optional[Iterable[str]] = None
                          ) -> List[Tuple[str, List[str]]]:
    """All candidate constraints (unary over each activity, binary over ordered pairs)."""
    acts = sorted(activities)
    tmpls = set(templates) if templates is not None else set(ALL_TEMPLATES)
    out: List[Tuple[str, List[str]]] = []
    for tmpl in sorted(tmpls & UNARY_TEMPLATES):
        for a in acts:
            out.append((tmpl, [a]))
    for tmpl in sorted(tmpls & BINARY_TEMPLATES):
        for a in acts:
            for b in acts:
                if a != b:
                    out.append((tmpl, [a, b]))
    return out


def _window_starts(n_traces: int, window_size: int, slide: int) -> List[int]:
    """Trace-index window starts: 0, slide, 2*slide, ... while start+window <= n."""
    if window_size >= n_traces:
        return [0]
    return list(range(0, n_traces - window_size + 1, slide))


def sliding_window_signals(traces: Sequence[Trace], window_size: int, slide: int, *,
                           measures: Iterable[str] = ("Confidence", "Coverage"),
                           activities: Optional[Iterable[str]] = None,
                           templates: Optional[Iterable[str]] = None,
                           scale: float = 100.0) -> Dict[str, List[list]]:
    """Compute per-window measure time series for every candidate constraint.

    Parameters
    ----------
    traces:
        Ordered list of traces (each a sequence of activity names). Order matters — the
        windows slide over this sequence, exactly as MINERful slides over the log.
    window_size, slide:
        Window length and step, in number of traces.
    measures:
        Which measures to return (default event ``Confidence`` + ``Coverage``, the two the
        X-PVI pipeline consumes).

    Returns
    -------
    ``{measure_name: [[template, act1, act2, v0, v1, ...], ...]}`` — one row per candidate
    constraint, values scaled by ``scale`` (default 100).
    """
    if activities is None:
        activities = {a for t in traces for a in t}
    candidates = candidate_constraints(activities, templates)
    starts = _window_starts(len(traces), window_size, slide)

    # Pre-compute each window's variants (dedup traces -> frequency) once.
    windows = []
    for s in starts:
        win = traces[s:s + window_size]
        counts: Dict[tuple, int] = {}
        for t in win:
            tt = tuple(t)
            counts[tt] = counts.get(tt, 0) + 1
        variants = list(counts.items())
        n_traces = len(win)
        n_events = sum(len(t) * f for t, f in variants)
        windows.append((variants, n_traces, n_events))

    wanted = list(measures)
    out: Dict[str, List[list]] = {m: [] for m in wanted}
    for tmpl, params in candidates:
        a = params[0]
        act1, act2 = a, (" " + params[1]) if len(params) > 1 else ""
        series = {m: [] for m in wanted}
        for variants, n_traces, n_events in windows:
            m = measure(tmpl, params, variants, n_traces, n_events)
            for name in wanted:
                series[name].append(m[_MEASURE_KEYS[name]] * scale)
        for name in wanted:
            out[name].append([tmpl, act1, act2] + series[name])
    return out


# Full measure set / order emitted by MINERful's slider CSV.
_CSV_MEASURES = ("Support", "Confidence", "Coverage",
                 "Trace support", "Trace confidence", "Trace coverage")


def write_slider_csv(traces: Sequence[Trace], window_size: int, slide: int,
                     out_path: str, *, activities: Optional[Iterable[str]] = None,
                     templates: Optional[Iterable[str]] = None) -> None:
    """Write a ``behavioral_signals.csv`` in MINERful's slider format.

    Drop-in for ``MinerFulMinerSlider -sliOut``: two header rows (constraint names, then
    the six measure names repeated) followed by one row per window
    (``From``; ``To``; values...). Consumers read the ``Confidence`` / ``Coverage``
    columns, which are bit-exact with MINERful.
    """
    import csv

    if activities is None:
        activities = {a for t in traces for a in t}
    candidates = candidate_constraints(activities, templates)
    starts = _window_starts(len(traces), window_size, slide)

    sig = sliding_window_signals(traces, window_size, slide, measures=_CSV_MEASURES,
                                 activities=activities, templates=templates)
    # sig[measure] rows are in `candidates` order; index by (template, act1, act2).
    per_constraint = {}
    for name in _CSV_MEASURES:
        for row in sig[name]:
            per_constraint.setdefault((row[0], row[1], row[2]), {})[name] = row[3:]

    def cname(tmpl, params):
        return f"{tmpl}({params[0]})" if len(params) == 1 else f"{tmpl}({params[0]}, {params[1]})"

    header1 = ["'From'", "'To'"]
    header2 = ["", ""]
    for tmpl, params in candidates:
        header1 += [f"'{cname(tmpl, params)}'"] + [""] * 5
        header2 += [f"'{m}'" for m in _CSV_MEASURES]

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=";")
        w.writerow(header1)
        w.writerow(header2)
        for wi, s in enumerate(starts):
            row = [str(s), str(min(s + window_size, len(traces)))]
            for tmpl, params in candidates:
                a = params[0]
                act2 = (" " + params[1]) if len(params) > 1 else ""
                vals = per_constraint[(tmpl, a, act2)]
                for name in _CSV_MEASURES:
                    row.append(f"{vals[name][wi]:.9f}")
            w.writerow(row)

