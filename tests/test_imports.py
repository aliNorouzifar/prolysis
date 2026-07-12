"""Import smoke tests for prolysis' public entry points.

These guard against structural regressions - e.g. an orphaned subpackage that
imports a module which is not a dependency (as the removed ``util/dfg`` did with
``local_pm4py``). They require the scientific stack (pm4py, numpy, ...), so they are
skipped automatically when it is not installed; run them in a full environment
(Python 3.11 with ``pip install -e .``).
"""
import importlib

import pytest

pytest.importorskip("pm4py", reason="prolysis' scientific stack is not installed")

PUBLIC_MODULES = [
    "prolysis",
    "prolysis.discovery.discovery",
    "prolysis.discovery.subtree_plain",
    "prolysis.discovery.candidate_search.search",
    "prolysis.discovery.candidate_search.is_allowed_2",
    "prolysis.discovery.base_case.check",
    "prolysis.analysis.EMD_based_framework",
    "prolysis.analysis.explainability_extraction",
    "prolysis.analysis.Optimzation_Goals",
    "prolysis.analysis.evaluation",
    "prolysis.calls.minerful_calls",
    "prolysis.rules_handling.declare_processing",
    "prolysis.rules_handling.utils",
]


@pytest.mark.parametrize("module", PUBLIC_MODULES)
def test_public_module_imports(module):
    importlib.import_module(module)


def test_discovery_entry_points_are_exposed():
    from prolysis.discovery.discovery import apply_bi, apply_tree

    assert callable(apply_bi)
    assert callable(apply_tree)


def test_version_is_defined():
    import prolysis

    assert isinstance(prolysis.__version__, str)
