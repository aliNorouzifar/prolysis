from enum import Enum, auto
import itertools
import sys
import time
import heapq
import copy

from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
from collections import defaultdict
from itertools import count

import numba as nb
import numpy as np
import numpy.typing as npt
import pm4py
from pm4py.objects.log.obj import EventLog
from pm4py.objects.petri_net.obj import PetriNet, Marking
from pm4py.objects.petri_net.semantics import ClassicSemantics

import ot

from prolysis.analysis.python_emsc.util import time_it


VariantLog = List[Tuple[Tuple[str], int]]

WFNet = Tuple[PetriNet, Marking, Marking]
PNTransition = PetriNet.Transition

class EMSCType(Enum):
    EMSC = auto()
    PEMSC = auto()


@dataclass
class StochasticPetriNet:
    wf_net: WFNet
    weights: Dict[PNTransition, float]


# EventLog implicitly defines a stochastic language
StochasticLanguageObject = Union[EventLog, StochasticPetriNet]


def get_variant_log(log: EventLog) -> VariantLog:
    variants = pm4py.get_variants_as_tuples(log)
    #print(variants)
    #return [(variant, len(traces)) for variant, traces in variants.items()]
    return [(variant, nbr_traces) for variant, nbr_traces in variants.items()]


def log_to_stochastic_language(log: EventLog) -> Dict[Tuple[str], float]:
    variant_log = get_variant_log(log)
    trace_count = sum(freq for _, freq in variant_log)
    return {trace: freq / trace_count for trace, freq in variant_log}


def view_stochastic_petri_net(stochastic_pn: StochasticPetriNet):
    # TODO it is quite painful that PM4Py returns a GraphViz object that cannot be further edited
    #   We work around that by adding the weight to the transition label. Taus become TAU:(weight)
    #   In the future, use my own visualization utilities
    net, im, fm = stochastic_pn.wf_net
    weights = {transition.name: weight for transition, weight in stochastic_pn.weights.items()}

    new_net = copy.deepcopy(net)
    for transition in new_net.transitions:
        if transition.label is None:
            transition.label = f'TAU:{weights[transition.name]}'
        else:
            transition.label = f'{transition.label}:{weights[transition.name]}'

    pm4py.view_petri_net(new_net, im, fm)


def sample_most_likely_net_paths(
        stochastic_net: StochasticPetriNet, prob_mass: Optional[float] = None,
        max_num_paths: Optional[int] = None, max_total_time_ms: Optional[int] = None
) -> List[Tuple[List[PNTransition], float]]:
    """
    Samples the most likely paths from the Petri net until one of the stop criteria is reached.
    If a parameter is None, then it's not considered.
    """
    assert prob_mass is not None or max_num_paths is not None or max_total_time_ms is not None, \
        "At least one of the stopping criteria must be specified."

    if prob_mass is None:
        prob_mass = 1.0
    if max_num_paths is None:
        max_num_paths = sys.maxsize
    if max_total_time_ms is None:
        max_total_time_ms = sys.maxsize
    max_total_time_ns = max_total_time_ms * 1_000_000

    net, im, fm = stochastic_net.wf_net
    weights = stochastic_net.weights

    paths = []
    sampled_mass = 0.0
    semantics = ClassicSemantics()
    marking_transition_cache: Dict[Marking, List[Tuple[PNTransition, Marking, float]]] = {}

    # heapq is one of the few situations where C++ is easier to use than Python
    # We use tuples because the overhead of class
    to_explore = []
    # To break ties in the queue
    id_counter = itertools.count(0)

    # path_id -> (incoming_transition, prev_path)
    path_transition_function = {}

    def _push_heap(_prob, _prev_path_id, _input_transition, _marking):
        path_id = next(id_counter)
        # Keep the heap compact by storing other data outside
        path_transition_function[path_id] = _prev_path_id, _input_transition
        # Must invert the probability to get a max-heap
        heapq.heappush(to_explore, (-_prob, path_id, _marking))

    start_time = time.time_ns()

    def _should_terminate():
        return sampled_mass >= prob_mass or len(paths) >= max_num_paths \
            or (time.time_ns() - start_time) > max_total_time_ns

    _push_heap(1.0, None, None, im)
    it_counter = itertools.count(0)
    while to_explore:
        it_count = next(it_counter)
        if it_count % 100 == 0 and _should_terminate():
            break
        #if it_count % 1000 == 0:
        #    print(sampled_mass, len(paths), len(to_explore))

        cur_prob, cur_path_id, cur_marking = heapq.heappop(to_explore)
        cur_prob = -cur_prob

        if cur_marking not in marking_transition_cache:
            enabled_transitions = semantics.enabled_transitions(net, cur_marking)
            weight_sum = sum(weights[transition] for transition in enabled_transitions)
            marking_transition_cache[cur_marking] = [
                (transition, semantics.weak_execute(transition, net, cur_marking),
                 weights[transition] / weight_sum)
                for transition in enabled_transitions
            ]

        transition_function = marking_transition_cache[cur_marking]
        if len(transition_function) == 0:
            # Reconstruct path
            path = []
            prev_path_id, prev_transition = path_transition_function[cur_path_id]
            while prev_transition is not None:
                path.append(prev_transition)
                prev_path_id, prev_transition = path_transition_function[prev_path_id]
            path.reverse()
            paths.append((path, cur_prob))
            sampled_mass += cur_prob
            continue

        for transition, tgt_marking, transition_prob in transition_function:
            tgt_prob = cur_prob * transition_prob
            assert tgt_prob != 0.0
            _push_heap(tgt_prob, cur_path_id, transition, tgt_marking)

    distinct_paths = {''.join(transition.name for transition in path) for path, _ in paths}
    assert len(distinct_paths) == len(paths), \
        "Cannot sample the same path twice. This indicates an error in the logic"

    path_probs = [path_prob for _, path_prob in paths]
    assert all(prev_prob >= next_prob for prev_prob, next_prob in zip(path_probs, path_probs[1:])), \
        "Paths should be generated from most to least likely"

    return paths


def path_language_to_trace_language(path_language: List[Tuple[List[PNTransition], float]]) \
        -> Dict[Tuple[str], float]:
    # A path language is a sampling over paths. The language maps the path to its labels
    language = defaultdict(lambda: 0.0)
    for path, path_prob in path_language:
        trace = tuple(transition.label for transition in path if transition.label is not None)
        language[trace] += path_prob

    #print(f"PATHS SIZE {len(path_language)} \t LANGUAGE SIZE {len(language)}")
    return dict(language)


@time_it
def compare_languages_levenshtein(lhs: Dict[Tuple[str], float], rhs: Dict[Tuple[str], float], emsc_type: EMSCType=EMSCType.EMSC) \
        -> float:
    # "Assign" traces to IDs that can be used as matrix indexes
    lhs = [*lhs.items()]
    rhs = [*rhs.items()]

    # Compute distance matrix
    #print(f'BUILDING MATRIX')
    if len(lhs) == 0 or len(rhs) == 0:
        raise Exception(f'Cannot compute distance matrix for size {len(lhs)} x {len(rhs)}!')
    dist_matrix = build_distance_matrix(lhs, rhs)

    #print(f'SOLVING SYSTEM')
    prob_src_dist = np.array([src_prob for _, src_prob in lhs])
    # Target distribution is SPN
    prob_tgt_dist = np.array([tgt_prob for _, tgt_prob in rhs])
    # Check if fully unrolled
    diff_mass = np.sum(prob_src_dist) - np.sum(prob_tgt_dist)
    assert diff_mass > -0.0001

    # Add "residual" target, cost function depends on EMSC type
    if diff_mass > 0.000001:
        prob_tgt_dist = np.append(prob_tgt_dist, diff_mass)
        if emsc_type == EMSCType.EMSC:
            # Extend by column that contains minimum distance for each source
            dist_matrix = np.concatenate((dist_matrix, np.min(dist_matrix, axis=1)[:, np.newaxis]), axis=1)
            pass
        elif emsc_type == EMSCType.PEMSC:
            # Extend by column of ones
            dist_matrix = np.concatenate((dist_matrix, np.ones((len(prob_src_dist), 1))), axis=1)
        else:
            raise Exception('Unknown EMSC Type')

    # Build and solve the EMD problem (old lemon wrapper)
    #start_time_wasserstein = time.time_ns()
    #emd_cost = solve_emd(prob_src_dist, prob_tgt_dist, dist_matrix)
    #time_wasserstein = time.time_ns() - start_time_wasserstein

    # OT library (also wrapper to lemon but seems to be more stable)
    #start_time_ot = time.time_ns()
    emd_cost_ot = solve_emd_ot(prob_src_dist, prob_tgt_dist, dist_matrix)

    return emd_cost_ot
    #time_ot = time.time_ns() - start_time_ot

    #if abs(emd_cost - emd_cost_ot) > 0.005:
    #    print(f'PROBLEM EMD MISMATCH: {emd_cost} | {emd_cost_ot} | mass diff {diff_mass}')
    #assert abs(emd_cost - emd_cost_ot) < 0.005

    #if (emd_cost < 0) and (emd_cost_ot < 0):
    #    return (-1, time_wasserstein, time_ot)
    #elif emd_cost_ot >= 0:
    #    return (1.0 - (emd_cost_ot - diff_mass), time_wasserstein, time_ot)
    #else:
    #    return (1.0 - (emd_cost - diff_mass), time_wasserstein, time_ot)



@nb.njit('uint32(uint8[::1], uint8[::1], uint32[::1])', inline='always', cache=True)
def _compiled_distance_lvs_pref_suff_aware(s0: np.array, s1: np.array, previous_row_tmp: np.array):
    start_index = 0
    source_end = len(s0)
    target_end = len(s1)

    while start_index < source_end and start_index < target_end and s0[start_index] == s1[start_index]:
        start_index += 1


    while start_index < source_end and start_index < target_end and s0[source_end - 1] == s1[target_end - 1]:
        source_end -= 1
        target_end -= 1

    if start_index == source_end: 
        return target_end - start_index
    elif start_index == target_end:
        return source_end - start_index


    for i in range(start_index, target_end):
        previous_row_tmp[i] = i + 1 - start_index

    for i in range(start_index, source_end):
        prev_col = i + 1 - start_index
        prev_diag = i - start_index

        for j in range(start_index, target_end):
            prev_row = previous_row_tmp[j]

            insertOrDelete = (prev_col if prev_col < prev_row else prev_row) + 1
            edit = prev_diag + (0 if (s0[i] == s1[j]) else 1)

            prev_diag = prev_row
            prev_col = edit if edit < insertOrDelete else insertOrDelete
            previous_row_tmp[j] = prev_col

    return prev_col


@nb.njit('uint32(uint8[::1], uint8[::1], uint32[::1])', inline='always', cache=True)
def _compiled_distance(s0: np.array, s1: np.array, previous_row_tmp: np.array):
    for i in range(len(s1)):
        previous_row_tmp[i] = i + 1
    
    for i in range(len(s0)):
        prev_col = i + 1
        prev_diag = i

        for j in range(len(s1)):
            prev_row = previous_row_tmp[j]

            insertOrDelete = (prev_col if prev_col < prev_row else prev_row) + 1
            edit = prev_diag + (0 if (s0[i] == s1[j]) else 1)

            prev_diag = prev_row
            prev_col = edit if edit < insertOrDelete else insertOrDelete
            previous_row_tmp[j] = prev_col

    return prev_col


@nb.njit("uint32(uint8[::1], uint8[::1], uint32[::1], uint32[::1])", inline='always', cache=True)
def _compiled_distance_two_rows(s0: np.array, s1: np.array, v0: np.array, v1: np.array):
    # This is a copy of strsimpy's implementation but optimized
    for idx in range(len(s1) + 1):
        v0[idx] = idx
        v1[idx] = 0


    for i in range(len(s0)):
        v1[0] = i + 1
        a = v1[0]
        for j in range(len(s1)):
            cost = 1
            if s0[i] == s1[j]:
                cost = 0
            a = v1[j] + 1
            b = v0[j + 1] + 1
            c = v0[j] + cost
            min_el = a if a <= b else b
            v1[j + 1] = min_el if min_el <= c else c
            # ===========
            # # Even more optimized
            # b = v0[j + 1]
            # min_el = a if a <= b else b
            # min_el += 1
            # if s0_i == s1[j]:
            #     c = v0[j]
            #     min_el = c if c < min_el else min_el
            # v1[j + 1] = min_el
            # a = min_el

        v0, v1 = v1, v0

    return v0[len(s1)]  # Example: Return the final cost (edit distance)


@nb.njit('float64(uint8[::1], uint8[::1], uint32[::1])', inline='always', cache=True)
def normalized_lev_distance(s0: np.array, s1: np.array, np_previous_row_tmp: np.array):
    if len(s0) == 0 and len(s1) == 0:
        return 0.0

    lev_distance = _compiled_distance(s0, s1, np_previous_row_tmp)
    #lev_distance = _compiled_distance_lvs_pref_suff_aware(s0, s1, np_previous_row_tmp)

    return lev_distance / max(1, max(len(s0), len(s1)))


@nb.njit('float64(uint8[::1], uint8[::1], uint32[::1], uint32[::1])', inline='always', cache=True)
def normalized_lev_distance_2(s0: np.array, s1: np.array, v0: np.array, v1):
    if len(s0) == 0 and len(s1) == 0:
        return 0.0

    lev_distance = _compiled_distance_two_rows(s0, s1, v0, v1)

    return lev_distance / max(1, max(len(s0), len(s1)))

def build_distance_matrix(
        src_dist: List[Tuple[Tuple[str], float]],
        tgt_dist: List[Tuple[Tuple[str], float]],
) -> np.ndarray:
    # Map labels to np_arrays
    counter = count(0)
    label_map = defaultdict(lambda: next(counter))

    src_dist = [np.array([label_map[label] for label in trace], dtype=np.uint8) for trace, _ in src_dist]
    tgt_dist = [np.array([label_map[label] for label in trace], dtype=np.uint8) for trace, _ in tgt_dist]

    # To numba
    src_dist_numba = nb.typed.List(src_dist)
    tgt_dist_numba = nb.typed.List(tgt_dist)

    d = fill_distance_matrix(src_dist_numba, tgt_dist_numba, normalized_lev_distance)
    #d = fill_distance_matrix_two_row(src_dist_numba, tgt_dist_numba, normalized_lev_distance_2)
    return d


@time_it
@nb.njit(parallel=True)
def fill_distance_matrix(
        src_language: nb.typed.List[np.array],
        tgt_language: nb.typed.List[np.array],
        distance_function
):
    # We do this ugly implementation because we want Numba to effectively compile it
    dist_matrix = np.empty(shape=(len(src_language), len(tgt_language)), dtype=np.float32)

    for src_id in nb.prange(len(src_language)):
        src_trace = src_language[src_id]
        np_previous_row_tmp = np.zeros(len(src_trace), dtype=np.uint32)
        for tgt_id, tgt_trace in enumerate(tgt_language):
            dist_matrix[src_id, tgt_id] = distance_function(tgt_trace, src_trace, np_previous_row_tmp)

    return dist_matrix


@time_it
@nb.njit(parallel=True)
def fill_distance_matrix_two_row(
        src_language: nb.typed.List[np.array],
        tgt_language: nb.typed.List[np.array],
        distance_function
):
    # We do this ugly implementation because we want Numba to effectively compile it
    dist_matrix = np.empty(shape=(len(src_language), len(tgt_language)), dtype=np.float32)

    for src_id in nb.prange(len(src_language)):
        src_trace = src_language[src_id]
        # Will be filled inside
        v0 = np.zeros(len(src_trace) + 1, dtype=np.uint32)
        v1 = np.zeros(len(src_trace) + 1, dtype=np.uint32)
        for tgt_id, tgt_trace in enumerate(tgt_language):
            dist_matrix[src_id, tgt_id] = distance_function(tgt_trace, src_trace, v0, v1)

    return dist_matrix


@time_it
def solve_emd_ot(
        prob_src_dist: npt.NDArray[np.float32],
        prob_tgt_dist: npt.NDArray[np.float32],
        distance_matrix: npt.NDArray[np.float32]
):
    try:
        # With transport matrix
        #(gamma, log) = ot.emd(prob_src_dist, prob_tgt_dist, distance_matrix, log=True, numItermax=500_000)
        # return log['cost']
        # Without transport matrix
        emd = ot.emd2(prob_src_dist, prob_tgt_dist, distance_matrix, log=False, numItermax=500_000)
        return emd
    except Exception as e:
        print(e)
        return -1
