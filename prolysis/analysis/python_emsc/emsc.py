from ast import literal_eval
from collections.abc import Iterable
from dataclasses import dataclass
import dataclasses
import functools
import json
from multiprocessing import Pool
import time
import logging

from collections import defaultdict
from itertools import chain
from pathlib import Path
from typing import Dict, Tuple
import argparse

from pandas import DataFrame
import pm4py
from pm4py.objects.log.obj import EventLog
from pm4py.objects.log.importer.xes import importer as xes_importer
from pm4py.objects.petri_net.importer import importer as pnml_importer
from prolysis.analysis.python_emsc import stochastic as stc
from prolysis.analysis.python_emsc.util import snake_to_camel

logger = logging.getLogger(__name__)

def _process_petri_net(df_ev, net, im):
    # Consider all activities not contained in log silent
    # (Current PM4Py import makes ProM tau transitions visible)
    for t in net.transitions:
        if t.properties['stochastic_distribution'].random_variable.weight == 0:
            logger.warning(f"Transition {str(t)} has weight zero -> set this to a very small value")
            t.properties['stochastic_distribution'].random_variable.weight = 0.0000001


    # No initial marking found
    # Set initial marking assuming that the Petri net
    # is a workflow net
    if len(im) == 0:
        logger.warning("Initial marking of SPN is empty. Creating one under workflow net assumption.")
        for p in net.places:
            # No incoming arcs
            if len(p.in_arcs) == 0:
                im[p] += 1
        logger.warning(f"Created initial marking {str(im)}")


def time_callback(callback, *args, **kwargs):
    start_time = time.perf_counter()
    value = callback(*args, **kwargs)
    end_time = time.perf_counter()
    run_time = end_time - start_time
    return value, run_time


def timer(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return time_callback(func, *args, **kwargs)

    return wrapper


DATASETS_SCRIPT_FOLDER = Path(__file__).parent.resolve()

DATA_FOLDER = DATASETS_SCRIPT_FOLDER / 'data'
MODELS_FOLDER = DATA_FOLDER / 'models'
LOGS_FOLDER = DATA_FOLDER / 'logs'

@dataclass
class MySPNEvalResult:
    log_name: str
    name_spn: str
    nbr_path: int
    sampled_prob: float
    emsc: float
    pemsc: float

    def to_json(self, include_null=False) -> dict:
        """Converts this to json. Assumes variables are snake cased, converts to camel case.

        Args:
            include_null (bool, optional): Whether null values are included. Defaults to False.

        Returns:
            dict: Json dictionary
        """
        return dataclasses.asdict(
            self,
            dict_factory=lambda fields: {
                snake_to_camel(key): value
                for (key, value) in fields
                if value is not None or include_null
            },
        )


def load_log(log_name: str):
    # Tries XES first, else CSV
    log_path_xes = LOGS_FOLDER / log_name / (log_name + ".xes")

    variant = xes_importer.Variants.ITERPARSE
    parameters = {variant.value.Parameters.TIMESTAMP_SORT: True}
    return xes_importer.apply(str(log_path_xes),
                              variant=variant, parameters=parameters)


def read_gspn_file(model: str, variant: str = '') -> stc.StochasticPetriNet:
    model_path = MODELS_FOLDER / model / (model + f'{variant}' + ".pnml")
    wf_net = pm4py.read_pnml(str(model_path), auto_guess_final_marking=True)

    # ProM exports the model wrongly.
    # We must set the taus and estimate the initial marking
    for transition in wf_net[0].transitions:
        if transition.label is not None:
            if 'tau' in transition.label or transition.label in {'start', 'end'}:
                transition.label = None

    if len(wf_net[1]) == 0:
        for place in wf_net[0].places:
            if len(place.in_arcs) == 0:
                wf_net[1][place] = 1

    return get_stochastic_information(wf_net)


def get_stochastic_information(wf_net: stc.WFNet) -> stc.StochasticPetriNet:
    net, im, fm = wf_net

    weights = {}
    for transition in net.transitions:
        stochastic_information = transition.properties['stochastic_distribution']
        weights[transition] = stochastic_information.random_variable.weight

    return stc.StochasticPetriNet(wf_net=wf_net, weights=weights)


def deduplicate_log_labels(log: EventLog):
    # TODO fix this hack
    alphabet = 'abcdefghijklmnopqrstuvxwyz'
    digits = '0123456789'
    tokens = chain(alphabet, alphabet.upper(), digits, [None])
    label_map = defaultdict(lambda: next(tokens))

    for trace in log:
        for event in trace:
            activity = event['concept:name']
            event['concept:name'] = label_map[activity]

    return log, label_map


@timer
def init_model_language(model, max_num_paths=15_000):
    path_language = stc.sample_most_likely_net_paths(model, 0.99, 15_000, 100_000)
    sampled_language = stc.path_language_to_trace_language(path_language)
    return path_language, sampled_language
    

@timer
def compute_emsc(stochastic_log_language, truncated_model_language, emsc_type: stc.EMSCType=stc.EMSCType.EMSC) -> Tuple[float, int, float]:
    try:
        emd = stc.compare_languages_levenshtein(stochastic_log_language, truncated_model_language, emsc_type=emsc_type)
        if emd < 0:
            return -1
        else:
            return 1 - emd
    except Exception as e:
        logger.error(e)
        return -1


def compute_my_metrics(stochastic_log_language, model, max_num_paths=15_000):
    (path_language, truncated_model_language), time_sampling = init_model_language(model, max_num_paths=max_num_paths)
    sampled_prob = sum(prob for _, prob in truncated_model_language.items())
    # Original EMSC
    emsc, time_emsc = compute_emsc(stochastic_log_language, truncated_model_language, emsc_type=stc.EMSCType.EMSC)
    # 1-Penalized EMSC
    pemsc, time_pemsc = compute_emsc(stochastic_log_language, truncated_model_language, emsc_type=stc.EMSCType.PEMSC)
    return emsc, pemsc, len(path_language), sampled_prob


def list_log_path_pn_directory_type(s):
    try:
        log_pn_dir_pairs = literal_eval(s)
    except:
        print(s)
        raise argparse.ArgumentTypeError("Cannot Evaluate (ast.literal_eval): (Log name, Log path, SPN directory), (...), ... should be a enumeration of triples that can be parsed as Python code")
    # Check enumeration of tuples
    if not isinstance(log_pn_dir_pairs, Iterable):
        raise argparse.ArgumentTypeError("Not an enumeration! (Log name, Log path, SPN directory), (...), ... should be a enumeration of triples.")
    if not all(isinstance(item, tuple) for item in log_pn_dir_pairs):
        raise argparse.ArgumentTypeError("Not an enumeration of tuples! (Log name, Log path, SPN directory), (...), ... should be a enumeration of triples.")
    if not all(len(item) == 3 for item in log_pn_dir_pairs):
        raise argparse.ArgumentTypeError("Not an enumeration of triples! (Log name, Log path, SPN directory), (...), ... should be a enumeration of triples.")
    if not all(isinstance(log_name, str) and isinstance(log_path, str) and isinstance(spn_directory, str) for  (log_name, log_path, spn_directory) in log_pn_dir_pairs):
        raise argparse.ArgumentTypeError("Not an enumeration of string triples! (Log name, Log path, SPN directory), (...), ... should be a enumeration of triples.")
    
    log_pn_dir_pairs_as_paths = tuple((log_name, Path(log_path_str), Path(spn_directory)) for (log_name, log_path_str, spn_directory) in log_pn_dir_pairs)
    return log_pn_dir_pairs_as_paths


def _main_emscc_ot_evaluation_full_cmd():
    logging.basicConfig(level=logging.INFO)
    ##############################
    # Arguments
    ##############################
    parser = argparse.ArgumentParser(
        prog='EMSCC on multiple log-SPN combinations',
        description=''
    )

    parser.add_argument('logSPNPairs', help='String pairs (log path, SPN directory) that can be evaluated using ast.literal_eval', 
                        type=list_log_path_pn_directory_type)
    parser.add_argument('--resultFile', type=str, help='Result path must be a json file (directory must exist)')
    parser.add_argument('--poolSize', type=int, help='Pool Size used for concurrent computation', default=10)

    args = parser.parse_args()
    log_SPN_pairs = args.logSPNPairs
    path_result = Path(args.resultFile)
    pool_size = args.poolSize

    logger.info(f"Config read from command line: {str(log_SPN_pairs)}")
    # Check if input files and directories exist
    for (log_name, path_log, path_spn_dir) in log_SPN_pairs:
        if not path_log.is_file():
            raise Exception(f'Event log does not exist: {str(path_log)}')
        if not path_spn_dir.is_dir():
            raise Exception(f'SPN directory does not exist: {str(path_spn_dir)}')

    # Write result header
    with open(path_result, 'a') as f:
        f.write("LogName,SPNName,EMSCC,EMSCCTime,TimeWasserstein,TimeOT,SampledProb,NbrPaths\n")

    emsc_spn_results_full = []
    for log_name, path_log, path_spn_dir in log_SPN_pairs:
        # Load log
        logger.info(f"Importing log {path_log}")
        log = pm4py.read_xes(str(path_log))
        df_ev = pm4py.convert_to_dataframe(log)

        # Stochastic log language
        logger.info(f"Creating stochast language for log {path_log}")
        stochastic_log_language = stc.log_to_stochastic_language(log)


        ##############################
        # Evaluate EMSCC for all SPNs in directory 
        ##############################
        logger.info(f"Start Parallel EMSC Computation with {pool_size} processes")
        # Iterate over Petri nets
        run_spn_evaluation_on_log_partial = functools.partial(run_spn_evaluation_on_log, log_name=log_name, df_ev=df_ev, stochastic_log_language=stochastic_log_language)
        with Pool(processes=pool_size) as pool:
            #emsc_spn_results = pool.map(run_spn_evaluation_on_log_partial, itertools.islice(path_spn_dir.glob('*.pnml'), 3))
            emsc_spn_results = pool.map(run_spn_evaluation_on_log_partial, path_spn_dir.glob('*.pnml'))
        emsc_spn_results_full +=  emsc_spn_results

    # Write results
    data = { 'SPNResultEMSC': list(map(MySPNEvalResult.to_json, emsc_spn_results_full)) }
    with open(path_result, 'w') as f:
        json.dump(data, f, indent=4)


def run_spn_evaluation_on_log(path_pn: Path, log_name: str, df_ev: DataFrame, stochastic_log_language: Dict[Tuple[str], float]) -> MySPNEvalResult:
    net, im, fm, stochastic_map = pnml_importer.apply(str(path_pn), parameters={"return_stochastic_map": True})
    _process_petri_net(df_ev, net, im)
    spn = get_stochastic_information((net, im, fm))

    emsc, pemsc, nbr_path, sampled_prob = compute_my_metrics(stochastic_log_language, spn)
    if emsc < 0 or pemsc < 0:
        logger.error(f'Failed to compute EMSC or PEMS for {path_pn}')

    name_spn = path_pn.name.replace('.pnml', '')
    return MySPNEvalResult(log_name=log_name, name_spn=name_spn,
                           nbr_path=nbr_path, sampled_prob=sampled_prob,
                           emsc=emsc, pemsc=pemsc)


if __name__ == '__main__':
    _main_emscc_ot_evaluation_full_cmd()
    # TODO Distinguish between timed and immediate transitions
    #print(_main_emsc_impl())
