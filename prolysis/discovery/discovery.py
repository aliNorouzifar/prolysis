from typing import Optional, Dict, Any, Tuple
from pm4py.objects.petri_net.obj import PetriNet, Marking
from pm4py.objects.conversion.process_tree import converter as tree_to_petri
from pm4py.objects.process_tree.utils import generic
from pm4py.objects.process_tree.utils.generic import tree_sort
from prolysis.discovery.subtree_plain import SubtreePlain
from pm4py.objects.process_tree.obj import ProcessTree
from pm4py.objects.process_tree.obj import Operator
from pm4py.util import constants
from enum import Enum
from collections import Counter
import json
import pandas as pd
import pm4py
from automata.fa.dfa import DFA
from pm4py.visualization.petri_net import visualizer as pn_visualizer
from pm4py.objects.petri_net.utils.reachability_graph import construct_reachability_graph
import ast
import os
import time
from pathlib import Path
from prolysis.rules_handling.utils import rules_from_json, preprocess, dfa_list_generator
import prolysis.rules_handling.declare_processing as declare_processing
from collections import Counter


class Parameters(Enum):
    ACTIVITY_KEY = constants.PARAMETER_CONSTANT_ACTIVITY_KEY
    START_TIMESTAMP_KEY = constants.PARAMETER_CONSTANT_START_TIMESTAMP_KEY
    TIMESTAMP_KEY = constants.PARAMETER_CONSTANT_TIMESTAMP_KEY
    CASE_ID_KEY = constants.PARAMETER_CONSTANT_CASEID_KEY



def apply_bi(Lp=pd.DataFrame(), Lm=pd.DataFrame(), parameters: Optional[Dict[Any, Any]] = None, sup= None, ratio = None, noise_thr =None, size_par = None, rules =None) -> Tuple[PetriNet, Marking, Marking]:
    file_path = r'output_files\discovery_log.json'
    with open(file_path, 'w') as file:
        json.dump([], file, indent=4)
    process_tree = apply_tree(Lp, Lm, parameters, sup=sup, ratio=ratio, noise_thr=noise_thr, size_par=size_par, rules=rules)
    net, initial_marking, final_marking = tree_to_petri.apply(process_tree)

    return net, initial_marking, final_marking



def apply_tree(logp,logm, parameters=None, sup= None, ratio = None, noise_thr =None, size_par = None, rules= None):
    if parameters is None:
        parameters = {}
    contains_empty_traces = False
    traces_length = [len(trace) for trace in logp]
    if traces_length:
        contains_empty_traces = min([len(trace) for trace in logp]) == 0

    recursion_depth = 0

    # logP_var = Counter(tuple([x['concept:name'] for x in t]) for t in logp)
    # logP_var = pm4py.stats.get_variants(logp)
    #
    # logM_var = pm4py.stats.get_variants(logm)
    if logm.empty:
        logm_var = Counter()
    else:
        logm_var = pm4py.stats.get_variants(logm)
    logp_var = pm4py.stats.get_variants(logp)
    sub = SubtreePlain(logp_var,logm_var, recursion_depth, noise_thr, parameters=parameters, sup= sup, ratio = ratio, size_par = size_par, rules= rules)

    process_tree = get_repr(sub, 0, contains_empty_traces=contains_empty_traces)
    # Ensures consistency to the parent pointers in the process tree
    fix_parent_pointers(process_tree)
    # Fixes a 1 child XOR that is added when single-activities flowers are found
    fix_one_child_xor_flower(process_tree)
    # folds the process tree (to simplify it in case fallthroughs/filtering is applied)
    process_tree = generic.fold(process_tree)
    # sorts the process tree to ensure consistency in different executions of the algorithm
    tree_sort(process_tree)

    return process_tree



def get_repr(spec_tree_struct, rec_depth, contains_empty_traces=False):



    base_cases = ('empty_log', 'single_activity')
    cut = ('xor', 'sequential', 'parallel', 'loopCut')


    # if a cut was detected in the current subtree:
    if spec_tree_struct.detected_cut in cut:
        if spec_tree_struct.detected_cut == "sequential":
            final_tree_repr = ProcessTree(operator=Operator.SEQUENCE)
        elif spec_tree_struct.detected_cut == "loopCut":
            final_tree_repr = ProcessTree(operator=Operator.LOOP)
        elif spec_tree_struct.detected_cut == "xor":
            final_tree_repr = ProcessTree(operator=Operator.XOR)
        elif spec_tree_struct.detected_cut == "parallel":
            final_tree_repr = ProcessTree(operator=Operator.PARALLEL)

        if not (spec_tree_struct.detected_cut == "loopCut" and len(spec_tree_struct.children) >= 3):
            for ch in spec_tree_struct.children:
                # get the representation of the current child (from children in the subtree-structure):
                child = get_repr(ch, rec_depth + 1)
                # add connection from child_tree to child_final and the other way around:
                final_tree_repr.children.append(child)
                child.parent = final_tree_repr

        else:
            child = get_repr(spec_tree_struct.children[0], rec_depth + 1)
            final_tree_repr.children.append(child)
            child.parent = final_tree_repr

            redo_child = ProcessTree(operator=Operator.XOR)
            for ch in spec_tree_struct.children[1:]:
                child = get_repr(ch, rec_depth + 1)
                redo_child.children.append(child)
                child.parent = redo_child

            final_tree_repr.children.append(redo_child)
            redo_child.parent = final_tree_repr

        if spec_tree_struct.detected_cut == "loopCut" and len(spec_tree_struct.children) < 3:
            while len(spec_tree_struct.children) < 2:
                child = ProcessTree()
                final_tree_repr.children.append(child)
                child.parent = final_tree_repr
                spec_tree_struct.children.append(None)

    if spec_tree_struct.detected_cut in base_cases:
        # in the base case of an empty log, we only return a silent transition
        if spec_tree_struct.detected_cut == "empty_log":
            return ProcessTree(operator=None, label=None)
        # in the base case of a single activity, we return a tree consisting of the single activity
        elif spec_tree_struct.detected_cut == "single_activity":
            act_a = spec_tree_struct.activitiesP.pop()
            return ProcessTree(operator=None, label=act_a)
    return final_tree_repr


def get_transition(label):
    """
    Create a node (transition) with the specified label in the process tree
    """
    return ProcessTree(operator=None, label=label)

class Counts(object):
    """
    Shared variables among executions
    """

    def __init__(self):
        """
        Constructor
        """
        self.num_places = 0
        self.num_hidden = 0
        self.num_visible_trans = 0
        self.dict_skips = {}
        self.dict_loops = {}

    def inc_places(self):
        """
        Increase the number of places
        """
        self.num_places = self.num_places + 1

    def inc_no_hidden(self):
        """
        Increase the number of hidden transitions
        """
        self.num_hidden = self.num_hidden + 1

    def inc_no_visible(self):
        """
        Increase the number of visible transitions
        """
        self.num_visible_trans = self.num_visible_trans + 1


def fix_parent_pointers(pt):
    """
    Ensures consistency to the parent pointers in the process tree

    Parameters
    --------------
    pt
        Process tree
    """
    for child in pt.children:
        child.parent = pt
        if child.children:
            fix_parent_pointers(child)


def fix_one_child_xor_flower(tree):
    """
    Fixes a 1 child XOR that is added when single-activities flowers are found

    Parameters
    --------------
    tree
        Process tree
    """
    if tree.parent is not None and tree.operator is Operator.XOR and len(tree.children) == 1:
        for child in tree.children:
            child.parent = tree.parent
            tree.parent.children.append(child)
            del tree.parent.children[tree.parent.children.index(tree)]
    else:
        for child in tree.children:
            fix_one_child_xor_flower(child)


def run_IMr(LPlus_LogFile,support,rules, activities,dim,abs_thr):
    event_log_xes = pm4py.read_xes(str(LPlus_LogFile), variant="rustxes")
    if rules !=[]:
        activities = set(activities)
        lookup_table_path = Path("assets") / "lookup_table.csv"
    
        # if rules_path=="":
        #     rules= {}
        #     activities = {}
        # else:
        #     rules, activities = rules_from_json(str(rules_path))
        rules_proccessed, absence_list = preprocess(rules, dim, abs_thr)

        event_log_xes = pm4py.filter_event_attribute_values(
            event_log_xes,
            attribute_key="concept:name",  # Default attribute for activity names in XES
            values=absence_list,
            retain=False  # Retain=False means we exclude the specified activities
        )

        print('rules are preprocessed')

        lookup_table = pd.read_csv(lookup_table_path, sep=';', index_col=0)
        lookup_table = lookup_table.map(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

        S_mapped = declare_processing.assign_alphabet(activities)

        print('conversion of rules to DFAs started')
        dfa_list, total_sup, total_conf = dfa_list_generator(rules_proccessed, S_mapped)
        print('conversion of rules to DFAs ended')

        print(f"_______________ support is {support}_____________________")
        print('process discovery started')
        start = time.time()
        net, initial_marking, final_marking = apply_bi(Lp=event_log_xes, sup=support, ratio=0, size_par=1,
                                                                 rules=(rules_proccessed, lookup_table, dim))
        end = time.time()
        print(end - start)
        print('process discovery ended')

        pm4py.write_pnml(net, initial_marking, final_marking, os.path.join(r"output_files", "model.pnml"))



        print('model_checking started')
        rg = construct_reachability_graph(net, initial_marking, use_trans_name=False, parameters=None)
        aa = declare_processing.reachability2NFA(rg, activities)
        model_dfa = DFA.from_nfa(aa)
        cond = []
        sup_cost = 0
        conf_cost = 0
        for dfa in dfa_list:
            constraint_dfa_complement = dfa[1].complement()
            intersection_dfa = model_dfa.intersection(constraint_dfa_complement)
            if not intersection_dfa.isempty():
                sup_cost += dfa[2]
                conf_cost += dfa[2]
            cond.append((dfa[0], intersection_dfa.isempty(), dfa[2], dfa[3]))
        print('model_checking ended')

        report = {}
        report['time'] = end - start
        report['N.rules'] = len(dfa_list)
        report['N.dev'] = len([(x[0][0], x[0][1]) for x in cond if x[1] == False])
        if total_sup!=0:
            report['support_cost'] = round(sup_cost / total_sup, 2)
        else:
            report['support_cost'] = 0
        if total_conf != 0:
            report['confidence_cost'] = round(conf_cost / total_conf, 2)
            report['confidence_cost'] = 0
        report['dev_list'] = [(x[0][0], x[0][1], round(x[2], 2), round(x[3], 2)) for x in cond if x[1] == False]

        with open(os.path.join(r"output_files/", "stats.json"), "w") as json_file:
            json.dump(report, json_file, indent=4)

        print('The report is generated')
    else:
        net, initial_marking, final_marking = apply_bi(Lp=event_log_xes, sup=support, ratio=0, size_par=1,
                                                       rules=([], [], ""))


    parameters = {pn_visualizer.Variants.WO_DECORATION.value.Parameters.FORMAT: "png"}
    gviz = pn_visualizer.apply(net, initial_marking, final_marking, parameters=parameters)
    gviz.attr("graph", bgcolor="transparent")
    return gviz



def run_IMr_norule(LPlus_LogFile,support):

    event_log_xes = pm4py.read_xes(LPlus_LogFile, variant="rustxes")

    print(f"_______________ support is {support}_____________________")
    print('process discovery started')
    start = time.time()
    net, initial_marking, final_marking = apply_bi(Lp=event_log_xes, sup=support, ratio=0, size_par=1,
                                                             rules=({}, {}))
    end = time.time()
    print(end - start)
    print('process discovery ended')

    pm4py.write_pnml(net, initial_marking, final_marking, os.path.join(r"output_files", "model.pnml"))



    return net, initial_marking, final_marking
