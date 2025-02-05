import prolysis.discovery.split_functions.split as split
from prolysis.discovery.candidate_search.search import find_possible_partitions
from prolysis.discovery.base_case.check import check_base_case
from prolysis.discovery.cut_quality.cost_functions import cost_functions
from prolysis.util.functions import generate_nx_graph_from_log, generate_nx_indirect_graph_from_log, get_input_output_B_indices
from prolysis.util.functions import convert_activities_to_array as act2ary
from prolysis.discovery.cut_quality.cost_functions.cost_functions import overal_cost
import networkx as nx
import numpy as np
from collections import Counter


class SubtreePlain(object):
    def __init__(self, logp,logm, rec_depth, noise_threshold= None,
                   parameters=None, sup= None, ratio = None, size_par = None, rules = None):

        nt = 0.000 * sum(logp.values())
        self.rec_depth = rec_depth
        self.noise_threshold = noise_threshold
        self.log = logp

        if sum(self.log.values())==0:
            self.activitiesP = set()
            self.start_activitiesP = {}
            self.end_activitiesP = {}
        else:
            self.activitiesP = set(y for x in logp for y in x if x!=())
            self.start_activities = set(x[0] for x in logp if x!=())
            self.end_activities = set(x[-1] for x in logp if x!=())


        if not logm:
            self.ignore_M = True
            self.activitiesM = set()
            self.start_activitiesM = {}
            self.end_activitiesM = {}
            self.size_adj = 1
        else:
            self.ignore_M = False
            self.logM = logm
            self.activitiesM = set(y for x in logm for y in x if x!=())
            self.start_activitiesM = set(x[0] for x in logm if x!=())
            self.end_activitiesM = set(x[-1] for x in logm if x!=())


        self.netP = generate_nx_graph_from_log(self.log, nt, self.activitiesM)
        self.all_acts = self.activitiesP.union(self.activitiesM)
        self.nodes_order = list(self.all_acts - {'start', 'end'})
        if 'start' in self.netP.nodes and 'end' in self.netP.nodes:
            self.nodes_order.append('start')
            self.nodes_order.append('end')
        self.mapping = {i: activity for i, activity in enumerate(self.nodes_order)}
        self.mapping_rev = {v: k for k, v in self.mapping.items()}
        self.adj_matrixP = nx.to_numpy_array(self.netP, nodelist=self.nodes_order)  # Ensure correct ordering
        self.adj_dict = {n: set(self.netP.neighbors(n)) for n in self.netP.nodes}

        if logm:
            self.netM = generate_nx_graph_from_log(self.logM,nt,self.activitiesP)
            self.adj_matrixM = nx.to_numpy_array(self.netM, nodelist=self.nodes_order)
            self.size_adj = sum(self.log.values())/sum(self.logM.values())


        self.original_log = logp
        self.parameters = parameters
        self.detected_cut = None
        self.children = []
        self.detect_cut(parameters=parameters, sup= sup, ratio = ratio, size_par = size_par, rules = rules)


    def detect_cut(self, parameters=None, sup= None, ratio = None, size_par = None, rules = None):

        if parameters is None:
            parameters = {}

        # check base cases:
        isbase, cut = check_base_case(rules , sup, ratio, self)

        if isbase==False:
            fP = generate_nx_indirect_graph_from_log(self.log,self.activitiesM)
            fp_adj = nx.to_numpy_array(fP, nodelist=self.nodes_order[:-2])
            if not self.ignore_M:
                fM = generate_nx_indirect_graph_from_log(self.logM,self.activitiesP)
                fm_adj = nx.to_numpy_array(fM, nodelist=self.nodes_order[:-2])

            possible_partitions, reserve2 = find_possible_partitions(rules, self.start_activities,
                                                                    self.end_activities, fP, self.adj_matrixP, self.nodes_order,self.mapping_rev,self.adj_dict)

            # print(len(possible_partitions))
            cut = []

            ratio_backup = ratio

            min_cost = 65534
            for pp in possible_partitions:
                A = pp[0] - {'start', 'end'}
                A_np = np.array([self.mapping_rev[s] for s in A], dtype=np.int16)
                B = pp[1] - {'start', 'end'}
                B_np = np.array([self.mapping_rev[s] for s in B], dtype=np.int16)

                start_A_P = self.start_activities & A
                end_A_P = self.end_activities & A
                input_B_P, output_B_P = get_input_output_B_indices(self.adj_matrixP,B_np)

                if not self.ignore_M:
                    start_A_M = self.start_activitiesM & A
                    end_A_M = self.end_activitiesM & A
                    input_B_M, output_B_M = get_input_output_B_indices(self.adj_matrixM, B_np)
                    if sum(self.logM.values()) == 0:
                        ratio = 0
                    else:
                        ratio = ratio_backup

                type = pp[2]


                #####################################################################
                # seq check
                if type=="seq":
                    min_cost, dev, mis, valid = cost_functions.cost_seq(self.adj_matrixP, A_np,
                    B_np, sup, fp_adj, min_cost, self.ignore_M)
                    cost_seq_P_dict = {'missing':mis, 'deviating':dev}
                    cost_seq_P = mis+dev
                    if not self.ignore_M:
                        min_cost, dev, mis, valid = cost_functions.cost_seq(self.adj_matrixM, A_np,
                                                                  B_np, sup, fm_adj, min_cost, self.ignore_M)
                        cost_seq_M_dict = {'missing':mis, 'deviating':dev}
                        cost_seq_M = mis+dev
                    else:
                        cost_seq_M = 0
                        cost_seq_M_dict = {}
                    if valid:
                        cut.append(((A, B), 'seq', cost_seq_P, cost_seq_M, overal_cost(cost_seq_P,cost_seq_M,ratio,self.size_adj),cost_seq_P_dict,cost_seq_M_dict))
                #####################################################################
                # xor check
                if type=="exc":
                    min_cost, dev, valid = cost_functions.cost_exc(self.adj_matrixP, A_np,
                     B_np,min_cost, self.ignore_M)
                    cost_exc_P_dict = {'missing':0, 'deviating':dev}
                    cost_exc_P = dev

                    if not self.ignore_M:
                        min_cost, dev, valid = cost_functions.cost_exc(self.adj_matrixM,A_np,
                                                                  B_np, min_cost, self.ignore_M)
                        cost_exc_M_dict = {'missing':0, 'deviating':dev}
                        cost_exc_M = dev
                    else:
                        cost_exc_M_dict = {}
                        cost_exc_M = 0
                    if valid:
                        cut.append(((A, B), 'exc', cost_exc_P, cost_exc_M, overal_cost(cost_exc_P,cost_exc_M,ratio,self.size_adj),cost_exc_P_dict,cost_exc_M_dict))
                #####################################################################
                # xor-tau check
                if type=="exc_tau":
                    mis = cost_functions.cost_exc_tau(self.adj_matrixP,sup,self.mapping_rev['start'],self.mapping_rev['end'])
                    cost_exc_tau_P_dict = {'missing':mis, 'deviating':0}
                    cost_exc_tau_P = mis
                    if not self.ignore_M:
                        mis = cost_functions.cost_exc_tau(self.adj_matrixM, sup,self.mapping_rev['start'],self.mapping_rev['end'])
                        cost_exc_tau_M_dict = {'missing':mis, 'deviating':0}
                        cost_exc_tau_M = mis
                    else:
                        cost_exc_tau_M_dict = {}
                        cost_exc_tau_M = 0
                    cut.append(((A.union(B), set()), 'exc_tau',cost_exc_tau_P , cost_exc_tau_M, overal_cost(cost_exc_tau_P,cost_exc_tau_M,ratio,self.size_adj),cost_exc_tau_P_dict,cost_exc_tau_M_dict))
                #####################################################################
                # parallel check
                if type=="par":
                    min_cost, mis, valid = cost_functions.cost_par(self.adj_matrixP, A_np,
                     B_np, sup, min_cost,self.ignore_M)
                    cost_par_P_dict = {'missing':mis, 'deviating':0}
                    cost_par_P = mis
                    if not self.ignore_M:
                        min_cost, mis, valid = cost_functions.cost_par(self.adj_matrixM, A_np,
                                                                  B_np, sup, min_cost,self.ignore_M)
                        cost_par_M_dict = {'missing':mis, 'deviating':0}
                        cost_par_M = mis
                    else:
                        cost_par_M_dict = {}
                        cost_par_M = 0
                    if valid:
                        cut.append(((A, B), 'par', cost_par_P, cost_par_M, overal_cost(cost_par_P,cost_par_M,ratio,self.size_adj),cost_par_P_dict,cost_par_M_dict))
                #####################################################################
                # loop check
                if type=="loop":
                    min_cost,dev, mis, valid = cost_functions.cost_loop(self.adj_matrixP, A_np, B_np, sup,
                                                                   act2ary(start_A_P, self.mapping_rev), act2ary(end_A_P, self.mapping_rev), input_B_P,
                                                                   output_B_P,self.mapping_rev['start'],self.mapping_rev['end'],min_cost,self.ignore_M)

                    cost_loop_P_dict = {'missing':mis, 'deviating':dev}
                    cost_loop_P = mis+dev
                    if not self.ignore_M:
                        min_cost,dev, mis, valid = cost_functions.cost_loop(self.adj_matrixM,
                                                                       A_np,
                                                                       B_np,
                                                                       sup,
                                                                       act2ary(start_A_M,
                                                                                                   self.mapping_rev),
                                                                       act2ary(end_A_M,
                                                                                                   self.mapping_rev),
                                                                       input_B_M,
                                                                       output_B_M,
                                                                       self.mapping_rev['start'], self.mapping_rev['end'],min_cost,self.ignore_M)
                        cost_loop_M_dict = {'missing':mis, 'deviating':dev}
                        cost_loop_M = mis+dev
                    else:
                        cost_loop_M_dict = {}
                        cost_loop_M = 0

                    if valid:
                        cut.append(((A, B), 'loop', cost_loop_P, cost_loop_M, overal_cost(cost_loop_P,cost_loop_M,ratio,self.size_adj),cost_loop_P_dict,cost_loop_M_dict))

                if type=="loop_tau":
                    mis = cost_functions.cost_loop_tau(self.adj_matrixP, sup, act2ary(start_A_P, self.mapping_rev),act2ary(end_A_P, self.mapping_rev),self.mapping_rev['start'],self.mapping_rev['end'])
                    cost_loop_P_dict = {'missing':mis, 'deviating': 0 }
                    cost_loop_P = mis
                    if not self.ignore_M:
                        mis = cost_functions.cost_loop_tau(self.adj_matrixM, sup, act2ary(start_A_M, self.mapping_rev),act2ary(end_A_M, self.mapping_rev),self.mapping_rev['start'],self.mapping_rev['end'])
                        cost_loop_M_dict = {'missing':mis, 'deviating': 0 }
                        cost_loop_M = mis
                    else:
                        cost_loop_M_dict = {}
                        cost_loop_M = 0
                    cut.append(((self.start_activities, self.end_activities), 'loop_tau', cost_loop_P,cost_loop_M, overal_cost(cost_loop_P,cost_loop_M,ratio,self.size_adj), cost_loop_P_dict,cost_loop_M_dict))


            if not cut:
                print("no good cut exists")
            sorted_cuts = sorted(cut, key=lambda x: (x[4], x[2],['exc_tau','exc','seq','par','loop','loop_tau'].index(x[1]), abs(len(x[0][0])-len(x[0][1]))))
            cut = sorted_cuts[0]
            # print(len(sorted_cuts))

        # print(cut[:-2])

        map_cut_op = {'par': 'parallel', 'seq': 'sequential', 'exc': 'xor', 'exc_tau':'xor', 'exc2': 'xor',
                      'loop': 'loopCut', 'loop1': 'loopCut', 'loop_tau': 'loopCut'}



        if cut[1] in map_cut_op.keys():
            self.detected_cut = map_cut_op[cut[1]]
            LAP, LBP = split.split(cut[1], [cut[0][0], cut[0][1]], self.log)
            if not self.ignore_M:
                LAM, LBM = split.split(cut[1], [cut[0][0], cut[0][1]], self.logM)
            else:
                LAM = Counter()
                LBM = Counter()
            new_logs = [[LAP, LAM], [LBP, LBM]]
            for l in new_logs:
                self.children.append(
                    SubtreePlain(l[0], l[1],
                                 self.rec_depth + 1,
                                 noise_threshold=self.noise_threshold,
                                 parameters=parameters, sup=sup, ratio=ratio, size_par=size_par, rules=rules))
        elif cut[1]!='single_activity' and cut[1]!='empty_log':
            print('It should not happen, if you see this error there could be a bug in the code!')

