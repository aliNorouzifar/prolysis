from prolysis.util.functions import add_SE, select_submatrix, bfs_descendants, bfs_ancestors, get_edge_weight_numba, n_edges_numba, convert_activities_to_array
from prolysis.discovery.candidate_search.is_allowed_2 import is_allowed_multi as is_allowed
import numpy as np
import statistics

def cut_thr(reserve,ind,precentage):
    values = [x[ind] for x in reserve]  # Extract the 4th column
    q1 = np.percentile(values, precentage)  # Compute the 25th percentile (Q1)
    return q1

def last_quartile(reserve):
    values = [x[4] for x in reserve]  # Extract the 4th column
    q3 = np.percentile(values, 75)  # Compute the 25th percentile (Q1)
    return q3

def find_possible_partitions(rules, st_net, en_net,fp, adj_matrix,nodes_order,mapping,adj_dict,surv_rate, actsP):
    # activity_list = set(net.nodes) - {'start', 'end'}
    
    # activity_list = set(mapping.keys()) - {'start', 'end'}
    activity_list = actsP - {'start', 'end'}
    queue = [(set(), {'start'})]
    visited = set()
    reserve = []

    nodes_order_np = np.array([mapping[s] for s in nodes_order], dtype=np.int16)

    block = False

    lowest_cost = 10000

    if get_edge_weight_numba(adj_matrix, mapping['start'], mapping['end']) > 0:
        if len(rules[0][0]) > 0 or len(rules[0][1]) > 0:
            na, block, penalty,penalty_des,lowest_cost, rules_dev_des = is_allowed(activity_list, set(), rules,{'exc_tau'},lowest_cost)
        else:
            na = set()
            penalty = {'seq': 0, 'exc': 0, 'par': 0, 'loop': 0, 'loop_tau': 0, 'exc_tau':0}
            penalty_des = {'seq': 0, 'exc': 0, 'par': 0, 'loop': 0, 'loop_tau': 0, 'exc_tau': 0}
            rules_dev_des = {'seq': [], 'exc': [], 'par': [], 'loop': [], 'loop_tau': [], 'exc_tau': []}
        possible_cuts = {'exc_tau'} - na
        reserve += [(activity_list, set(), ct, penalty[ct],penalty_des[ct], rules_dev_des[ct]) for ct in possible_cuts]


    if (n_edges_numba(adj_matrix, convert_activities_to_array(en_net,mapping), convert_activities_to_array(st_net,mapping)) > 0) and st_net==set([x for x in st_net if (x,x) in fp.edges() and fp[x][x]['weight'] > 0]):
        if len(rules[0][0]) > 0 or len(rules[0][1]) > 0:
            na, block, penalty,penalty_des,lowest_cost, rules_dev_des = is_allowed(activity_list, set(), rules,{'loop_tau'},lowest_cost)
        else:
            na = set()
            penalty = {'seq': 0, 'exc': 0, 'par': 0, 'loop': 0, 'loop_tau': 0, 'exc_tau':0}
            penalty_des = {'seq': 0, 'exc': 0, 'par': 0, 'loop': 0, 'loop_tau': 0, 'exc_tau': 0}
            rules_dev_des = {'seq': [], 'exc': [], 'par': [], 'loop': [], 'loop_tau': [], 'exc_tau': []}
        possible_cuts = {'loop_tau'} - na
        reserve += [(activity_list, set(), ct, penalty[ct],penalty_des[ct], rules_dev_des[ct]) for ct in possible_cuts]


    while len(queue) != 0:
        current = queue.pop()
        current_set, current_adj = current

        for x in current_adj:
            new_state = current_set | {x}
            new_state = add_SE(new_state, st_net, en_net)

            if frozenset(new_state) not in visited:
                new_adj = (current_adj | adj_dict[x]) - new_state

                visited.add(frozenset(new_state))

                B = activity_list - new_state

                if (len(B) == 0) or (len(B) == len(activity_list)):
                    queue.append((new_state, new_adj))
                    continue

                B = add_SE(B,st_net,en_net)
                new_state_np = np.array([mapping[s] for s in new_state], dtype=np.int16)
                B_np = np.array([mapping[s] for s in B], dtype=np.int16)

                adj_A = select_submatrix(adj_matrix, nodes_order_np, new_state_np)
                adj_B = select_submatrix(adj_matrix, nodes_order_np, B_np)


                # 'start' ~> netB
                if 'start' in B:
                    descendants = bfs_descendants(adj_B, np.where(B_np==mapping['start'])[0][0])
                    start2B = (descendants.shape[0] == (B_np.shape[0]-1))
                else:
                    start2B = False
                # 'end' ~> netA
                if 'end' in new_state:
                    ancestors = bfs_ancestors(adj_A, np.where(new_state_np==mapping['end'])[0][0])
                    A2end = (ancestors.shape[0] == (new_state_np.shape[0]-1))
                else:
                    A2end = False

                # 'end' ~> netB
                if 'end' in B:
                    ancestors = bfs_ancestors(adj_B, np.where(B_np==mapping['end'])[0][0])
                    B2end = (ancestors.shape[0] == (B_np.shape[0] - 1))
                else:
                    B2end = False


                possible_cuts = set()
                if B2end:
                    possible_cuts.add("seq")
                    if A2end:
                        if n_edges_numba(adj_matrix, convert_activities_to_array(B-{'start','end'}, mapping), convert_activities_to_array(new_state & st_net, mapping)) != 0 and n_edges_numba(adj_matrix,convert_activities_to_array(new_state & en_net,mapping), convert_activities_to_array(B-{'start','end'},mapping)) != 0:
                                possible_cuts.add("loop")
                        if start2B and (B not in visited):
                            possible_cuts.update(["exc", "par"])
                elif A2end:
                    if n_edges_numba(adj_matrix, convert_activities_to_array(B,mapping), convert_activities_to_array(new_state & st_net,mapping)) != 0 and n_edges_numba(adj_matrix, convert_activities_to_array(new_state & en_net,mapping),convert_activities_to_array(B,mapping)) != 0:
                            possible_cuts.add("loop")

                if possible_cuts:
                    if len(rules[0][0]) > 0 or len(rules[0][1]) > 0:
                        na, block, penalty,penalty_des,lowest_cost, rules_dev_des = is_allowed(new_state, B, rules,possible_cuts,lowest_cost)
                    else:
                        penalty = {'seq': 0, 'exc': 0, 'par': 0, 'loop': 0, 'loop_tau': 0, 'exc_tau': 0}
                        penalty_des = {'seq': 0, 'exc': 0, 'par': 0, 'loop': 0, 'loop_tau': 0, 'exc_tau': 0}
                        rules_dev_des = {'seq': [], 'exc': [], 'par': [], 'loop': [], 'loop_tau': [], 'exc_tau': []}
                    reserve += [(new_state, B, ct, penalty[ct],penalty_des[ct], rules_dev_des[ct]) for ct in possible_cuts]


                if not block:
                    queue.append((new_state, new_adj))
    # print(penalty)
    # min_value = min(reserve, key=lambda x: x[3])[3]
    if surv_rate!= "skip":
        print(f"reserve: {len(reserve)}")
        min_value = cut_thr(reserve,3,surv_rate)
        # buff_thr1 = 0.001 * sum([x['support'] for x in rules[0][0] if x['support'] > 0])
        # print(buff_thr1)
        buff_thr1 = 0

        # print(max([tup[3] for tup in reserve]))

        survived_l1 = [tup for tup in reserve if tup[3] <= min_value+buff_thr1]
        print(f"survived: {len(survived_l1)}")

        min_value2 = cut_thr(survived_l1,4,surv_rate)
        # buff_thr2 = 0.01 * (sum([x['desirability'] for x in rules[0][1] if x['desirability'] > 0])-sum([x['desirability'] for x in rules[0][1] if x['desirability'] < 0]))
        buff_thr2 = 0
        # print(buff_thr2)
        # print(max([tup[4] for tup in survived_l1]))
        # min_value2 = cut_thr(survived_l1, 4, surv_rate)

        result = [(tup[0],tup[1],tup[2]) for tup in survived_l1 if tup[4] <= min_value2+buff_thr2]

        print(f"final: {len(result)}")

        # print(f"reserve: {len(reserve)}")
        # k = 1
        #
        # mean_val_1 = statistics.mean([tup[3] for tup in reserve])
        # if len(reserve) !=1:
        #     std_val_1 = statistics.stdev([tup[3] for tup in reserve])
        # else:
        #     std_val_1 = 0
        # thr_1 = (std_val_1 + (k*mean_val_1))
        #
        # survived_l1 = [tup for tup in reserve if tup[3] <= thr_1]
        # print(f"survived: {len(survived_l1)}")
        #
        # mean_val_2 = statistics.mean([tup[4] for tup in survived_l1])
        # if len(survived_l1)!=1:
        #     std_val_2 = statistics.stdev([tup[4] for tup in survived_l1])
        # else:
        #     std_val_2 = 0
        # thr_2 = (mean_val_2 + (k * std_val_2))
        # # min_value2 = max(0, cut_thr(survived_l1, 4, surv_rate))
        # # min_value2 = cut_thr(survived_l1, 4, surv_rate)
        #
        # result = [(tup[0], tup[1], tup[2]) for tup in survived_l1 if tup[4] <= thr_2]
        #
        # print(f"final: {len(result)}")
    else:

        min_value = min(reserve, key=lambda x: x[3])[3]
        print(min_value)
        result = [(tup[0], tup[1], tup[2]) for tup in reserve if tup[3] == min_value]
    return result,reserve
