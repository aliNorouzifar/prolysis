from prolysis.util.functions import sum_out_degree_numba, sum_out_degree_single, get_edge_weight_numba, n_edges_numba
from numba import njit
import numpy as np



# @njit
def cost_seq(adj, A_indices, B_indices, sup, flow,min_cost,flag):
    # Compute total out-degree sums
    sum_out_degrees_A = sum_out_degree_numba(adj, A_indices)
    sum_out_degrees_B = sum_out_degree_numba(adj, B_indices)
    sum_out_degrees_total = sum_out_degrees_A + sum_out_degrees_B

    # Preallocate fixed-size storage for efficiency
    dev = 0
    mis = 0

    # Compute c1 and c2
    for i, x_idx in enumerate(A_indices):
        for j, y_idx in enumerate(B_indices):
            dev += get_edge_weight_numba(adj, y_idx, x_idx)

            if flag and dev> min_cost:
                return min_cost, dev, mis, False

            out_degree_x = sum_out_degree_single(adj, x_idx)
            out_degree_y = sum_out_degree_single(adj, y_idx)

            if sum_out_degrees_total > 0:
                mis += max(0, out_degree_x * sup * (out_degree_y / sum_out_degrees_total) -
                               get_edge_weight_numba(flow, x_idx, y_idx))
                if flag and dev > min_cost:
                    return min_cost, dev, mis, False
    if flag and dev < min_cost:
        min_cost = dev+mis

    return min_cost, dev, mis, True


# @njit
def cost_exc(adj, A_indices, B_indices,min_cost, flag):
    # Preallocate fixed-size storage for efficiency
    dev = 0

    # Compute c1
    for i, x_idx in enumerate(A_indices):
        for j, y_idx in enumerate(B_indices):
            dev += get_edge_weight_numba(adj, x_idx, y_idx)
            dev += get_edge_weight_numba(adj, y_idx, x_idx)
            if flag and dev > min_cost:
                return min_cost, dev, False
    if flag and dev < min_cost:
        min_cost = dev

    return min_cost, dev, True



def cost_exc_tau(adj, sup,start_index, end_index):
    # c = {}
    # if 'start' in net.nodes():
    mis = max(0, sup * sum_out_degree_single(adj, start_index) - get_edge_weight_numba(adj, start_index, end_index))
    return mis



# @njit
def cost_par(adj, A_indices, B_indices, sup,min_cost, flag):
    # Compute total out-degree sums
    sum_out_degrees_A = sum_out_degree_numba(adj, A_indices)
    sum_out_degrees_B = sum_out_degree_numba(adj, B_indices)
    sum_out_degrees_total = sum_out_degrees_A + sum_out_degrees_B

    # Preallocate fixed-size storage for efficiency
    mis = 0

    # Compute c1 and c2
    for i, a_idx in enumerate(A_indices):
        out_degree_a = sum_out_degree_single(adj, a_idx)
        for j, b_idx in enumerate(B_indices):
            out_degree_b = sum_out_degree_single(adj, b_idx)

            if sum_out_degrees_total > 0:
                mis += max(0, (out_degree_a * sup * out_degree_b) / sum_out_degrees_total -
                               get_edge_weight_numba(adj, a_idx, b_idx))
                mis += max(0, (out_degree_b * sup * out_degree_a) / sum_out_degrees_total -
                               get_edge_weight_numba(adj, b_idx, a_idx))
                if flag and mis > min_cost:
                    return min_cost, mis, False
    if flag and mis < min_cost:
        min_cost = mis

    return min_cost, mis, True


# @njit
def cost_loop(adj, A_indices, B_indices, sup, start_A_indices, end_A_indices, input_B_indices, output_B_indices, start_index, end_index,min_cost,flag):
    M_P = max(n_edges_numba(adj, output_B_indices, start_A_indices),
              n_edges_numba(adj, end_A_indices, input_B_indices))

    dev_cost = 0
    mis_cost = 0

    # Compute c1, c2, c3
    for x_idx in B_indices:
        dev_cost += get_edge_weight_numba(adj, start_index, x_idx) + get_edge_weight_numba(adj, x_idx, end_index)
        if flag and dev_cost > min_cost:
            return min_cost, dev_cost, mis_cost, False

        for y_idx in A_indices:
            if y_idx not in end_A_indices:
                dev_cost += get_edge_weight_numba(adj, y_idx, x_idx)
                if flag and dev_cost > min_cost:
                    return min_cost, dev_cost, mis_cost, False

            if y_idx not in start_A_indices:
                dev_cost += get_edge_weight_numba(adj, x_idx, y_idx)

            if flag and dev_cost > min_cost:
                    return min_cost, dev_cost, mis_cost, False

    # Compute c4 (Avoid np.array inside loops)
    if n_edges_numba(adj, output_B_indices, start_A_indices) > 0:
        for a_idx in start_A_indices:
            for b_idx in output_B_indices:
                start_edges = n_edges_numba(adj, np.array([start_index]), np.array([a_idx]))
                total_start_edges = n_edges_numba(adj, np.array([start_index]), start_A_indices)
                output_edges = n_edges_numba(adj, np.array([b_idx]), start_A_indices)
                total_output_edges = n_edges_numba(adj, output_B_indices, start_A_indices)

                if total_start_edges > 0 and total_output_edges > 0:
                    mis_cost += max(0, M_P * sup * (start_edges / total_start_edges) *
                                    (output_edges / total_output_edges) -
                                    get_edge_weight_numba(adj, b_idx, a_idx))

                if flag and dev_cost + mis_cost > min_cost:
                    return min_cost, dev_cost, mis_cost, False
    else:
        mis_cost += M_P * sup

    # Compute c5 (Avoid np.array inside loops)
    if n_edges_numba(adj, end_A_indices, input_B_indices) > 0:
        for a_idx in end_A_indices:
            for b_idx in input_B_indices:
                end_edges = n_edges_numba(adj, np.array([a_idx]), np.array([end_index]))
                total_end_edges = n_edges_numba(adj, end_A_indices, np.array([end_index]))
                input_edges = n_edges_numba(adj, end_A_indices, np.array([b_idx]))
                total_input_edges = n_edges_numba(adj, end_A_indices, input_B_indices)

                if total_end_edges > 0 and total_input_edges > 0:
                    mis_cost += max(0, M_P * sup * (end_edges / total_end_edges) *
                                    (input_edges / total_input_edges) -
                                    get_edge_weight_numba(adj, a_idx, b_idx))

                if flag and dev_cost + mis_cost > min_cost:
                    return min_cost, dev_cost, mis_cost, False
    else:
        mis_cost += M_P * sup

    if flag and dev_cost + mis_cost < min_cost:
        min_cost = dev_cost + mis_cost
    return min_cost,dev_cost, mis_cost, True




def cost_loop_tau(adj, sup, start_A_indices, end_A_indices,start_index, end_index):
    M_P = n_edges_numba(adj, end_A_indices, start_A_indices)


    start_sum = n_edges_numba(adj, np.array([start_index]), start_A_indices)
    end_sum = n_edges_numba(adj, end_A_indices, np.array([end_index]))
    mis_cost = 0

    for a_idx in start_A_indices:
        for b_idx in end_A_indices:
            if adj[b_idx,a_idx]>0:
                mis_cost += max(0, M_P * sup * (get_edge_weight_numba(adj,start_index,a_idx)/ start_sum)*
                                (get_edge_weight_numba(adj, b_idx, end_index)/ end_sum) -
                                get_edge_weight_numba(adj, b_idx, a_idx))
            else:
                mis_cost +=  M_P * sup * (get_edge_weight_numba(adj, start_index, a_idx) / start_sum) * (get_edge_weight_numba(adj, b_idx, end_index) / end_sum)
    return mis_cost


def overal_cost(costP,costM,ratio, size_par):
    ov_c= costP - ratio * size_par * costM
    return ov_c