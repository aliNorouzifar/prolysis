import networkx as nx
import json
import numpy as np
from numba import njit
from numba.typed import List


# @njit
def select_submatrix(adj_matrix, nodes_order, included_nodes):
    # Select submatrix
    new_adj_matrix = adj_matrix[included_nodes[:, None], included_nodes]

    return new_adj_matrix



@njit
def bfs_descendants(adj_matrix, start_index):
    n = adj_matrix.shape[0]

    # Using Numba-typed list for performance
    queue = List()
    queue.append(start_index)

    visited = np.zeros(n, dtype=np.uint8)
    visited[start_index] = 1

    # Use a fixed-size NumPy array (worst-case: all nodes reachable)
    descendants = np.full(n, -1, dtype=np.int32)  # Initialize with -1 (invalid index)
    count = 0  # Track number of valid descendants

    while len(queue) > 0:
        node = queue.pop(0)
        for neighbor in range(n):
            if adj_matrix[node, neighbor] > 0 and not visited[neighbor]:
                queue.append(neighbor)
                visited[neighbor] = 1
                descendants[count] = neighbor
                count += 1  # Move to next slot

    return descendants[:count]  # Return only valid indices


@njit
def bfs_ancestors(adj_matrix, start_index):
    return bfs_descendants(adj_matrix.T, start_index)  # BFS on reversed edges




def n_edges(net, S, T):
    edges_reweight = list(nx.edge_boundary(net, S, T, data='weight', default=1))
    return sum(weight for u, v, weight in edges_reweight if (u in S and v in T))


# def n_edges2(adj_matrix, index_to_activity, S, T):
#     # Convert activity names to row indices
#     S_indices = {i for i, activity in index_to_activity.items() if activity in S}
#     T_indices = {i for i, activity in index_to_activity.items() if activity in T}
#
#     # Sum up the edges from S to T in the adjacency matrix
#     edge_count = sum(adj_matrix[i, j] for i in S_indices for j in T_indices)
#
#     return edge_count

# @njit
def sum_out_degree_numba(adj_matrix, A_indices):
    total_degree = 0.0
    for i in A_indices:
        total_degree += adj_matrix[i, :].sum()
    return total_degree


# @njit
def sum_out_degree_single(adj_matrix, A_index):
    return adj_matrix[A_index, :].sum()

# @njit
def get_edge_weight_numba(adj, node1_idx, node2_idx):
    return adj[node1_idx, node2_idx]

# @njit
def n_edges_numba(adj, S_indices, T_indices):
    count = 0
    for s in S_indices:
        for t in T_indices:
            count += adj[s, t]
    return count


def get_edge_weight(graph, node1, node2):
    edge_data = graph.get_edge_data(node1, node2)
    if edge_data is None:
        return 0
    else:
        return edge_data['weight']

def get_edge_weight2(adj_matrix, index_to_activity, node1, node2):
    # Convert node names to indices
    node1_index = next((i for i, act in index_to_activity.items() if act == node1), None)
    node2_index = next((i for i, act in index_to_activity.items() if act == node2), None)

    # If nodes are not found, return 0
    if node1_index is None or node2_index is None:
        return 0

    # Retrieve the edge weight from the adjacency matrix
    return adj_matrix[node1_index, node2_index]

def sum_out_degree(adj_matrix, index_to_activity, A):
    # Convert activity names in A to row indices
    A_indices = {i for i, activity in index_to_activity.items() if activity in A}

    # Sum over the selected rows to get total out-degree
    return sum(adj_matrix[i, :].sum() for i in A_indices)


def convert_activities_to_array(activity_set, mapping):
    x = np.array([mapping[activity] for activity in activity_set], dtype=np.int32)
    return x

def add_SE(nodes, st,en):
    if nodes & st:
        nodes.add('start')
    if nodes & en:
        nodes.add('end')
    return nodes

def out_degree_activity(adj_matrix, index_to_activity, activity):
    # Get the index of the activity
    activity_index = next((i for i, act in index_to_activity.items() if act == activity), None)

    # If the activity is not found, return 0
    if activity_index is None:
        return 0

    # Return the sum of the row corresponding to the activity (out-degree)
    return adj_matrix[activity_index, :].sum()


@njit
def get_input_output_B_indices(adj_matrix, B_indices):
    n = adj_matrix.shape[0]

    # Create a boolean mask for nodes in B
    B_mask = np.zeros(n, dtype=np.uint8)
    for i in range(B_indices.shape[0]):
        B_mask[B_indices[i]] = 1  # Mark nodes in B

    # Lists to store results
    input_B_list = []
    output_B_list = []

    for node in range(n):
        if B_mask[node]:  # Only check nodes in B
            has_incoming = False
            has_outgoing = False

            for neighbor in range(n):
                if not B_mask[neighbor]:  # Node is outside B
                    if adj_matrix[neighbor, node] > 0:  # Incoming edge
                        has_incoming = True
                    if adj_matrix[node, neighbor] > 0:  # Outgoing edge
                        has_outgoing = True

                # If both incoming & outgoing edges are found, break early
                if has_incoming and has_outgoing:
                    break

            if has_incoming:
                input_B_list.append(node)
            if has_outgoing:
                output_B_list.append(node)

    return np.array(input_B_list, dtype=np.int32), np.array(output_B_list, dtype=np.int32)

def aggregate_dictionaries(deviating, missing):
    aggregated_dict = {}

    for d in deviating:
        for key, value in d.items():
            if key in aggregated_dict:
                aggregated_dict[key]['deviating'] += value
            else:
                aggregated_dict[key] = {'deviating':value,'missing':0}

    for d in missing:
        for key, value in d.items():
            if key in aggregated_dict:
                aggregated_dict[key]['missing'] += value
            else:
                aggregated_dict[key] = {'deviating':0,'missing':value}
    return aggregated_dict

def aggregate_dictionaries2(deviating, missing):
    aggregated_dict = {'deviating':0,'missing':0}

    for d in deviating:
        aggregated_dict['deviating'] += np.sum(d)


    for d in missing:
        aggregated_dict['missing'] += np.sum(d)

    return aggregated_dict

def generate_nx_graph_from_log(log,nt, dummy_nodes):
    G = nx.DiGraph()
    for trace in log:
        tr_art = ('start',) + trace + ('end',)
        for i in range(len(tr_art) - 1):
            if G.has_edge(tr_art[i], tr_art[i + 1]):
                G[tr_art[i]][tr_art[i + 1]]['weight'] += log[trace]
            else:
                G.add_edge(tr_art[i], tr_art[i + 1], weight=log[trace])

    st = {}
    en = {}
    # Filter outgoing edges from start_node based on weight
    outgoing_edges = list(G.edges('start', data=True))
    for u, v, data in outgoing_edges:
        if data['weight'] < nt:
            G.remove_edge(u, v)
            if not nx.has_path(G, 'start', v):
               G.add_edge(u, v, weight=data['weight'])
               st[v] = data['weight']
        else:
            st[v] = data['weight']

    # Filter incoming edges to end_node based on frequency
    incoming_edges = list(G.in_edges('end', data=True))
    for u, v, data in incoming_edges:
        if data['weight'] < nt:
            G.remove_edge(u, v)
            if not nx.has_path(G, u, 'end'):
               G.add_edge(u, v, weight=data['weight'])
               en[u] = data['weight']
        else:
            en[u] = data['weight']
    for n in dummy_nodes:
        if n not in G.nodes:
            G.add_node(n)
    return G




def generate_nx_indirect_graph_from_log(log, dummy_nodes):
    # print("**************")
    G = nx.DiGraph()
    for trace in log:
        # print(trace)
        for i in range(0,len(trace)):
            if trace[i] not in G.nodes:
                G.add_node(trace[i])
            visited = set()
            for j in range(i+1,len(trace)):
                if trace[j] not in visited:
                    visited.add(trace[j])
                    if G.has_edge(trace[i], trace[j]):
                        G[trace[i]][trace[j]]['weight'] += log[trace]
                    else:
                        G.add_edge(trace[i], trace[j], weight=log[trace])
    for n in dummy_nodes:
        if n not in G.nodes:
            G.add_node(n)

    # print(G.nodes)
    return G

def read_append_write_json(file_path, new_item):
    # Read the current content of the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)

    # Append the new item to the list
    data.append(new_item)

    # Write the updated list back to the JSON file
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)

