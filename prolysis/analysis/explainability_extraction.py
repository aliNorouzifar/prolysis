import matplotlib
matplotlib.use('agg')
import pandas as pd
import pm4py
import shutil
import os
from prolysis.calls.minerful_calls import mine_minerful_for_declare_constraints, prune_constraints_minerful
import ruptures as rpt
import csv
from sklearn.metrics import mean_squared_error
from collections import defaultdict
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from scipy.spatial import KDTree
import json
from prolysis.util.redis_connection import redis_client
from prolysis.util.logging import log_command

linkage_method = 'ward'
linkage_metric = 'euclidean'

def export_constraints_per_cluster(constraints, constraints_json_path):
    dict_out = {}
    dict_out["constraints"] = []
    for constraint in constraints:
        temp_dict = {}
        temp_dict["template"] = constraint[0]
        temp_dict['parameters'] = []
        for par in constraint[1:3]:
            if par != '':
                temp_dict['parameters'].append([par.lstrip()])

        dict_out["constraints"].append(temp_dict)

    for key in range(len(dict_out['constraints'])):
        dict_out['constraints'][key]['template'] = dict_out['constraints'][key]['template'].rstrip().lstrip()

        dict_out['constraints'][key]['support'] = 0.99
        dict_out['constraints'][key]['confidence'] = 0.99

        for el in dict_out['constraints'][key]['parameters']:
            el[0] = el[0]

    with open(constraints_json_path, 'w') as fp:
        json.dump(dict_out, fp)
        fp.close()


def generate_features(w,kpi,n_bin):
    case_table = pd.read_csv("output_files/out.csv").sort_values(by=[kpi])
    ordered_case_ids = case_table['case_id']
    ordered_case_ids = ordered_case_ids.to_list()
    bin_size = round(len(case_table) / n_bin)
    event_table = pd.read_csv("output_files/out_event.csv")
    event_table['case:concept:name'] = event_table['case_id'].astype(str)
    event_table['concept:name'] = event_table['activity_name'].astype(str)
    event_table['time:timestamp'] = pd.to_datetime(event_table['timestamp'],format='mixed')
    event_table = event_table.sort_values('time:timestamp')
    event_table['case_id'] = pd.Categorical(event_table['case_id'], categories=ordered_case_ids, ordered=True)
    event_table = event_table.sort_values(by=['case_id', 'time:timestamp'], ascending=[True, True])
    event_table = event_table[['case:concept:name','concept:name','time:timestamp']]
    log_xes = pm4py.convert_to_event_log(event_table)
    pm4py.write_xes(log_xes, f"output_files/log_ordered.xes")

    window_size = 2 * w * bin_size
    sliding_window_size = bin_size
    mine_minerful_for_declare_constraints(window_size,sliding_window_size)

def find_all_supersets(graph,node, visited):
    for neighbor in graph[node]:
        if neighbor not in visited:
            visited.add(neighbor)
            find_all_supersets(graph,neighbor, visited)
    return visited

def correlation_calc(peaks,w,constraints,clusters_dict):
    n = len(peaks)
    peakmodif = [p-(w-1) for p in peaks]
    segments_sig = {}

    for i in range(n + 1):
        l = [0] * len(constraints[0])

        if i == 0:
            l[0:peakmodif[0]] = [1] * (peakmodif[0] - 0)

        elif i == n:
            l[peakmodif[n - 1]:] = [1] * (len(constraints[0]) - peakmodif[n - 1])

        else:
            l[peakmodif[i - 1]:peakmodif[i]] = [1] * (peakmodif[i] - peakmodif[i - 1])

        segments_sig[f'segment_{i}']=l




    corr_mat = []
    for seg in segments_sig.keys():
        target_array = np.array(segments_sig[seg])
        average_correlations = []

        for cluster in sorted([i for i in clusters_dict.keys()]):
            correlations = []
            for series in clusters_dict[cluster]:
                # Convert the current series to a NumPy array
                series_array = np.array(series[3:])

                # Calculate the correlation coefficient between target_series and this series
                correlation = np.corrcoef(target_array, series_array)[0, 1]

                # Append the correlation to the list
                correlations.append(correlation)
            average_correlation = np.mean(correlations)
            average_correlations.append(average_correlation)
        corr_mat.append(average_correlations)
    return corr_mat

def sort_by_closest_neigbour_HEADER(data):

    print('There were: ' + str(len(data)) + " values")

    # Convert data to numpy array for efficient operations
    data = np.array(data, dtype=object)

    # Initialize sorted data with the starting point (first point) and track remaining indices
    new_data = [data[0].tolist()]  # Start with the first point
    index_set = set(range(1, len(data)))  # Skip the first index, it's already in new_data

    # Track the original indices of points in the KDTree
    kd_tree_indices = list(index_set)
    kd_tree_data = data[kd_tree_indices]
    kd_tree = KDTree(kd_tree_data[:, 3:])

    while index_set:
        # Find nearest point in KDTree to the last point in `new_data`
        last_point = new_data[-1][3:]
        _, nearest_in_kd = kd_tree.query(last_point)

        # Map the KDTree index to the original data index
        min_ind = kd_tree_indices[nearest_in_kd]
        new_data.append(data[min_ind].tolist())


        # Update index_set and KDTree
        index_set.remove(min_ind)

        # Rebuild kd_tree for remaining points
        kd_tree_indices = list(index_set)
        kd_tree_data = data[kd_tree_indices]
        kd_tree = KDTree(kd_tree_data[:, 3:])

    return new_data

def import_minerful_constraints_timeseries_data(minerful_constraints_path,constraint_type_used):
    csvfile = open(minerful_constraints_path, 'r')
    csv_reader = csv.reader(csvfile, delimiter=';', quotechar='|')

    hea = next(csv_reader, None)
    hea2 = next(csv_reader, None)

    hea2 = hea2[2:]
    hea = hea[2:]

    header_output = list()

    for i in range(len(hea)):
        if i % 6 == 0:
            tem_h = [hea2[i][1:-1]]
            temp = hea[i]
            if temp[0] == '\'':
                temp = temp[1:]
            if temp[-1] == '\'':
                temp = temp[:-1]
            if temp[-1] == ')':
                temp = temp[:-1]
            # now we split the string
            name_of_constraint_end_index = temp.find('(')
            tem_h.append(temp[:name_of_constraint_end_index])
            temp = temp[name_of_constraint_end_index+1:]
            #find if we have two events or one
            separated_constraints_index = temp.find(', ')
            if not separated_constraints_index == -1:
                tem_h.append(temp[:separated_constraints_index])
                tem_h.append(temp[separated_constraints_index+1:])
            else:
                tem_h.append(temp)
                tem_h.append('')
        else:
            tem_h = [hea2[i][1:-1]] + tem_h[1:]

        header_output.append(tem_h)

    sequences = list()

    for i in range(len(hea)):
        sequences.append(list())

    corresponding_number_of_traces = []
    n_lines =0
    for r in csv_reader:
        corresponding_number_of_traces.append(r[:2])
        n_lines += 1
        counter = 0
        for i in range(len(r)):
            if counter > 1:
                sequences[i-2].append(float(r[i]))
            else:
                counter += 1


    constraints = []
    for i, j in zip(sequences, header_output):
        if j[0] == constraint_type_used:
            constraints.append(j[1:] + i)

    return constraints


def group_signals_by_type_activity(signals):
    # Create a nested defaultdict structure to handle the required dictionary format
    grouped_signals = defaultdict(lambda: defaultdict(set))

    for signal in signals:
        # Extract components
        signal_type = signal[0]  # The type of the signal
        activities = (signal[1], signal[2])  # Tuple of (activity1, activity2)
        signal_samples = tuple(signal[3:])  # Convert samples from position 3 onward to a tuple

        # Populate the dictionary
        if signal_type in {'RespondedExistence'}:
            grouped_signals[signal_samples][(activities[1].lstrip(), f' {activities[0]}')].add(
                'RespondedExistence_r')
            grouped_signals[signal_samples][activities].add(signal_type)
        else:
            grouped_signals[signal_samples][activities].add(signal_type)

    return grouped_signals

def prune_signals(theta_cvg):
    constraints_conf = import_minerful_constraints_timeseries_data(r"output_files/behavioral_signals.csv", 'Confidence')
    constraints_cov = import_minerful_constraints_timeseries_data(r"output_files/behavioral_signals.csv", 'Coverage')
    s_thr = theta_cvg * 100
    main_dim = constraints_conf
    filt_dim = constraints_cov

    not_include_templates = {'NotRespondedExistence',
                             'NotCoExistence',
                             'NotPrecedence',
                             'NotChainPrecedence',
                             'NotResponse',
                             'NotChainResponse',
                             'NotSuccession',
                             'NotChainSuccession'}

    filtered = []
    filtered_colored = []
    for i, c in enumerate(main_dim):
        filtered_sig = []
        filtered_colored_sig = []
        for j, v in enumerate(c):
            if j >= 3:
                if filt_dim[i][j] >= s_thr:
                    filtered_sig.append(v)
                    filtered_colored_sig.append(False)
                else:
                    filtered_sig.append(0)
                    if filt_dim[i][j] > 0:
                        filtered_colored_sig.append(True)
                    else:
                        filtered_colored_sig.append(False)
            else:
                filtered_sig.append(v)
                filtered_colored_sig.append(v)
        filtered.append(filtered_sig)
        filtered_colored.append(filtered_colored_sig)


    data = []
    data_color = []
    not_include = 0
    zeros_removed = 0
    hundreds_removed = 0
    nonchaning_removed = 0


    def all_list_in_interval(L):
        for i in L[1:]:
            if i < 1.05 * L[0] and i > L[0] * 0.95:
                continue
            else:
                return False
        return True

    for ind, j in enumerate(filtered):
        i = j[3:]
        if j[0] in not_include_templates:
            not_include += 1
        elif (mean_squared_error(i, [0] * len(i)) < 10):
            zeros_removed += 1
        elif mean_squared_error(i, [100] * len(i)) < 1:
            # print(j[0:3])
            hundreds_removed += 1
        elif all_list_in_interval(i):
            nonchaning_removed += 1
        else:
            data.append(j)
            data_color.append(filtered_colored[ind])

    print('there are : ' + str(len(data)) + " values left after deleting 100, and 0s")
    print("there were: " + str(not_include) + " vectors not included")
    print("there were: " + str(zeros_removed) + " vectors with zeros removed")
    print("there were: " + str(hundreds_removed) + " vectors with hundreds removed")
    print("there were: " + str(nonchaning_removed) + " vectors with non changing values removed")

    data_uncut = data

    ddd = group_signals_by_type_activity(data_uncut)
    subset_relations = [('RespondedExistence', 'Response'),
                        ('Response', 'AlternateResponse'),
                        ('AlternateResponse', 'ChainResponse'),
                        ('RespondedExistence', 'CoExistence'),
                        ('Response', 'Succession'),
                        ('AlternateResponse', 'AlternateSuccession'),
                        ('ChainResponse', 'ChainSuccession'),
                        ('CoExistence', 'Succession'),
                        ('Succession', 'AlternateSuccession'),
                        ('AlternateSuccession', 'ChainSuccession'),
                        ('RespondedExistence_r', 'Precedence'),
                        ('Precedence', 'AlternatePrecedence'),
                        ('AlternatePrecedence', 'ChainPrecedence'),
                        ('RespondedExistence_r', 'CoExistence'),
                        ('Precedence', 'Succession'),
                        ('AlternatePrecedence', 'AlternateSuccession'),
                        ('ChainPrecedence', 'ChainSuccession'),
                        ('AtLeast1', 'AtLeast2'),
                        ('AtLeast2', 'AtLeast3'),
                        ('AtMost3', 'AtMost2'),
                        ('AtMost2', 'AtMost1'),
                        ('AtMost1', 'Absence'),
                        ('AtLeast1', 'Init'),
                        ('AtLeast1', 'End')]
    subset_relations = [(a[1], a[0]) for a in subset_relations]

    templates = {'RespondedExistence',
                 'RespondedExistence_r',
                 'CoExistence',
                 'Precedence',
                 'AlternatePrecedence',
                 'ChainPrecedence',
                 'Response',
                 'AlternateResponse',
                 'ChainResponse',
                 'Succession',
                 'AlternateSuccession',
                 'ChainSuccession',
                 'NotRespondedExistence',
                 'NotCoExistence',
                 'NotPrecedence',
                 'NotChainPrecedence',
                 'NotResponse',
                 'NotChainResponse',
                 'NotSuccession',
                 'NotChainSuccession',
                 'Absence',
                 'AtLeast1',
                 'AtLeast2',
                 'AtLeast3',
                 'AtMost1',
                 'AtMost2',
                 'AtMost3',
                 'End',
                 'Init'}

    # Create a dictionary to store direct subset relationships
    graph = {letter: set() for letter in templates}

    # Populate the graph based on subset_relations

    for subset, superset in subset_relations:
        graph[subset].add(superset)

    pruned_list = []
    co_exist_list = set()
    for x in ddd.keys():
        for y in ddd[x].keys():
            maximal_letters = remove_subsets(ddd[x][y], subset_relations, graph)
            for dcl in maximal_letters - {'RespondedExistence_r'}:
                if dcl == 'CoExistence':
                    if (y[1].lstrip(), f' {y[0]}') not in co_exist_list:
                        pruned_list.append([dcl, y[0], y[1]] + list(x))
                        co_exist_list.add((y[0], y[1]))
                else:
                    pruned_list.append([dcl, y[0], y[1]] + list(x))
    return pruned_list, data_color

def remove_subsets(letters, subset_relations, graph):
    # Find all letters that are subsets of other letters
    subsets_to_remove = set()
    for letter in letters:
        # Find all supersets of the current letter
        supersets = find_all_supersets(graph, letter, set())
        subsets_to_remove.update(supersets)

    # The final set of letters excluding any that are subsets of another
    maximal_letters = letters - subsets_to_remove

    return maximal_letters

def clustering(pruned_list, linkage_method, linkage_metric, best_n_clusters):

    data_cut = []
    for data_point in pruned_list:
        data_cut.append(data_point[3:])

    '''build the clustering method'''
    Z = linkage(data_cut, method=linkage_method, metric=linkage_metric)  # metric='correlation'


    a = fcluster(Z, best_n_clusters, 'maxclust')
    clusters_dict = dict()
    for cluster_n, data in zip(a, pruned_list):
        if not cluster_n in clusters_dict:
            clusters_dict[cluster_n] = [data]
        else:
            clusters_dict[cluster_n].append(data)

    constraints = []
    cluster_bounds = []
    count = 0
    clusters_with_declare_names = {}

    order_cluster = [(i, len(clusters_dict[i])) for i in clusters_dict.keys()]
    # order_cluster = sorted(order_cluster, key=lambda x : -x[1])
    order_cluster = sorted(order_cluster, key=lambda x: -x[0])
    order_cluster = [key for key, _ in order_cluster]

    for key in order_cluster:
        # preprocess the clusters for plotting better (sorting them by similarity)
        clusters_dict[key] = sort_by_closest_neigbour_HEADER(clusters_dict[key])

        clusters_with_declare_names[key] = [clusters_dict[key][0][:3]]
        constraints.append(clusters_dict[key][0][3:])
        for i in clusters_dict[key][1:]:
            clusters_with_declare_names[key].append(i[:3])
            constraints.append(i[3:])
        count += len(clusters_dict[key])
        cluster_bounds.append(count)
        print("Cluster size: " + str(len(clusters_dict[key])))

    print('number of clusters: ' + str(len(clusters_dict)))
    print("clusters were ordered in the : ")
    print(order_cluster)

    return order_cluster, clusters_with_declare_names, cluster_bounds, clusters_dict, constraints
def constraints_export(clusters_with_declare_names, peaks, w,clusters_dict):
    file_path = r"output_files/"

    import statistics
    def generate_natural_language(constraint):
        c, a, b = constraint
        if c == "RespondedExistence":
            description = f"If {a} occurs, {b} occurs as well."
        elif c == "CoExistence":
            description = f"{a} and {b} always occur together."
        elif c == "Response":
            description = f" If {a} occurs, then {b} occurs after it."
        elif c == "AlternateResponse":
            description = f"If {a} occurs, then {b} occurs afterwards before {a} recurs."
        elif c == "ChainResponse":
            description = f"If {a} occurs, then {b} occurs immediately after it."
        elif c == "Precedence":
            description = f"{b} occurs only if preceded by {a}."
        elif c == "AlternatePrecedence":
            description = f"{b} occurs only if preceded by {a} with no other {b} in between."
        elif c == "ChainPrecedence":
            description = f"{b} occurs only if {a} occurs immediately before it. "
        elif c == "Succession":
            description = f"{a} occurs if and only if it is followed by {b}."
        elif c == "AlternateSuccession":
            description = f"{a} and {b} occur if and only if they follow one another, alternating."
        elif c == "ChainSuccession":
            description = f"{a} and {b} occurs if and only if {b} immediately follows {a}."
        elif c == "Init":
            description = f"{a} is the first to occur."
        elif c == "End":
            description = f"{a} is the last to occur."
        elif c == "Absence":
            description = f"{a} must never occur."
        elif c == "AtMost1":
            description = f"{a} occurs at most once."
        elif c == "AtMost2":
            description = f"{a} occurs at most two times."
        elif c == "AtMost3":
            description = f"{a} occurs at most three times."
        elif c == "AtLeast1":
            description = f"{a} occurs at least once."
        elif c == "AtLeast2":
            description = f"{a} occurs at least two times."
        elif c == "AtLeast3":
            description = f"{a} occurs at least three times."
        else:
            description = f"The constraint is un known!"

        return description

    def stat_extract(value, peaks, w):
        new_list = []
        for const in value:
            stat_dic = {'constraint type': const[0], 'first parameter': const[1], 'second parameter': const[2],
                        'description': generate_natural_language(const[0:3])}
            c = const[3:]
            seg_num = 1
            for p in range(len(peaks) + 1):
                if seg_num == 1:
                    stat_dic[f'segment_{seg_num}'] = round(statistics.mean(c[0:(peaks[p] - (w - 1) + 1)]), 2)
                    stat_dic[f'~segment_{seg_num}'] = round(statistics.mean(c[(peaks[p] - (w - 1)):]), 2)
                elif seg_num <= len(peaks):
                    stat_dic[f'segment_{seg_num}'] = round(
                        statistics.mean(c[(peaks[p - 1] - (w - 1)):(peaks[p] - (w - 1) + 1)]), 2)
                    stat_dic[f'~segment_{seg_num}'] = round(
                        statistics.mean(c[0:(peaks[p - 1] - (w - 1) + 1)] + c[(peaks[p] - (w - 1)):]), 2)
                else:
                    stat_dic[f'segment_{seg_num}'] = round(statistics.mean(c[(peaks[p - 1] - (w - 1)):]), 2)
                    stat_dic[f'~segment_{seg_num}'] = round(statistics.mean(c[0:(peaks[p - 1] - (w - 1) + 1)]), 2)
                stat_dic[f'delta_segment_{seg_num}'] = round(
                    stat_dic[f'segment_{seg_num}'] - stat_dic[f'~segment_{seg_num}'], 2)
                seg_num += 1
            new_list.append(stat_dic)
        return new_list


    data_str_keys = {str(key): stat_extract(value, peaks, w) for key, value in clusters_dict.items()}
    # Write the dictionary to a JSON file
    with open(file_path+"data1.json", 'w') as file:
        json.dump(data_str_keys, file, indent=4)

    for cl in clusters_with_declare_names.keys():
        export_constraints_per_cluster(clusters_with_declare_names[cl], file_path + f'constraints_{cl}.json')
        # prune_constraints_minerful(file_path + f'constraints_{cl}.json', file_path + f'constraints_{cl}_pruned.csv')




    print(f"Dictionary successfully saved to {file_path}")

def PELT_change_points(order_cluster,clusters_dict):
    def dingOptimalNumberOfPoints(algo):
        point_detection_penalty = 50
        x_lines = algo.predict(pen=point_detection_penalty)

        while point_detection_penalty >= len(x_lines):
            point_detection_penalty -= 1
            x_lines = algo.predict(pen=point_detection_penalty)

        if len(x_lines) > 15:
            x_lines = x_lines[-1:]
        return x_lines

    horisontal_separation_bounds_by_cluster = {}
    # in ths case we want to detect the drifts in the whole range of constrains at the same time

    dd = []
    for dk in order_cluster:
        for i in clusters_dict[dk]:
            dd.append(i[3:])

    sig = np.array(dd)
    signal = np.transpose(sig)
    # algo = rpt.Pelt(model="rbf", custom_cost=c).fit(signal)
    algo = rpt.Pelt(model="rbf").fit(signal)
    x_lines = dingOptimalNumberOfPoints(algo)
    horisontal_separation_bounds_by_cluster[0] = x_lines
    # pen - penalizing the number of change points

    return x_lines


import plotly.graph_objects as go
import numpy as np


def plot_figures(df, masks, n_bin, map_range, peaks, constraints, w, cluster_bounds, clusters_with_declare_names,
                 data_color, corr_mat, WINDOWS):
    # every = 1
    color_theme_drift_map = 'Blues'


    ################## Figure 3: Control Flow Features and Change Points ####################
    L1 = []
    for k in clusters_with_declare_names.keys():
        for s in clusters_with_declare_names[k]:
            L1.append(s)
    L1_index = {tuple(sublist[0:3]): idx for idx, sublist in enumerate(L1)}

    # Sort and filter L2 based on L1
    L2_ordered = sorted(
        [sublist for sublist in data_color if tuple(sublist[0:3]) in L1_index],
        key=lambda x: L1_index[tuple(x[0:3])]
    )
    L2_ordered = [x[3:] for x in L2_ordered]

    # Prepare data for the heatmap
    data_c = []
    for i in range(len(constraints)):
        new_list = [np.nan] * (w - 1) + constraints[i] + [np.nan] * (w - 1)
        data_c.append(new_list)

    # Create the first subplot (Heatmap with clusters and peaks)
    fig3 = go.Figure()
    ordered_keys = sorted(L1_index, key=L1_index.get)
    ordered_keys =[str(x) for x in ordered_keys]
    # Add heatmap
    fig3.add_trace(go.Heatmap(
        z=data_c,
        y = ordered_keys,
        colorscale=color_theme_drift_map,
        zmin=0,
        zmax=100,
        colorbar=dict(
            title=dict(
                text = "Confidence",
            font=dict(size=18)),
            tickfont=dict(size=16)
        )
    ))

    # Add cluster boundary lines
    for bound in cluster_bounds:
        fig3.add_shape(
            type="line",
            x0=0,
            x1=n_bin-2,
            y0=bound-0.5,
            y1=bound-0.5,
            line=dict(color="black", dash="dash", width=2)
        )

    # Add vertical peak lines
    for pp in peaks:
        fig3.add_shape(
            type="line",
            x0=pp,
            x1=pp,
            y0=0,
            y1=len(data_c),
            line=dict(color="red", dash="dash", width=4)
        )

    # Set x-axis ticks and labels
    x_tick = []
    x_labels = []
    for i in range(len(peaks) + 1):
        if i == 0:
            x_tick.append(round(peaks[i] / 2, 0))
        elif i == len(peaks):
            x_tick.append(round((n_bin + peaks[i - 1]) / 2, 0))
        else:
            x_tick.append(round((peaks[i] + peaks[i - 1]) / 2, 0))
        x_labels.append(f'segment {i + 1}')

    fig3.update_layout(
        plot_bgcolor="gray",  # Gray plot background
        title="Control Flow Features and Change Points Overview",
        xaxis=dict(
            title="Segments",
            tickvals=x_tick,
            ticktext=x_labels,
            tickfont=dict(size=16)
        ),
        yaxis=dict(
            title="Clusters",
            tickvals=[round((cb + cluster_bounds[i - 1]) / 2, 0) if i > 0 else round(cb / 2, 0) for i, cb in
                      enumerate(cluster_bounds)],
            ticktext=[f'cluster {len(cluster_bounds) - i}' for i in range(len(cluster_bounds))],
            tickfont=dict(size=16)
        ),
        template="plotly_white",
        height=500,
        width=750
    )

    ##################### Figure 4: Correlation Heatmap #####################
    corr_mat_transposed = np.array(corr_mat).T  # Transpose the matrix

    # Create correlation heatmap
    fig4 = go.Figure()

    fig4.add_trace(go.Heatmap(
        z=corr_mat_transposed,
        colorscale=[[0.0, '#8B0000'], [0.5, 'white'], [1.0, '#00008B']],
        zmin=-1,
        zmax=1,
        colorbar=dict(
            title=dict(
                text = "Correlation",
            font=dict(size=18)),
            tickfont=dict(size=16)
        )
    ))

    fig4.update_layout(
        title="Correlation Between Features and Segments",
        xaxis=dict(
            title="Segments",
            tickvals=np.arange(len(x_labels)),
            ticktext=x_labels,
            tickfont=dict(size=16)
        ),
        yaxis=dict(
            title="Features",
            tickvals=np.arange(len(cluster_bounds)),
            ticktext=[f'cluster {i+1}' for i in range(len(cluster_bounds))],
            tickfont=dict(size=16)
        ),
        template="plotly_white",
        height=600,
        width=600
    )

    return fig3, fig4


# We made the visualization in Plotly later, to allow for intractive visualization
# def plot_figures(df, masks, n_bin, map_range, peaks, constraints, w, cluster_bounds,clusters_with_declare_names, data_color, corr_mat,WINDOWS):
#     every = 2
#     color_theme_drift_map = 'Blues'
#
#     def insert_and_clean_np(array, new_number):
#         """
#         Insert a new number into a sorted NumPy array, maintain the order,
#         and remove any numbers within a distance of 1 from the new number.
#
#         Parameters:
#         - array (np.ndarray): A sorted NumPy array of numbers.
#         - new_number (int or float): The number to insert.
#
#         Returns:
#         - np.ndarray: The updated NumPy array.
#         """
#         if new_number not in array:
#             # Find the correct insertion point
#             insertion_index = np.searchsorted(array, new_number)
#
#             # Insert the new number into the array
#             array = np.insert(array, insertion_index, new_number)
#
#             # Identify elements to keep (distance > 1 from the new number)
#             mask = np.abs(array - new_number) != 1
#
#             # Return the updated array
#             return array[mask]
#             # return mask
#         else:
#             return array
#
#     ################## Figure 3############################################
#     L1 = []
#     for k in clusters_with_declare_names.keys():
#         for s in clusters_with_declare_names[k]:
#             L1.append(s)
#     L1_index = {tuple(sublist[0:3]): idx for idx, sublist in enumerate(L1)}
#
#     # Sort and filter L2 based on L1
#     L2_ordered = sorted([sublist for sublist in data_color if tuple(sublist[0:3]) in L1_index],
#                         key=lambda x: L1_index[tuple(x[0:3])])
#     L2_ordered = [x[3:] for x in L2_ordered]
#     data_c = []
#     data_c_color_1 = []
#     data_c_color_2 = []
#     for i in range(len(constraints)):
#         new_list = [0] * (w - 1) + constraints[i] + [0] * (w - 1)
#         new_list_color_1 = [True] * (w - 1) + L2_ordered[i] + [True] * (w - 1)
#         new_list_color_2 = [True] * (w - 1) + [False] * len(L2_ordered[i]) + [True] * (w - 1)
#         data_c.append(new_list)
#         data_c_color_1.append(new_list_color_1)
#         data_c_color_2.append(new_list_color_2)
#
#     # Create a new figure with two subplots with different heights
#     fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 9),
#                                    gridspec_kw={'height_ratios': [10, 1]})  # First plot 3 times bigger
#
#     fig3.suptitle('Control Flow Features and Change Points Overview', fontsize=20)
#
#     # sns.heatmap(data_c, mask=list(mask_ax1[windows.index(w)]), linewidth=0, cmap=color_theme_drift_map,  vmin=0, vmax=100, ax=ax1)
#     light_gray_cmap = ListedColormap(['#d3d3d3'])
#     sns.heatmap(data_c_color_1, mask=np.array(data_c_color_2), cmap=light_gray_cmap, cbar=False, ax=ax1)
#     ax1.set_facecolor("gray")
#     sns.heatmap(np.array(data_c), mask=np.array(data_c_color_1), linewidth=0, cmap=color_theme_drift_map, cbar=True,
#                 vmin=0, vmax=100, ax=ax1)
#
#     cbar1 = ax1.collections[1].colorbar
#     cbar1.ax.tick_params(labelsize=16)
#     cbar1.set_label('Confidence', fontsize=18)
#
#     original_ticks = cbar1.get_ticks()  # Get the current ticks (e.g., [0, 20, 40, ..., 100])
#     normalized_ticks = np.linspace(0, 1, len(original_ticks))  # Map the ticks to 0–1
#
#     # Update the color bar with normalized ticks
#     cbar1.set_ticks(original_ticks)  # Retain original tick positions
#     cbar1.set_ticklabels([f"{tick:.1f}" for tick in normalized_ticks])  # Set new labels (0–1)
#
#     for bound in cluster_bounds:
#         ax1.axhline(y=bound, color='black', linestyle='--', linewidth=2)
#
#     for pp in peaks:
#         ax1.axvline(x=0.5 + pp - 0.05, color='red', linestyle='--', linewidth=4)
#
#
#     ax1.set_xticks([])
#     y_tick = []
#     y_labels = []
#     for i, cb in enumerate(cluster_bounds):
#         if i == 0:
#             y_tick.append(round(cb / 2, 0))
#         else:
#             y_tick.append(round((cluster_bounds[i] + cluster_bounds[i - 1]) / 2, 0))
#         y_labels.append(f'cluster {len(cluster_bounds) - i}')
#     ax1.set_yticks(y_tick)
#     ax1.set_yticklabels(y_labels)
#
#     x_tick = []
#     x_labels = []
#     for i in range(len(peaks) + 1):
#         if i == 0:
#             x_tick.append(round(peaks[i] / 2, 0))
#         elif i == len(peaks):
#             x_tick.append(round((n_bin + peaks[i - 1]) / 2, 0))
#         else:
#             x_tick.append(round((peaks[i] + peaks[i - 1]) / 2, 0))
#         x_labels.append(f'segment {i + 1}')
#     ax1.set_xticks(x_tick)
#     ax1.set_xticklabels(x_labels)
#     ax1.tick_params(axis='x', labelsize=16, rotation=0)
#     ax1.tick_params(axis='y', labelsize=16, rotation=0)
#     ax1.grid(False)
#
#
#     sns.heatmap(df.loc[w].to_frame().T, mask=masks[WINDOWS.index(w)], cmap="Reds", linewidth=0,
#                 ax=ax2)
#
#     ticks = np.arange(0, n_bin - 1, 2)
#     pk_id = []
#     for new_number in peaks:
#         ticks = insert_and_clean_np(ticks, new_number)
#         pk_id.append(np.where(ticks == new_number)[0][0])
#     ticks_labels = [str(round(x * (100 / n_bin))) + "% (" + str(round(map_range[str(x)], 1)) + ")" for x in (ticks + 1)]
#
#     print(pk_id)
#     ax2.set_xticks(0.5 + ticks)
#     ax2.set_xticklabels(ticks_labels)
#
#     tick_labels = ax2.get_xticklabels()  # Get all x-axis tick labels
#     print(tick_labels)
#     for pk in pk_id:
#         # Change the color of specific ticks
#         tick_labels[pk].set_color('red')  # Change tick at position 2 to red
#
#     ax2.set_facecolor("gray")
#     ax2.grid(False)
#     # ax2.set_title('sliding window analysis', fontsize=20)
#     ax2.set_xlabel('traces', fontsize=18)
#     ax2.set_ylabel('window size', fontsize=18)
#     ax2.tick_params(axis='x', labelsize=16, rotation=90)
#     ax2.tick_params(axis='y', labelsize=16)
#
#     # Adjust colorbar font size
#     cbar2 = ax2.collections[0].colorbar
#     cbar2.ax.tick_params(labelsize=16)
#     cbar2.set_label('ldist', fontsize=18)
#
#     ##################### Figure 4####################
#     from matplotlib.colors import LinearSegmentedColormap
#
#     corr_mat_transposed = np.array(corr_mat).T  # Transpose the matrix
#
#     # Output the list of correlations
#     fig4 = plt.figure(figsize=(8, 7))
#     fig4.suptitle('Correlation Between Features and Segments', fontsize=20)
#
#     colors = ['#8B0000', 'white', '#00008B']  # Dark blue, white, dark red
#     # Create a colormap
#     custom_cmap = LinearSegmentedColormap.from_list("custom_cmp", colors)
#
#     # Create a heatmap with black grid lines
#     ax3 = sns.heatmap(
#         corr_mat_transposed,
#         annot=False,
#         cmap=custom_cmap,
#         cbar=True,
#         vmin=-1,
#         vmax=1,
#         linewidths=0.5,
#         linecolor='gray'
#     )
#
#     ax3.set_xticks(0.5 + np.arange(0, len(peaks) + 1))
#     ax3.set_xticklabels(x_labels)
#
#     ax3.set_yticks(0.5 + np.arange(0, len(cluster_bounds)))
#     ax3.set_yticklabels(y_labels[::-1])
#
#     # Rotate the y-axis tick labels instead of x-axis
#     ax3.tick_params(axis='y', rotation=0, labelsize=16)
#     ax3.tick_params(axis='x', labelsize=16)
#
#     ax3.invert_yaxis()
#
#     cbar3 = ax3.collections[0].colorbar
#     cbar3.ax.tick_params(labelsize=16)
#     cbar3.set_label('correlation', fontsize=18)
#
#     buf = BytesIO()
#     fig3.savefig(buf, format="png", bbox_inches='tight')
#     fig_data3 = base64.b64encode(buf.getbuffer()).decode("ascii")
#
#     buf = BytesIO()
#     fig4.savefig(buf, format="png", bbox_inches='tight')
#     fig_data4 = base64.b64encode(buf.getbuffer()).decode("ascii")
#
#     return f'data:image/png;base64,{fig_data3}', f'data:image/png;base64,{fig_data4}'



def report(data, cluster, segment):
    if f'segment_{segment}' not in data[str(cluster)][0].keys():
        print('segment does not exist')
    else:
        file_path = r"output_files/"
        prune_constraints_minerful(file_path + f'constraints_{cluster}.json', file_path + f'constraints_{cluster}_pruned.csv')
        fp = f"output_files/constraints_{cluster}_pruned.csv"
        df = pd.read_csv(fp, sep=";")
        df = df.applymap(lambda x: x.strip("'\"") if isinstance(x, str) else x)
        df.columns = df.columns.str.strip("'\"")
        df.set_index('Constraint', inplace=True)
        list = []
        for d in data[str(cluster)]:
            if d['second parameter'] == '':
                const = f"{d['constraint type']}({d['first parameter']})"
                if df.loc[f'{const}', 'Redudant'] == False:
                    list.append((d['description'], round(d[f'segment_{segment}'] - d[f'~segment_{segment}'], 2)))

            else:
                const = f"{d['constraint type']}({d['first parameter']},{d['second parameter']})"
                if df.loc[f'{const}', 'Redudant'] == False:
                    list.append((d['description'], round(d[f'segment_{segment}'] - d[f'~segment_{segment}'], 2)))

        list_sorted = sorted(list, key=lambda x: x[1], reverse=False)[0:]
        list_sorted_reverse = sorted(list, key=lambda x: x[1], reverse=True)[0:]

        return [f'rank {k+1}: {x[0]}, with score {x[1]}' for k,x in enumerate(list_sorted)], [f'rank {k+1}: {x[0]}, with score {x[1]}' for k,x in enumerate(list_sorted_reverse)]
def decl2NL(cluster, segment):
    file_path = r"output_files/data1.json"
    with open(file_path, 'r') as file:
        data = json.load(file)
    list_sorted, list_sorted_reverse = report(data, cluster, segment)
    return list_sorted, list_sorted_reverse



def apply_feature_extraction(theta_cvg):
    """Main function to apply the analysis."""
    pruned_list, data_color = prune_signals(theta_cvg)
    redis_client.set("pruned_list",json.dumps(pruned_list))
    redis_client.set("data_color",json.dumps(data_color))
    return "Feature Generation Done!"

def apply_X(n_bin, w, n_clusters):
    """Main function to apply the analysis."""
    # if w not in WINDOWS:
    #     WINDOWS.append(w)
    #     WINDOWS.sort(reverse=True)
    WINDOWS = [w]

    pruned_list = json.loads(redis_client.get("pruned_list"))
    data_color = json.loads(redis_client.get("data_color"))

    log_command("clustering started!")
    order_cluster, clusters_with_declare_names, cluster_bounds, clusters_dict, constraints = clustering(pruned_list,
                                                                                                        linkage_method,
                                                                                                        linkage_metric,
                                                                                                        n_clusters)
    log_command("clustering done!")

    df = pd.read_json(redis_client.get("df"), orient="split")
    masks = json.loads(redis_client.get("masks"))
    map_range = json.loads(redis_client.get("map_range"))
    peaks = json.loads(redis_client.get("peaks"))


    redis_client.set("segments_count",len(peaks) + 1)
    redis_client.set("clusters_count", len(clusters_with_declare_names.keys()))

    # PELT_change_points(order_cluster, clusters_dict)
    log_command("exporting the constraints started!")
    constraints_export(clusters_with_declare_names, peaks, w, clusters_dict)
    log_command("exporting the constraints done!")

    log_command("correlation matrix is being calculated!")
    corr_mat = correlation_calc(peaks, w, constraints, clusters_dict)
    log_command("correlation matrix calculated!")

    log_command("figure 3 and 4 are being generated!")
    fig3_path, fig4_path = plot_figures(df, masks, n_bin, map_range, peaks,
                                        constraints, w, cluster_bounds,
                                        clusters_with_declare_names, data_color, corr_mat,
                                        WINDOWS)
    log_command("figure 3 and 4 are generated!")
    return fig3_path, fig4_path