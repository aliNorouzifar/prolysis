import pandas as pd
import matplotlib
matplotlib.use('agg')
import scipy.signal as sci_sig
import os
import itertools
import pm4py
import json
import ast
from prolysis.analysis.python_emsc import stochastic
import plotly.graph_objects as go
import numpy as np
from prolysis.util.redis_connection import redis_client
from prolysis.util.logging import log_command

# # Constants
OUTPUT_DIR = "output_files"

os.makedirs(OUTPUT_DIR, exist_ok=True)  # Ensure output directory exists
color_theme_drift_map = 'Blues'


def emd_dist(bin1, bin2):
    if isinstance(bin1, list):
        bin1 = pd.Series(bin1)
    if isinstance(bin2, list):
        bin2 = pd.Series(bin2)
    lang1 = bin1.value_counts(normalize=True)
    lang1 = {ast.literal_eval(key): value for key, value in lang1.items()}
    lang2 = bin2.value_counts(normalize=True)
    lang2 = {ast.literal_eval(key): value for key, value in lang2.items()}
    dist = round(stochastic.compare_languages_levenshtein(lang1, lang2),2)
    return dist

def bins_generation(kpi, n_bin):
    """Generate bins and map ranges for a given KPI."""
    case_table = pd.read_csv("output_files/out.csv").sort_values(by=[kpi])

    map_range = {i: float(case_table[kpi].iloc[round((i / n_bin) * len(case_table[kpi]))]) for i in range(n_bin)}
    map_range[n_bin] = float(case_table[kpi].iloc[-1])

    bin_size = round(len(case_table) / n_bin)
    bins = [
        ((min(case_table[point:point + bin_size][kpi]), max(case_table[point:point + bin_size][kpi])), idx,
         case_table[point:point + bin_size]['trace'])
        if point + bin_size < len(case_table) else
        ((min(case_table[point:][kpi]), max(case_table[point:][kpi])), idx, case_table[point:]['trace'])
        for idx, point in enumerate(range(0, len(case_table), bin_size))
    ]

    return bins, map_range, case_table


def sliding_window(bins, n_bin,WINDOWS):
    """Perform sliding window analysis for change detection."""
    df = pd.DataFrame(0.0, index=WINDOWS, columns=[i for i in range(1, n_bin)])
    for window_size in WINDOWS:
        for mid in range(window_size, n_bin - window_size + 1):
            left = [item for b in bins[mid - window_size:mid] for item in b[2]]
            # left = combine_dics([b[2] for b in bins[mid - window_size:mid]])
            right = [item for b in bins[mid:mid + window_size] for item in b[2]]
            # right = combine_dics([b[2] for b in bins[mid:mid + window_size]])
            df.at[window_size, mid] = emd_dist(left,right)
    masks = [
        [True] * (window - 1) + [False] * (n_bin - 2 * window + 1) + [True] * (window - 1)
        for window in WINDOWS
    ]
    return df, masks


def segmentation(df,bins,n_bin,w,sig):
    peaks, _ = sci_sig.find_peaks(df.loc[w], height=[sig])
    peaks = [int(p) for p in peaks]
    segments = []
    segments_ids = []
    last_p = -1
    x_state = 0
    for p in peaks:
        new = (x_state, x_state, pd.Series(list(itertools.chain.from_iterable([x[2] for x in bins[last_p + 1:p + 1]]))))
        # new_ids = [item for b in bins[last_p + 1:p + 1] for item in b[2].index]
        new_ids = [item for b in bins[last_p + 1:p + 1] for item in b[2].index]
        segments.append(new)
        segments_ids.append(new_ids)
        last_p = p
        x_state += 1

    new = (x_state, x_state, pd.Series(list(itertools.chain.from_iterable([x[2] for x in bins[last_p + 1:]]))))
    new_ids = [item for b in bins[last_p + 1:] for item in b[2].index]
    segments.append(new)
    segments_ids.append(new_ids)

    state = n_bin
    m_dend = []
    cal_list = {}

    itr = 1
    data_points = []

    matrices = []
    labels = []
    mins_vec = []

    state_dic = {}

    while len(segments) > 1:
        dist_matrix = np.ones((len(segments), len(segments)))
        for i in range(0, len(segments)):
            for j in range(i + 1, len(segments)):
                if (segments[i][0], segments[j][0]) not in cal_list:
                    cal_list[(segments[i][0], segments[j][0])] = emd_dist(segments[i][2], segments[j][2])
                dist_matrix[i, j] = cal_list[(segments[i][0], segments[j][0])]
                data_points.append((itr, dist_matrix[i, j]))
        matrices.append(dist_matrix)
        labels.append([b[1] for b in segments])
        min_dist_ind = np.unravel_index(np.argmin(dist_matrix, axis=None), dist_matrix.shape)
        mins_vec.append(dist_matrix[min_dist_ind])
        m_dend.append([segments[min_dist_ind[0]][1], segments[min_dist_ind[1]][1], dist_matrix[min_dist_ind],
                       len(segments[min_dist_ind[0]][2]) + len(segments[min_dist_ind[1]][2])])
        state_dic[state] = ((segments[min_dist_ind[0]][1], len(segments[min_dist_ind[0]][2])),
                            (segments[min_dist_ind[1]][1], len(segments[min_dist_ind[1]][2])))
        segments = [segments[k] for k in range(0, len(segments)) if (k != min_dist_ind[0] and k != min_dist_ind[1])] + [
            ((segments[min_dist_ind[0]][0], segments[min_dist_ind[1]][0]), state,
             pd.concat([segments[min_dist_ind[0]][2], segments[min_dist_ind[1]][2]]))]
        state = state + 1
        itr += 1


    ittr = 0

    new_m = matrices[ittr]
    for i in range(0, len(matrices[ittr])):
        for j in range(i, len(matrices[ittr][i])):
            if i == j:
                new_m[j, i] = 0
            else:
                new_m[j, i] = new_m[i, j]
    return segments, segments_ids, new_m, peaks


def plot_figures_EMD(df, masks, n_bin, map_range,WINDOWS):

    every = 1

    # Replace masked values with NaN for transparency
    heatmap_data = np.ma.array(df, mask=masks).filled(np.nan)

    # Generate x-tick labels
    x_labels = [
        f"{round(x * (100 / n_bin))}% ({round(map_range[x], 1)})"
        for x in range(1, n_bin, every)
    ]

    # Generate y-tick labels
    y_labels = WINDOWS

    # Define the colorscale (exclude gray, focus on valid data range)
    custom_colorscale = "Reds"  # Default Plotly Reds colorscale

    # Create the Plotly Heatmap
    fig1 = go.Figure(
        data=go.Heatmap(
            x=x_labels,
            z=heatmap_data,  # Heatmap values
            y=y_labels,
            colorscale=custom_colorscale,  # Valid data colorscale
            zmin=np.nanmin(df),  # Min value for color scaling
            zmax=np.nanmax(df),  # Max value for color scaling
            colorbar=dict(title="ldist"),  # Colorbar title
            showscale=True,  # Show the scale
        )
    )


    fig1.update_layout(
        # paper_bgcolor="gray",  # Set the figure's background color to gray
        plot_bgcolor="gray",  # Set the plot area background color to gray
        title=dict(
            text="Sliding Window Analysis",  # Title text
            x=0.5,  # Horizontal alignment (0=left, 0.5=center, 1=right)
            xanchor="center",  # Ensures proper anchoring
            font=dict(size=18),  # Optional: Set title font size
        ),
        xaxis=dict(
            title="Traces",
        ),
        yaxis=dict(
            title="Window Size",
        ),
        template="plotly_white",
    )
    fig1.update_xaxes(tickangle=90)
    return fig1

def plot_figures_segments(dist_matrix, peaks):
    labels = [f"segment{i}" for i in range(1, dist_matrix.shape[0] + 1)]
    # Create the Plotly Heatmap
    fig2 = go.Figure(
        data=go.Heatmap(
            z=dist_matrix,  # Heatmap values
            x=labels,  # X-axis labels
            y=labels,  # Y-axis labels
            colorscale="Reds",  # Color scheme
            colorbar=dict(
                title=dict(
                    text = "ldist",  # Colorbar title
                font=dict(size=18)),  # Font size for the colorbar title
                tickfont=dict(size=18)  # Font size for the colorbar ticks
            ),
            zmin=np.min(dist_matrix),  # Minimum value for color scaling
            zmax=np.max(dist_matrix),  # Maximum value for color scaling
        )
    )

    # Update layout for the figure
    fig2.update_layout(
        title=dict(
            text="Segments Comparison",
            font=dict(size=20),  # Title font size
            x=0.5  # Center the title
        ),
        xaxis=dict(
            title="",
            tickfont=dict(size=18)  # Font size for x-axis ticks
        ),
        yaxis=dict(
            title="",
            tickfont=dict(size=18)  # Font size for y-axis ticks
        ),
        template="plotly_white",
        width = 580,  # Set the width of the figure
        height = 500  # Set the height of the figure
    )

    return fig2


### previous version of the figures with mathplotlib
# def plot_figures(df, masks, n_bin, map_range, dist_matrix, peaks, w,WINDOWS):
#     every = 2
#     """Generate heatmaps and comparison plots."""
#     # Sliding Window Heatmap
#     fig1, ax1 = plt.subplots(figsize=(15, 3))
#     sns.heatmap(df, cmap="Reds", mask=np.array(masks), ax=ax1)
#     ax1.set_xticks(0.5 + np.arange(0, n_bin - 1, 3))
#     ax1.set_xticklabels(
#         [f"{round(x * (100 / n_bin))}% ({round(map_range[x], 1)})" for x in range(1, n_bin, 3)]
#     )
#     ax1.set_facecolor("gray")
#     ax1.set_title("Sliding Window Analysis")
#     ax1.set_xlabel("Traces")
#     ax1.set_ylabel("Window Size")
#     ax1.set_xticks(0.5 + np.arange(0, n_bin - 1, every))
#     ax1.set_xticklabels(
#         [str(round(x * (100 / n_bin))) + "% (" + str(round(map_range[x], 1)) + ")" for x in np.arange(1, n_bin, every)])
#     ax1.set_yticks([x+0.5 for x in range(0,len(WINDOWS))], labels=WINDOWS)
#     plt.xticks(rotation=90)
#     plt.close(fig1)
#
#     cmap = plt.cm.Reds
#     fig2 = plt.figure(figsize=(7, 7))
#     ax = sns.heatmap(dist_matrix, cmap=cmap, xticklabels=['segment' + str(i) for i in range(1, dist_matrix.shape[0] + 1)],
#                      yticklabels=['segment' + str(i) for i in range(1, dist_matrix.shape[0] + 1)])
#
#     fig2.suptitle('segments comparison', fontsize=20)
#     plt.xticks(fontsize=18)
#     plt.yticks(fontsize=18)
#     cbar = ax.collections[0].colorbar
#     cbar.ax.tick_params(labelsize=18)
#     plt.xlabel(' ', fontsize=18)
#     cbar.set_label('ldist', fontsize=18)
#     plt.close(fig2)
#
#     buf = BytesIO()
#     fig1.savefig(buf, format="png", bbox_inches = 'tight')
#     # Embed the result in the html output.
#     fig_data1 = base64.b64encode(buf.getbuffer()).decode("ascii")
#
#     buf = BytesIO()
#     fig2.savefig(buf, format="png", bbox_inches='tight')
#     # Embed the result in the html output.
#     fig_data2 = base64.b64encode(buf.getbuffer()).decode("ascii")
#     return f'data:image/png;base64,{fig_data1}', f'data:image/png;base64,{fig_data2}'
#
#
def export_logs(segments_ids):
    """Export logs for each segment."""
    case_table = pd.read_csv("output_files/out.csv")
    event_file = "output_files/out_event.csv"
    event_table = pd.read_csv(event_file)

    os.makedirs("event_logs/exported_logs", exist_ok=True)
    for idx, segment_id_set in enumerate(segments_ids, start=1):
        segment_cases = case_table.loc[segment_id_set, 'case_id']
        segment_log = event_table[event_table['case_id'].isin(segment_cases)]
        segment_log = pm4py.format_dataframe(segment_log, case_id="case_id", activity_key="activity_name",
                                             timestamp_key="timestamp")
        event_log = pm4py.convert_to_event_log(segment_log)
        pm4py.write_xes(event_log, f"event_logs/exported_logs/segment_{idx}.xes")



def apply_EMD(n_bin, w, kpi):
    """Main function to apply the analysis."""
    # if w not in WINDOWS:
    #     WINDOWS.append(w)
    #     WINDOWS.sort(reverse=True)
    WINDOWS = [w]

    log_command("binning started!")
    bins, map_range, case_table = bins_generation(kpi, n_bin)
    log_command("binning done!")

    log_command("sliding window started!")
    df, masks = sliding_window(bins, n_bin,WINDOWS)
    log_command("sliding window done!")

    redis_client.set("df", df.to_json(orient="split"))
    redis_client.set("masks", json.dumps(masks))
    redis_client.set("map_range", json.dumps(map_range))
    redis_client.set("bins",json.dumps([(b[0],b[1],list(b[2].items())) for b in bins]))
    redis_client.set("max_dist", df.iloc[0].max())

    # fig1_path, fig2_path = plot_figures(df, masks, n_bin, map_range, dist_matrix, peaks, w,WINDOWS)
    log_command("figure 1 is being generated!")
    fig1 = plot_figures_EMD(df, masks, n_bin, map_range,WINDOWS)
    log_command("figure 1 generated!")
    return fig1

def apply_segmentation(n_bin, w, signal_threshold):
    """Main function to apply the analysis."""
    df = pd.read_json(redis_client.get("df"), orient="split")


    bins = [(b[0], b[1], pd.Series(dict(b[2]))) for b in json.loads(redis_client.get("bins"))]
    log_command("segmentation started!")
    segments, segments_ids, dist_matrix, peaks = segmentation(df, bins, n_bin, w, signal_threshold)
    log_command("segmentation done!")


    redis_client.set("peaks",json.dumps(peaks))
    redis_client.set("segments_ids", json.dumps(segments_ids))



    # fig1_path, fig2_path = plot_figures(df, masks, n_bin, map_range, dist_matrix, peaks, w,WINDOWS)
    log_command("figure 2 is being generated!")
    fig2 = plot_figures_segments(dist_matrix, peaks)
    log_command("figure 2 generated!")


    map_range = json.loads(redis_client.get("map_range"))

    peak_explanations = []
    if peaks:
        for i, p in enumerate(peaks):
            if i == 0:
                peak_explanations.append(f"Segment {i + 1}: From the beginning to {round(map_range[str(peaks[0]+1)],2)}")
            else:
                peak_explanations.append(f"Segment {i + 1}: From {round(map_range[str(peaks[i - 1]+1)],2)} to {round(map_range[str(peaks[i]+1)],2)}")

        peak_explanations.append(f"Segment {len(peaks)+1}: From {round(map_range[str(peaks[-1]+1)],2)} to the end")
    else:
        peak_explanations.append("We have only one Segment!")

    return fig2,peak_explanations
