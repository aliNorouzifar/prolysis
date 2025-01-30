import os
import subprocess

def mine_minerful_for_declare_constraints(window_size,sliding_window_size):
    # input_log_path = os.getcwd() + r"/output_files/log_ordered.xes"
    # output_log_path = os.getcwd()+ r"/output_files/behavioral_signals.csv"
    input_log_path = r"output_files/log_ordered.xes"
    output_log_path = r"output_files/behavioral_signals.csv"
    env = dict(os.environ)
    env = dict(os.environ)
    env['JAVA_OPTS'] = 'foo'
    subprocess.call(['java', '-version'])
    file_input = input_log_path
    subprocess.call([
        'java', "-Xmx16G",
        '-cp', f'MINERful.jar',
        'minerful.MinerFulMinerSlider',
        "-iLF", file_input,
        "-iLStartAt", "0",
        "-iLSubLen", str(window_size),
        "-sliBy", str(sliding_window_size),
        '-para', '4',
        '-s', '0.0',
        '-c', '0.0',
        # '-g', '0.0',
        '-prune', 'none',
        '-sliOut', output_log_path,
        '-vShush'
    ], env=env
        , cwd=os.getcwd())
# , cwd=r"/")

def prune_constraints_minerful(output_constraint_path,output_constraint_path_pruned):
    env = dict(os.environ)
    env['JAVA_OPTS'] = 'foo'
    subprocess.call(['java', "-Xmx16G", '-cp', f'MINERful.jar',
                     'minerful.MinerFulSimplificationStarter',
                     "-iSF",
                     output_constraint_path,
                     "-iSE", 'json',
                     "-oCSV", output_constraint_path_pruned,
                     "-keep",
                     "-prune", "hierarchyconflictredundancy"],
                    env=env,
                    cwd=os.getcwd())

def discover_declare(input_log_path,output_log_path,support, confidence):
    env = dict(os.environ)
    env['JAVA_OPTS'] = 'foo'
    subprocess.call(['java', '-version'])
    file_input = input_log_path
    subprocess.call([
        'java', "-Xmx16G",
        '-cp', f'MINERful.jar',
        'minerful.MinerFulMinerStarter',
        "-iLF", file_input,
        "-s", str(support),
        "-c", str(confidence),
        "-g", "0.0",
        "-sT", "0.00",
        "-cT", "0.00",
        "-gT", "0.0",
        '-prune', 'hierarchy',
        '-oJSON', output_log_path
    ], env=env
        , cwd=os.getcwd())