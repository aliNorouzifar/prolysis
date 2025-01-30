import os
import subprocess


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