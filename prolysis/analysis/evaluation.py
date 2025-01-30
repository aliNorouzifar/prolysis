import pm4py

def conformance_checking(log_path, model_path):
    log = pm4py.read_xes(str(log_path))
    net,im,fm= pm4py.read.read_pnml(model_path)

    fitness_dict = pm4py.conformance.fitness_alignments(log,net,im,fm)
    precision = pm4py.conformance.precision_alignments(log,net,im,fm)

    return round(fitness_dict["log_fitness"],2),round(precision,2)

def extract_significant_dev(dev_list):
    list_dev = []
    for x in dev_list:
        if isinstance(x[1],str):
            list_dev.append((f"{x[0]}({x[1]})",x[2],x[3]))
        else:
            list_dev.append((f"{x[0]}({x[1][0]},{x[1][1]})",x[2],x[3]))
    sorted_IMr = sorted(list_dev, key=lambda x: x[1],reverse=True)
    return sorted_IMr[0:10]

