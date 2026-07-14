import pm4py

def conformance_checking(log_path, model_path):
    log = pm4py.read_xes(str(log_path))
    net,im,fm= pm4py.read.read_pnml(model_path)

    fitness_dict = pm4py.conformance.fitness_alignments(log,net,im,fm)
    precision = pm4py.conformance.precision_alignments(log,net,im,fm)

    return round(fitness_dict["log_fitness"],2),round(precision,2)

def conformance_checking_bi(log_path_Lp,log_path_Lm, model_path):
    Lp = pm4py.read_xes(str(log_path_Lp), variant="rustxes")
    Lm = pm4py.read_xes(str(log_path_Lm), variant="rustxes")

    net,im,fm= pm4py.read.read_pnml(model_path)

    fitness_dict_Lp = pm4py.conformance.fitness_alignments(Lp,net,im,fm)
    fitness_dict_Lm = pm4py.conformance.fitness_alignments(Lm, net, im, fm)

    align_acc = fitness_dict_Lp["log_fitness"] - fitness_dict_Lm["log_fitness"]
    align_F1 = (2*fitness_dict_Lp["log_fitness"]*(1-fitness_dict_Lm["log_fitness"]))/(fitness_dict_Lp["log_fitness"]+(1-fitness_dict_Lm["log_fitness"]))
    trace_acc = (fitness_dict_Lp["percentage_of_fitting_traces"]/100) - (fitness_dict_Lm["percentage_of_fitting_traces"]/100)
    trace_F1 = (2 * (fitness_dict_Lp["percentage_of_fitting_traces"]/100) * (1 - (fitness_dict_Lm["percentage_of_fitting_traces"]/100))) / (
                (fitness_dict_Lp["percentage_of_fitting_traces"]/100) + (1 - (fitness_dict_Lm["percentage_of_fitting_traces"]/100)))

    precision_Lp = pm4py.conformance.precision_alignments(Lp,net,im,fm)

    rep_dict = {"align_fit_Lp": fitness_dict_Lp["log_fitness"],
                "align_fit_Lm": fitness_dict_Lm["log_fitness"],
                "trace_fit_Lp": fitness_dict_Lp["percentage_of_fitting_traces"]/100,
                "trace_fit_Lm": fitness_dict_Lm["percentage_of_fitting_traces"]/100,
                "precision_Lp": precision_Lp,
                "align_acc":align_acc,
                "align_F1":align_F1,
                "trace_acc":trace_acc,
                "trace_F1":trace_F1,
                }

    return rep_dict

def conformance_checking_bi_cl(log_path_Lp,log_path_Lm,lob_path_cl, model_path):
    Lp = pm4py.read_xes(str(log_path_Lp), variant="rustxes")
    Lm = pm4py.read_xes(str(log_path_Lm), variant="rustxes")
    L_cl = pm4py.read_xes(str(lob_path_cl), variant="rustxes")

    net,im,fm= pm4py.read.read_pnml(model_path)

    fitness_dict_Lp = pm4py.conformance.fitness_alignments(Lp,net,im,fm)
    xx = pm4py.conformance.conformance_diagnostics_alignments(Lp,net,im,fm)
    fitness_dict_Lm = pm4py.conformance.fitness_alignments(Lm, net, im, fm)
    fitness_dict_cl = pm4py.conformance.fitness_alignments(L_cl,net,im,fm)

    align_acc = fitness_dict_Lp["log_fitness"] - fitness_dict_Lm["log_fitness"]
    align_F1 = (2*fitness_dict_Lp["log_fitness"]*(1-fitness_dict_Lm["log_fitness"]))/(fitness_dict_Lp["log_fitness"]+(1-fitness_dict_Lm["log_fitness"]))
    trace_acc = (fitness_dict_Lp["percentage_of_fitting_traces"]/100) - (fitness_dict_Lm["percentage_of_fitting_traces"]/100)
    trace_F1 = (2 * (fitness_dict_Lp["percentage_of_fitting_traces"]/100) * (1 - (fitness_dict_Lm["percentage_of_fitting_traces"]/100))) / (
                (fitness_dict_Lp["percentage_of_fitting_traces"]/100) + (1 - (fitness_dict_Lm["percentage_of_fitting_traces"]/100)))

    precision_Lp = pm4py.conformance.precision_alignments(Lp,net,im,fm)
    precision_cl = pm4py.conformance.precision_alignments(L_cl, net, im, fm)

    rep_dict = {r'$a\text{-} fit(L\vert_{r},M)$': fitness_dict_cl["log_fitness"],
                r'$prc(L\vert_{r},M)$': precision_cl,
                # "trace_fit_cl": fitness_dict_cl["percentage_of_fitting_traces"]/100,
                r'$a\text{-} fit(L^+,M)$': fitness_dict_Lp["log_fitness"],
                r'$a\text{-} fit(L^-,M)$': fitness_dict_Lm["log_fitness"],
                r'$t\text{-} fit(L^+,M)$': fitness_dict_Lp["percentage_of_fitting_traces"]/100,
                r'$t\text{-} fit(L^-,M)$': fitness_dict_Lm["percentage_of_fitting_traces"]/100,
                # "precision_Lp": precision_Lp,
                r'$a\text{-} acc(L^+,L^-,M)$':align_acc,
                r'$a\text{-} F1(L^+,L^-,M)$':align_F1,
                r'$t\text{-} acc(L^+,L^-,M)$':trace_acc,
                r'$t\text{-} F1(L^+,L^-,M)$':trace_F1,
                }
    return rep_dict

def extract_significant_dev(dev_list):
    list_dev = []
    for x in dev_list:
        if isinstance(x[1],str):
            list_dev.append((f"{x[0]}({x[1]})",x[2],x[3]))
        else:
            list_dev.append((f"{x[0]}({x[1][0]},{x[1][1]})",x[2],x[3]))
    sorted_IMr = sorted(list_dev, key=lambda x: x[1],reverse=True)
    return sorted_IMr[0:10]

