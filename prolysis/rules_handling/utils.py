import json
import prolysis.rules_handling.declare_processing as declare_processing


def rules_from_json(file_path):
    with open(file_path, 'r') as file:
        declare_file = json.load(file)
    return declare_file['constraints'],set(declare_file['tasks'])
def preprocess(rules):
    rules_processes = []
    co_exist_list = []
    absence_list = set([r['parameters'][0][0] for r in rules if r['template'] == "Absence"])
    for r in rules:
        if r['template'] == 'AtLeast2' or r['template'] == 'AtLeast3':
            r_new = r.copy()
            r_new['template'] = 'AtLeast1'
            rules_processes.append(r_new)
        elif r['template'] =="Absence":
            continue
        elif r['template'] == 'AtMost2' or r['template'] == 'AtMost3':
            continue
        elif r['template'] == 'CoExistence':
            co_exist_list.append((r['parameters'][0][0],r['parameters'][1][0]))
            if (r['parameters'][1][0],r['parameters'][0][0]) not in co_exist_list:
                rules_processes.append(r)
        else:
            rules_processes.append(r)
    return rules_processes,absence_list

def dfa_list_generator(rules,S_mapped):
    dfa_list = []
    for c in rules:
        if "Absence" in c['template'] or "Init" in c['template'] or "End" in c['template'] or "AtMost" in c[
            'template'] or "AtLeast" in c['template']:
            dfa_to_add = declare_processing.gen_reg_dfa(c['template'], [c['parameters'][0][0]], S_mapped)
            dfa_list.append(((c['template'], (c['parameters'][0][0])), dfa_to_add, c['support'],c['confidence']))
        else:
            dfa_to_add = declare_processing.gen_reg_dfa(c['template'],
                                                        [c['parameters'][0][0], c['parameters'][1][0]],
                                                        S_mapped)
            dfa_list.append((
                            (c['template'], (c['parameters'][0][0], c['parameters'][1][0])),
                            dfa_to_add, c['support'],c['confidence']))

    total_sup = sum([el[2] for el in dfa_list])
    total_conf = sum([el[3] for el in dfa_list])
    return dfa_list, total_sup, total_conf

