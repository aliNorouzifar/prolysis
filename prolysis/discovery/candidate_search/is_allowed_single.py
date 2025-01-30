import functions.functions.declare_processing as declare_processing


def is_allowed(activity,rules,low_pen):
    penalty = {'single':0, 'xor':0, 'loop':0}

    da_list = rules[1]
    S_mapped = rules[2]
    exclude = []


    for cut_type in {'single', 'xor', 'loop'}:
        for dfa in da_list:
            if (declare_processing.check_all(S1,S2,S1&{x for x in dfa[0][1]}, S2&{x for x in dfa[0][1]}, dfa[1], cut_type, S_mapped)):
                print(dfa[0])
                exclude.append(cut_type)
                if penalty[cut_type]+dfa[2]>low_pen:
                    penalty[cut_type] += dfa[2]
                    break
                else:
                    penalty[cut_type] += dfa[2]
                    low_pen = min(low_pen,penalty[cut_type])
    # if sum(penalty.values())>0:
    #     print(penalty)
    print(low_pen)
    return set(exclude), set(), penalty,low_pen
