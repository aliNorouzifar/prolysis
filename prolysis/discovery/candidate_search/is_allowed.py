import prolysis.rules_handling.declare_processing as declare_processing


def is_allowed(SS1,SS2,rules,low_pen,to_be_checked):
    dfa_intersect = rules[0]

    da_list = rules[1]
    S_mapped = rules[2]
    exclude = []
    S1 = {rules[2][x] for x in SS1}
    S2 = {rules[2][x] for x in SS2}
    penalty = {'seq': 0, 'exc': 0, 'par': 0, 'loop': 0, 'loop_tau': 0, 'exc_tau':0}
    for cut_type in to_be_checked:
        # if (declare_processing.check_all(S1, S2, dfa_intersect, cut_type, S_mapped)):
        # if (
        # declare_processing.check_all(S1, S2, S1 , S2 , dfa_intersect, cut_type,
        #                              S_mapped)):
        #     # print(dfa[0])
        #     exclude.append(cut_type)
        #     penalty[cut_type] += +1
        #
        for dfa in da_list:
            if (declare_processing.check_all(S1,S2,S1&{x for x in dfa[0][1]}, S2&{x for x in dfa[0][1]}, dfa[1], cut_type, S_mapped)):

                # print(dfa[0])
                exclude.append(cut_type)
                if penalty[cut_type]+dfa[2]>low_pen:
                    penalty[cut_type] += dfa[2]
                    break
                else:
                    penalty[cut_type] += dfa[2]
                    low_pen = min(low_pen,penalty[cut_type])
        if penalty[cut_type]==0:
            low_pen = 0

    # if sum(penalty.values())>0:
    #     print(penalty)
    # print(low_pen)
    return set(exclude), set(), penalty,low_pen


def is_allowed_single(activity, rules):
    penalty = {'single_single': 0, 'xor_single': 0, 'loop_single': 0}

    da_list = rules[1]
    S_mapped = rules[2]
    exclude = []

    for cut_type in {'single_single', 'xor_single', 'loop_single'}:
        for dfa in da_list:
            # if ("Absence" in dfa[0][0] or "Init" in dfa[0][0] or "End" in dfa[0][0] or "AtMost" in dfa[0][0] or "AtLeast" in dfa[0][0]):
            if dfa[0][1]==S_mapped[activity]:
                if (declare_processing.check_all_single(S_mapped[activity], dfa[1],
                                                 cut_type, S_mapped)):
                    print(dfa[0])
                    exclude.append(cut_type)
                    penalty[cut_type] += dfa[2]

    return set(exclude), penalty