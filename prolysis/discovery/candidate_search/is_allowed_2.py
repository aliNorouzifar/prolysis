import prolysis.rules_handling.declare_processing as declare_processing


def is_allowed(SS1,SS2,rules_lookup,to_be_checked, lowest_cost):
    # dim = 'support'
    dim = rules_lookup[2]
    rules =  rules_lookup[0]
    lookup_table = rules_lookup[1]
    SS = SS1.union(SS2)


    exclude = set()
    penalty = {'seq': 0, 'exc': 0, 'par': 0, 'loop': 0, 'loop_tau': 0, 'exc_tau': 0}
    for r in rules:
        x =r["template"]
        if len(r['parameters'])==2:
            if (r['parameters'][0][0] not in SS) or (r['parameters'][1][0] not in SS):
                continue
            if SS2:
                y = f"{int(r['parameters'][0][0] in SS1)}{int(r['parameters'][1][0] in SS1)}00"
            else:
                y = f"{int(r['parameters'][0][0] in SS1)}{int(r['parameters'][1][0] in SS1)}10"
        else:
            if (r['parameters'][0][0] not in SS):
                continue
            else:
                if SS2:
                    y = f"{int(r['parameters'][0][0] in SS1)}000"
                else:
                    y = f"{int(r['parameters'][0][0] in SS1)}110"
        exclude_new = lookup_table.loc[x][y]
        # if not SS2 and ('loop_tau' in exclude_new):
        #     print('wait')
        for ex in exclude_new:
            penalty[ex] += r[dim]
        exclude = exclude | exclude_new

        if all([penalty[ct]>lowest_cost for ct in to_be_checked]):
            # print(f'visited')
            break
    # print(to_be_checked)
    lowest_cost = min(lowest_cost,min([penalty[ct] for ct in to_be_checked]))
    # print(f'{lowest_cost}')
    return exclude, set(), penalty,lowest_cost



def is_allowed_single(activity,rules_lookup):
    rules =  rules_lookup[0]
    lookup_table = rules_lookup[1]
    # dim = 'support'
    dim = rules_lookup[2]
    exclude = set()
    penalty = {'single_single': 0, 'xor_single': 0, 'loop_single': 0}
    for r in rules:
        if len(r['parameters']) == 2 or activity!=r['parameters'][0][0]:
            continue
        elif activity==r['parameters'][0][0]:
            x =r["template"]
            y = f"0001"
            exclude_new = lookup_table.loc[x][y]
            for ex in exclude_new:
                penalty[ex] += r[dim]
            exclude = exclude | exclude_new
            # if penalty>=low_pen:
            #     break
    return exclude, penalty

