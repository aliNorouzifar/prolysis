from automata.fa.nfa import NFA
from automata.fa.dfa import DFA
from frozendict import frozendict

def extract_templates():
    return{
        "AtLeast1": (lambda a, S: f"{g(S-{a})}*({a}{g(S-{a})}*){{1,}}{g(S-{a})}*",1),
        "AtLeast2": (lambda a, S: f"{g(S-{a})}*({a}{g(S-{a})}*){{2,}}{g(S-{a})}*",1),
        "AtLeast3": (lambda a, S: f"{g(S-{a})}*({a}{g(S-{a})}*){{3,}}{g(S-{a})}*",1),
        "AtMost1": (lambda a, S: f"{g(S-{a})}*({a}{g(S-{a})}*){{0,1}}{g(S-{a})}*",1),
        "AtMost2": (lambda a, S: f"{g(S-{a})}*({a}{g(S-{a})}*){{0,2}}{g(S-{a})}*",1),
        "AtMost3": (lambda a, S: f"{g(S-{a})}*({a}{g(S-{a})}*){{0,3}}{g(S-{a})}*",1),
        "Absence": (lambda a, S: f"{g(S-{a})}*",1),
        "Init": (lambda a, S: f"{a}.*",1),
        "End": (lambda a, S: f".*{a}",1),
        "CoExistence": (lambda a, b, S: f"{g(S-{a,b})}*(({a}.*{b}.*)|({b}.*{a}.*))*{g(S-{a,b})}*",2),
        "RespondedExistence": (lambda a, b, S: f"{g(S-{a})}*(({a}.*{b}.*)|({b}.*{a}.*))*{g(S-{a})}*",2),
        "Response": (lambda a, b, S: f"{g(S-{a})}*({a}.*{b})*{g(S-{a})}*",2),
        "AlternateResponse": (lambda a, b, S: f"{g(S-{a})}*({a}{g(S-{a})}*{b}{g(S-{a})}*)*{g(S-{a})}*",2),
        "ChainResponse": (lambda a, b, S: f"{g(S-{a})}*({a}{b}{g(S-{a})}*)*{g(S-{a})}*",2),
        "Precedence": (lambda a, b, S: f"{g(S-{b})}*({a}.*{b})*{g(S-{b})}*",2),
        "AlternatePrecedence": (lambda a, b, S: f"{g(S-{b})}*({a}{g(S-{b})}*{b}{g(S-{b})}*)*{g(S-{b})}*",2),
        "ChainPrecedence": (lambda a, b, S: f"{g(S-{b})}*({a}{b}{g(S-{b})}*)*{g(S-{b})}*",2),
        "Succession": (lambda a, b, S: f"{g(S-{a,b})}*({a}.*{b})*{g(S-{a,b})}*",2),
        "AlternateSuccession": (lambda a, b, S: f"{g(S-{a,b})}*({a}{g(S-{a,b})}*{b}{g(S-{a,b})}*)*{g(S-{a,b})}*",2),
        "ChainSuccession": (lambda a, b, S: f"{g(S-{a,b})}*({a}{b}{g(S-{a,b})}*)*{g(S-{a,b})}*",2),
        "NotRespondedExistence": (lambda a, b, S: f"{g(S-{a,b})}*(({a}{g(S-{b})}*)|({b}{g(S-{a})}*))?",2),
        "NotCoExistence": (lambda a, b, S: f"{g(S-{a,b})}*(({a}{g(S-{b})}*)|({b}{g(S-{a})}*))?",2),
        "NotResponse": (lambda a, b, S: f"{g(S-{a})}*({a}{g(S-{b})}*)*{g(S-{a,b})}*",2),
        "NotChainResponse": (lambda a, b, S: f"{g(S-{a})}*({a}{a}*{g(S-{a,b})}{g(S-{a})}*)*({g(S-{a})}*|{a}*)",2),
        "NotPrecedence": (lambda a, b, S: f"{g(S-{a})}*({a}{g(S-{b})}*)*{g(S-{a,b})}*",2),
        "NotChainPrecedence": (lambda a, b, S: f"{g(S-{a})}*({a}{a}*{g(S-{a,b})}{g(S-{a})}*)*({g(S-{a})}*|{a}*)",2),
        "NotSuccession": (lambda a, b, S: f"{g(S-{a})}*({a}{g(S-{b})}*)*{g(S-{a,b})}*",2),
        # "NotChainSuccession": (lambda a, b, S: f"{g(S-{b})}*({b}{b}*{g(S-{a,b})}{g(S-{b})}*)*({g(S-{b})}*|{b}*)",2),
        "NotChainSuccession": (lambda a, b, S: f"{g(S-{a})}*({a}{a}*{g(S-{a,b})}{g(S-{a})}*)*({g(S-{a})}*|{a}*)",2),
    }
templates = extract_templates()

def g(input_set):
    """
    Converts a set of strings into a single string formatted as 'a|b|c'.

    Parameters:
    - input_set (set): A set of strings.

    Returns:
    - str: A single string in the format 'a|b|c'.
    """
    if not input_set:
        return ""  # Return an empty string for an empty set
    return f'({"|".join(sorted(input_set))})'

def gen_reg_dfa(template,parameters,S_mapped):
    S = set(S_mapped.values())
    if template in templates:

        if "Absence" in template or "Init" in template or "End" in template or "AtMost" in template or "AtLeast" in template:
            a = S_mapped[parameters[0]]
            reg = templates[template][0](a, S)
        else:
            a, b = S_mapped[parameters[0]], S_mapped[parameters[1]]
            reg = templates[template][0](a, b, S)
        # print(f"Generated regex for {template}: {reg}")
        nfa = NFA.from_regex(reg, input_symbols=S)
        dfa = rename_dfa_transitions(DFA.from_nfa(nfa), S_mapped)
        # Make the DFA complete
        complete_dfa = make_dfa_complete(dfa)
        return complete_dfa
    else:
        return ValueError("Invalid template")

def rename_dfa_transitions(dfa, mapping):
    """
    Renames the transition symbols in a DFA using a provided mapping dictionary.

    Parameters:
    - dfa (DFA): The input DFA whose transitions need renaming.
    - mapping (dict): A dictionary mapping current symbols to new symbols.

    Returns:
    - DFA: A new DFA with renamed transitions.
    """
    from frozendict import frozendict
    # from automata.fa.dfa import DFA
    mapping = {value: key for key, value in mapping.items()}
    # Rename the transitions
    renamed_transitions = {}
    for state, transition_map in dfa.transitions.items():
        renamed_transitions[state] = {
            mapping.get(symbol, symbol): dest_state  # Rename symbols based on mapping
            for symbol, dest_state in transition_map.items()
        }

    # Create a new DFA with renamed transitions
    return DFA(
        states=frozenset(dfa.states),
        input_symbols=frozenset(mapping.get(symbol, symbol) for symbol in dfa.input_symbols),
        transitions=frozendict({
            state: frozendict(trans) for state, trans in renamed_transitions.items()
        }),
        initial_state=dfa.initial_state,
        final_states=frozenset(dfa.final_states), allow_partial=True
    )

def make_dfa_complete(dfa):
    """
    Ensure the given DFA is complete by adding a dead state for missing transitions.
    """
    # Convert frozenset of states to a mutable set
    states = set(dfa.states)  # Convert frozenset to set for modification
    input_symbols = dfa.input_symbols  # Input symbols remain the same
    transitions = {state: dict(trans) for state, trans in dfa.transitions.items()}  # Convert frozendict to dict
    initial_state = dfa.initial_state  # Initial state remains the same
    final_states = dfa.final_states  # Final states remain the same

    # Add a dead state
    dead_state = "dead"
    states.add(dead_state)

    # Add missing transitions for each state
    for state in states:
        if state not in transitions:
            transitions[state] = {}
        for symbol in input_symbols:
            if symbol not in transitions[state]:
                transitions[state][symbol] = dead_state

    # Add self-loops on the dead state
    transitions[dead_state] = {symbol: dead_state for symbol in input_symbols}

    # Return the completed DFA with the updated states and transitions
    return DFA(
        states=frozenset(states),  # Convert back to frozenset for DFA consistency
        input_symbols=input_symbols,
        transitions=frozendict({state: frozendict(trans) for state, trans in transitions.items()}),  # Re-freeze transitions
        initial_state=initial_state,
        final_states=final_states,
    )


def assign_alphabet(S):
    """
    Assigns a unique alphabet character to each string in the list S.

    Parameters:
    - S (list): List of strings.

    Returns:
    - dict: A dictionary mapping each string in S to a unique alphabet character.
    """
    # List of single-character alphabet symbols
    alphabet = [chr(i) for i in range(ord('a'), ord('z') + 1)]  # 'a' to 'z'

    if len(S) > len(alphabet):
        raise ValueError("Not enough unique alphabet characters to assign to all strings in S.")

    # Map each string in S to a unique alphabet character
    mapping = {string: alphabet[i] for i, string in enumerate(S)}
    return mapping

def symbol_extract(tuple_str):
    tuple_str = tuple_str.strip("()")
    elements = tuple_str.split(", ")
    second_element = elements[1].strip("'")
    return second_element

def reachability2NFA(rg,activities):
    start_state = [s.name for s in rg.states if not s.incoming][0]
    final_state = set([s.name for s in rg.states if not s.outgoing])
    states = set([s.name for s in rg.states])
    input_symbols = activities
    # input_symbols = set()
    # for t in rg.transitions:
    #     if t.name=='None':
    #         input_symbols.add("")
    #     else:
    #         input_symbols.add(symbol_extract(t.name))
    transitions = {}
    for state in rg.states:
        tr_to_state = {}
        for tr in state.outgoing:
            if symbol_extract(tr.name) == 'None':
                if "" in tr_to_state:
                    tr_to_state[""].add(tr.to_state.name)
                else:
                    tr_to_state[""] = {tr.to_state.name}
            else:
                tr_to_state[symbol_extract(tr.name)] = {tr.to_state.name}
        transitions[state.name] = tr_to_state

    nfa = NFA(
            states=states,
            input_symbols=input_symbols,
            transitions=transitions,
            initial_state=start_state,
            final_states=final_state,
        )
    return nfa

def check_all(S1,S2,SS1, SS2, dfa_decl, cut_type, S_mapped):
    flag = False
    if cut_type == 'seq':
        for p in SS1:
            for q in SS2:
                reg_to_check1 = f"{g(S1 - {p})}*{p}{g(S1 - {p})}*{g(S2 - {q})}*{q}{g(S2 - {q})}*"
                nfa1 = NFA.from_regex(reg_to_check1, input_symbols=set(S_mapped.values()))
                dfa_to_check1 = rename_dfa_transitions(DFA.from_nfa(nfa1),S_mapped)
                if dfa_to_check1.issubset(dfa_decl.complement()):
                    flag = True
                    # block = True
                    break

    elif cut_type == 'exc':
        for p in SS1:
            reg_to_check2 = f"{g(S1 - {p})}*{p}{g(S1 - {p})}*"
            nfa2 = NFA.from_regex(reg_to_check2, input_symbols=set(S_mapped.values()))
            dfa_to_check2 = rename_dfa_transitions(DFA.from_nfa(nfa2),S_mapped)
            if dfa_to_check2.issubset(dfa_decl.complement()):
                flag = True
                # block = True
                break

        for q in SS2:
            reg_to_check3 = f"{g(S2 - {q})}*{q}{g(S2 - {q})}*"
            nfa3 = NFA.from_regex(reg_to_check3, input_symbols=set(S_mapped.values()))
            dfa_to_check3 = rename_dfa_transitions(DFA.from_nfa(nfa3),S_mapped)
            if dfa_to_check3.issubset(dfa_decl.complement()):
                # print(f"q is {q}")
                flag = True
                break



    elif cut_type == 'par':
        for p in SS1:
            for q in SS2:
                reg_to_check4 = f"{g(S1.union(S2) - {p, q})}*{p}{g(S1.union(S2) - {p, q})}*{q}{g(S1.union(S2) - {p, q})}*"
                reg_to_check5 = f"{g(S1.union(S2) - {p, q})}*{q}{g(S1.union(S2) - {p, q})}*{p}{g(S1.union(S2) - {p, q})}*"
                nfa4 = NFA.from_regex(reg_to_check4, input_symbols=set(S_mapped.values()))
                dfa_to_check4 = rename_dfa_transitions(DFA.from_nfa(nfa4),S_mapped)
                nfa5 = NFA.from_regex(reg_to_check5, input_symbols=set(S_mapped.values()))
                dfa_to_check5 = rename_dfa_transitions(DFA.from_nfa(nfa5),S_mapped)
                if dfa_to_check4.issubset(dfa_decl.complement()):
                    flag = True
                    # block = True
                    break
                if dfa_to_check5.issubset(dfa_decl.complement()):
                    flag = True
                    break

    elif cut_type == 'loop':
        for p in SS1:
            reg_to_check6 = f"{g(S1)}*{p}{g(S1)}*"
            nfa6 = NFA.from_regex(reg_to_check6, input_symbols=set(S_mapped.values()))
            dfa_to_check6 = rename_dfa_transitions(DFA.from_nfa(nfa6), S_mapped)
            if dfa_to_check6.issubset(dfa_decl.complement()):
                flag = True
                break
            for pp in SS1-{p}:
                reg_to_check9 = f"({g(S1-{p})}*{p}{g(S1-{p})}*)({g(S2)}*({g(S1-{pp})}*{pp}{g(S1-{pp})}))+"
                nfa9 = NFA.from_regex(reg_to_check9, input_symbols=set(S_mapped.values()))
                dfa_to_check9 = rename_dfa_transitions(DFA.from_nfa(nfa9), S_mapped)
                if dfa_to_check9.issubset(dfa_decl.complement()):
                    flag = True
                    break
            for q in SS2:
                reg_to_check7 = f"({g(S1-{p})}*{p}{g(S1-{p})}*)(({g(S2-{q})}*{q}{g(S2-{q})}*){g(S1)}*)+"
                nfa7 = NFA.from_regex(reg_to_check7, input_symbols=set(S_mapped.values()))
                dfa_to_check7 = rename_dfa_transitions(DFA.from_nfa(nfa7),S_mapped)
                if dfa_to_check7.issubset(dfa_decl.complement()):
                    flag = True
                    break
                reg_to_check8 = f"{g(S1)}*(({g(S2 - {q})}*{q}{g(S2 - {q})}*)({g(S1-{p})}*{p}{g(S1-{p})}*))+"
                nfa8 = NFA.from_regex(reg_to_check8, input_symbols=set(S_mapped.values()))
                dfa_to_check8 = rename_dfa_transitions(DFA.from_nfa(nfa8), S_mapped)
                if dfa_to_check8.issubset(dfa_decl.complement()):
                    flag = True
                    break
                for qq in SS1 - {q}:
                    reg_to_check10 = f"{g(S1)}*({g(S2-{q})}*{q}{g(S2-{q})}*){g(S1)}*(({g(S2-{qq})}*{qq}{g(S2-{qq})}*){g(S1)}*)+"
                    nfa10 = NFA.from_regex(reg_to_check10, input_symbols=set(S_mapped.values()))
                    dfa_to_check10 = rename_dfa_transitions(DFA.from_nfa(nfa10), S_mapped)
                    if dfa_to_check10.issubset(dfa_decl.complement()):
                        flag = True
                        break



    elif cut_type == 'loop_tau':
        for p in SS1:
            # p = S_mapped[pp]
            reg_to_check = f"{g(S1)}*{p}{g(S1)}"
            reg_to_check2 = f"({g(S1)}*{p}{g(S1)}){{2,}}"

            # reg_to_check2 = f"{g(S1)}*(({g(S2)}*{q}{g(S2)}*)({g(S1)}*{p}{g(S1)}*)){{2,}}"
            nfa = NFA.from_regex(reg_to_check, input_symbols=set(S_mapped.values()))
            dfa_to_check = rename_dfa_transitions(DFA.from_nfa(nfa),S_mapped)
            if dfa_to_check.issubset(dfa_decl.complement()):
                flag = True
                # block = True
                break
            nfa2 = NFA.from_regex(reg_to_check2, input_symbols=set(S_mapped.values()))
            dfa_to_check2 = rename_dfa_transitions(DFA.from_nfa(nfa2), S_mapped)
            if dfa_to_check2.issubset(dfa_decl.complement()):
                flag = True
                # block = True
                break
    elif cut_type == 'exc_tau':
        for p in SS1:
            reg_to_check = f"{g(S1 - {p})}*{p}{g(S1 - {p})}*"
            nfa = NFA.from_regex(reg_to_check, input_symbols=set(S_mapped.values()))
            dfa_to_check = rename_dfa_transitions(DFA.from_nfa(nfa), S_mapped)
            if dfa_to_check.issubset(dfa_decl.complement()):
                # print(f"p is {p}")
                # print(reg_to_check)
                flag = True
                # block = True
                break

        reg_to_check2 = ""
        nfa2 = NFA.from_regex(reg_to_check2, input_symbols=set(S_mapped.values()))
        dfa_to_check2 = rename_dfa_transitions(DFA.from_nfa(nfa2), S_mapped)
        if dfa_to_check2.issubset(dfa_decl.complement()):
            # print(f"p is {p}")
            # print(reg_to_check)
            flag = True
    return flag


def check_all_single(activity, dfa_decl, cut_type, S_mapped):
    flag = False
    if cut_type == 'single_single':
        reg_to_check = f"{activity}"
        nfa = NFA.from_regex(reg_to_check, input_symbols=set(S_mapped.values()))
        dfa_to_check = rename_dfa_transitions(DFA.from_nfa(nfa), S_mapped)
        if dfa_to_check.issubset(dfa_decl.complement()):
            flag = True
    elif cut_type == 'xor_single':
        reg_to_check = f"{activity}"
        nfa = NFA.from_regex(reg_to_check, input_symbols=set(S_mapped.values()))
        dfa_to_check = rename_dfa_transitions(DFA.from_nfa(nfa), S_mapped)
        if dfa_to_check.issubset(dfa_decl.complement()):
            flag = True
        reg_to_check2 = ""
        nfa2 = NFA.from_regex(reg_to_check2, input_symbols=set(S_mapped.values()))
        dfa_to_check2 = rename_dfa_transitions(DFA.from_nfa(nfa2), S_mapped)
        if dfa_to_check2.issubset(dfa_decl.complement()):
            # print(f"p is {p}")
            # print(reg_to_check)
            flag = True
    elif cut_type == 'loop_single':
        reg_to_check = f"{activity}"
        nfa = NFA.from_regex(reg_to_check, input_symbols=set(S_mapped.values()))
        dfa_to_check = rename_dfa_transitions(DFA.from_nfa(nfa), S_mapped)
        if dfa_to_check.issubset(dfa_decl.complement()):
            flag = True
        reg_to_check2 = f"{activity}{activity}"
        nfa2 = NFA.from_regex(reg_to_check2, input_symbols=set(S_mapped.values()))
        dfa_to_check2 = rename_dfa_transitions(DFA.from_nfa(nfa2), S_mapped)
        if dfa_to_check2.issubset(dfa_decl.complement()):
            flag = True
        reg_to_check3 = f"{activity}{activity}{activity}"
        nfa3 = NFA.from_regex(reg_to_check3, input_symbols=set(S_mapped.values()))
        dfa_to_check3 = rename_dfa_transitions(DFA.from_nfa(nfa3), S_mapped)
        if dfa_to_check3.issubset(dfa_decl.complement()):
            flag = True


    return flag


