import medicc
import fstlib


# define useful variables
separator = 'X'
max_cn = 8
SYMBOL_TABLE = medicc.create_symbol_table(max_cn=max_cn, separator=separator)
SYMBOL_DICT = dict([(i, s) for i, s in SYMBOL_TABLE])


def get_fst_events(s1fsa, s2fsa, pre_fst, post_fst=None, path_limit=1000):
    if post_fst is not None:
        tmp = fstlib.intersect((s1fsa*pre_fst).project('output').arcsort('olabel'), (post_fst * s2fsa).project('input'))
    else:
        tmp = fstlib.intersect((s1fsa*pre_fst).project('output').arcsort('olabel'), (s2fsa).project('input'))
    cur = fstlib.disambiguate(fstlib.prune(tmp, weight=0))
    cur = fstlib.shortestpath(cur, nshortest=path_limit)
    return cur


def create_1step_forced_WGD_fst(symbol_table, separator='X', wgd_cost=1, minimize=True, wgd_x2=False,
                                total_cn=False):

    if total_cn:
        wgd_distance = 2.
    else:
        wgd_distance = 1.

    cns = medicc.factory._get_int_cns_from_symbol_table(symbol_table, separator)

    W = fstlib.Fst()
    W.set_input_symbols(symbol_table)
    W.set_output_symbols(symbol_table)
    W.add_states(2)
    W.set_start(0)
    W.set_final(0, 0)
    W.set_final(1, 0)
    # W.set_final(2, 0)
    if wgd_x2:
        W.add_arcs(0, [(s, t, wgd_cost, 1) for s in cns.keys()
                       for t in cns.keys() if (s != '0') and ((cns[t]/cns[s]) == 2.)])
    else:
        W.add_arcs(0, [(s, t, wgd_cost, 1) for s in cns.keys()
                       for t in cns.keys() if (s != '0') and ((cns[t]-cns[s]) == wgd_distance)])

    W.add_arc(1, ('0', '0', 0, 1))
    if wgd_x2:
        W.add_arcs(1, [(s, t, 0, 1) for s in cns.keys()
                       for t in cns.keys() if (s != '0') and ((cns[t]/cns[s]) == 2.)])
    else:
        W.add_arcs(1, [(s, t, 0, 1) for s in cns.keys()
                       for t in cns.keys() if (s != '0') and ((cns[t]-cns[s]) == wgd_distance)])
    W.add_arc(0, ('0', '0', 0, 0))
    if separator is not None and separator != '':
        W.add_arc(0, (separator, separator, 0, 0))
        W.add_arc(1, (separator, separator, 0, 1))
    W.add_arc(1, ('0', '0', 0, 1))

    # W.add_arcs(0, [(s, s, 0, 2) for s in cns.keys() if s != '0'])
    # W.add_arcs(2, [(s, s, 0, 2) for s in cns.keys()])
    # W.add_arc(2, ('0', '0', 0, 2))
    # if separator is not None and separator != '':
    #     W.add_arc(2, (separator, separator, 0, 2))
    W.arcsort('olabel')
    if minimize:
        W = fstlib.encode_determinize_minimize(W)

    return W


def fsa_from_string(cur_str, symbol_table=SYMBOL_TABLE):
    return fstlib.factory.from_string(
        cur_str, isymbols=symbol_table, osymbols=symbol_table, arc_type=fstlib.Semiring.TROPICAL)


def fsa_to_string(cur, start=None, symbol_dict=SYMBOL_DICT):

    if start is None:
        start = cur.start()
    cur_state = start
    out = []
    while float(cur.final(cur_state)) != 0:
        arc = next(cur.arcs(cur_state))
        cur_label = symbol_dict[arc.olabel]
        cur_state = arc.nextstate
        if cur_label == 'gap':
            continue
        out.append(cur_label)

    return ''.join(out)


def fsa_to_string_old(fsa, start=None):
    if start is None:
        start = fsa.start()
    fsa_df = fsa.to_dataframe()
    cur_i = start
    output = ''
    for _ in range(len(fsa_df)):
        cur_row = fsa_df.loc[fsa_df['state_from'] == cur_i]
        if len(cur_row) == 0:
            break
        cur_i = cur_row['state_to'].values[0]
        if cur_row['ilabel'].values[0] == 'gap':
            continue
        output += cur_row['ilabel'].values[0]
    return output


def aggregate_copy_number_profile(cnp, separator='X'):
    return separator.join([separator.join(["".join(x.astype('str'))
                                           for _, x in cnp[allele].groupby('chrom')]) for allele in cnp.columns])


def get_all_paths_from_fst(fst, symbol_table=SYMBOL_TABLE):
    starting_nodes = [arc.nextstate for arc in fst.arcs(fst.start())]
    paths = []
    for node in starting_nodes:
        output = fsa_to_string(fst, start=node)
        fsa = fstlib.factory.from_string(output, isymbols=symbol_table, osymbols=symbol_table, arc_type=fstlib.Semiring.TROPICAL)
        paths.append((output, fsa))
    return paths
