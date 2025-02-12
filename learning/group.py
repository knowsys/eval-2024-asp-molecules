import sys
import re
import pprint
import networkx as nx
import bisect
import random

random.seed(42) # Set the random number generator to a fixed sequence

check_arg = lambda short, long: len(sys.argv) > 1 and ('-' + short in sys.argv[1:] or '--' + long in sys.argv[1:])

if check_arg("h", "help"):
    exit("Usage: clingo 0 input.lp <params> | python3.10 group.py {lex-leader,min-adj} [-h | -c | -d] 2> >(tee examples.lasp)\n\n"
       + "Convert answer-sets of logic program into Networkx graphs and find isomprphism classes\n"
       + "Mark the chosen representative of each isomorphism class with #pos and the others with #neg\n\n"
       + "criterium='lex-leader' --> choose lexicographic leader as representative\n"
       + "criterium='min-adj'    --> choose graph with minimum adjacency matrix as representative (default)\n\n"
       + "Options:\n"
       + "-h, --help   : Print this help message and exit\n"
       + "-c, --check  : Check for single-swap pruning correctness (based on `prune/3` predicate)"
       + "-d, --debug  : Print the models\n"
       + "<params>     : Either supply sum formula using c, h, o, n via Clingo's --const option,\n"
       + "             : or use Genmol's to-factbase command to generate input\n\n"
       + "Examples:\n"
       + "clingo 0 input.lp --const c=2 --const h=6 --const o=1 | python3.10 group.py lex-leader 2> >(tee lex-leader/examples.lasp)\n"
       + "clingo 0 input.lp --const c=2 --const h=6 --const o=1 | python3.10 group.py min-adj 2> >(tee min-adj/examples.lasp)\n"
       + "clingo 0 input.lp <(genmol to-factbase -f C2H6O) | python3.10 group.py lex-leader -d 2> >(tee lex-leader/examples.lasp)\n"
       + "clingo 0 input.lp <(genmol to-factbase -f C2H6O) | python3.10 group.py min-adj -d 2> >(tee min-adj/examples.lasp)\n"
       + "clingo 0 input.lp lex-leader/lex.lp --const c=5 --const h=10 | python3.10 group.py lex-leader -c\n"
       + "clingo 0 input.lp min-adj/lex.lp --const c=5 --const h=10 | python3.10 group.py min-adj -c")

criterium = next(filter(lambda x: not x.startswith('-'), sys.argv[1:]), 'min-adj')
if criterium not in ['lex-leader', 'min-adj']:
    exit(f'Invalid criterium "{criterium}", choose from {{lex-leader,min-adj}}.')

DEBUG = check_arg("d", "debug")
CHECK = check_arg("c", "check")

isomorphism_classes = []

try:
    i = 0
    for line in sys.stdin:
        if "type" in line:
            i += 1

            model = dict()
            for symbol in line.split():
                group = re.match("^([^(]+)\\(([^)]*)\\)$", symbol)
                pred = group[1]
                args = group[2].split(",")
                if pred not in model.keys():
                    model[pred] = []
                model[pred].append([arg[1:-1] if arg.startswith('"') else int(arg) if arg[0].isnumeric() else arg for arg in args])

            def adj(sorted_model, linebreak=False):
                if len(sorted_model[0][1]) == 0:
                    return '0'
                last_i = 1
                last_j = 0
                w = sorted_model[0][1][-1][0]
                mat = ''
                for k, (i, j, m) in enumerate(sorted_model[0][1]):
                    #print(f'w={w}, k={k}, i={i}, j={j}, m={m}, last_i={last_i}, last_j={last_j}', end='')
                    if i > last_i:
                        last_i = i
                        #print(f', w-last_j={w-last_j}', end='')
                        mat += '0'*(w-last_j) + ('\n' if linebreak else '')
                        last_j = 0
                    #print(f', j-last_j-1={j-last_j-1}', end='')
                    mat += '0'*(j-last_j-1) + str(m)
                    last_j = j
                    if k == len(sorted_model[0][1])-1:
                        #print(f', w-last_j={w-last_j}', end='')
                        mat += '0'*(w-last_j)
                return mat

            if "edge" not in model:
                model["edge"] = []
            sorted_model = [(x[0], sorted(x[1])) for x in sorted(list({"edge": model["edge"]}.items()))]
            order = sorted_model if criterium=='lex-leader' else adj(sorted_model)

            if DEBUG:
                print(f"{i}:")
                pprint.pp(model)
                pprint.pp(sorted_model)
                if criterium=='min-adj':
                    print(adj(sorted_model, linebreak=True))
            elif not CHECK:
                print(f'{i}:', order)

            G = nx.MultiGraph()
            G.add_nodes_from([(typ[0], {"color": typ[1]}) for typ in model["type"]])
            G.add_edges_from(sum([[(edge[0], edge[1])] * edge[2] for edge in model["edge"]], []))
            if DEBUG:
                print(f"{G}: {G.edges()}")
            is_isomorphic = False
            for iso in isomorphism_classes:
                if nx.is_isomorphic(G, iso[0][0], node_match=lambda n1, n2: n1['color']==n2['color']):
                    if not CHECK:
                        print(f"--> isomorphic to {iso[0][2]}")
                    bisect.insort(iso, (G, order, i, model), key=lambda x: x[1])
                    is_isomorphic = True
            if not is_isomorphic:
                isomorphism_classes.append([(G, order, i, model)])

        elif not line.startswith("Answer"):
            print(line, end="")

    print(f"\n{len(isomorphism_classes)} isomorphy classes")
    if DEBUG:
        pprint.pp(isomorphism_classes)

    cnt_pos_pruned = 0
    cnt_pos_not_pruned = 0
    cnt_pos_prunable = 0
    cnt_neg_pruned = 0
    cnt_neg_not_pruned = 0
    cnt_neg_prunable = 0
    for j, iso in enumerate(isomorphism_classes):
        selected_neg_example_k = int(random.uniform(1, len(iso)))
        for k, (_G, order, i, model) in enumerate(iso):
            c = sum(map(lambda x: x[1]=='C', model['type']))
            n = sum(map(lambda x: x[1]=='N', model['type']))
            o = sum(map(lambda x: x[1]=='O', model['type']))
            h = 4*c+3*n+2*o - sum(map(lambda e: e[2], model["edge"]))
            if CHECK:
                if k == 0:
                    print("----")
                    if "prune" in model:
                        cnt_pos_pruned += 1
                        print("ERROR:", order, "-->", model["prune"])
                    else:
                        cnt_pos_not_pruned += 1
                        print(order)
                else:
                    if "prune" in model:
                        cnt_neg_pruned += 1
                    else:
                        cnt_neg_not_pruned += 1
                # check whether there is any single swap which would lead to smaller model (according to criterium)
                nodes = list(set(sum(map(lambda e: e[:2], model["edge"]), [])))
                for (a,b) in [(a, b) for idx, a in enumerate(nodes) for b in nodes[idx + 1:]]:
                    tmp = [('edge',sorted([list(map(lambda x: a if x==b else (b if x==a else x), e[:2])) + [e[2]] for e in model["edge"]]))]
                    if criterium=='min-adj':
                        tmp = adj(tmp)
                    if tmp < order:
                        if "prune" not in model:
                            print("WARNING:", order, (a,b), tmp)
                        if k == 0:
                            cnt_pos_prunable += 1
                        else:
                            cnt_neg_prunable += 1
                        break
            elif k in [0, selected_neg_example_k]:
                # only generate the positive and a single negative example per isomorphy class
                print(f"#{'pos' if k == 0 else 'neg'}({f'id_C{c}H{h}N{n}O{o}_{i}_iso{j}_'+ ('pos' if k == 0 else 'neg@100')},{{{','.join([pred + '(' + ','.join([str(v) for v in val]) + ')' for pred, vals in {'edge': model['edge']}.items() for val in vals])}}},{{}},{{molecular_formula(\"C\",{c}). molecular_formula(\"H\",{h}). molecular_formula(\"N\",{n}). molecular_formula(\"O\",{o}).}}).", file=sys.stderr)
    if CHECK:
        print(f"cnt_pos_pruned = {cnt_pos_pruned}, cnt_pos_not_pruned = {cnt_pos_not_pruned}, cnt_pos_prunable = {cnt_pos_prunable},\ncnt_neg_pruned = {cnt_neg_pruned}, cnt_neg_not_pruned = {cnt_neg_not_pruned}, cnt_neg_prunable = {cnt_neg_prunable}")

except KeyboardInterrupt:
    print()
