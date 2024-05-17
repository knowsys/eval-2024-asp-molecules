from sys import stderr

with open("validation.csv", "r") as f:
    num_timeout = 0
    num_bad_molecule_formula = 0
    num_skipped = 0
    num_unconnected = 0
    num_bad_smiles = 0
    num_non_standard_valences = 0
    num_unsatisfiable = 0
    num_grounding_failed = 0
    num_satisfiable = 0
    total_num_models = 0
    for l in f.readlines():
        parts = l.split(",")
        val = int(parts[0])
        idx = int(parts[1])
        if val < 0:
            if val == -420:
                num_timeout += 1
            else:
                print(f"old timeout value in: {l}", end="", file=stderr)
        elif val == 1:
            num_bad_molecule_formula += 1
        elif val == 2:
            num_skipped += 1
        elif val == 3:
            num_unconnected += 1
        elif val == 4 or val == 5 or val == 6:
            num_bad_smiles += 1
        elif val == 7:
            num_non_standard_valences += 1
        elif val == 8:
            num_unsatisfiable += 1
        elif val == 9:
            num_grounding_failed += 1
        elif val > 100:
            num_satisfiable += 1
            total_num_models += val - 100
        else:
            print(f"unkown value in: {l}", end="", file=stderr)
    print('num_timeout =', num_timeout)
    print('num_bad_molecule_formula =', num_bad_molecule_formula)
    print('num_skipped =', num_skipped)
    print('num_unconnected =', num_unconnected)
    print('num_bad_smiles =', num_bad_smiles)
    print('num_non_standard_valences =', num_non_standard_valences)
    print('num_unsatisfiable =', num_unsatisfiable)
    print('num_grounding_failed =', num_grounding_failed)
    print('num_satisfiable =', num_satisfiable)
    print('total_num_models =', total_num_models)
    print('avg_num_isomorphic_models =', float(total_num_models) / num_satisfiable)
    print()
    num_all = num_timeout + num_bad_molecule_formula + num_skipped + num_unconnected + num_bad_smiles + num_non_standard_valences + num_unsatisfiable + num_grounding_failed + num_satisfiable
    print("all: ", num_all)
    data_set_size = num_all - num_bad_molecule_formula - num_unconnected - num_non_standard_valences
    print("data set size: ", data_set_size)
    print("considered: ", data_set_size - num_bad_smiles)
