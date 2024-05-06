import sys
import re
import pprint

try:
    i = 0
    for line in sys.stdin:
        if "edge" in line:
            i += 1
            model = dict()
            for symbol in line.split():
                group = re.match("^([^(]+)\\(([^)]*)\\)$", symbol)
                pred = group[1]
                args = group[2].split(",")
                if pred not in model.keys():
                    model[pred] = []
                model[pred].append([arg[1:-1] if arg.startswith('"') else int(arg) if arg[0].isnumeric() else arg for arg in args])
            def _model(pred, default = []):
                return model[pred] if pred in model else default
            edges = list(filter(lambda x: len(x) == 3 and x[0] < x[1], model["edge"]))
            size = max([e[0] for e in edges] + [e[1] for e in edges])
            adj = [[0 for j in range(size)] for i in range(size)]
            for e in edges:
                adj[e[0]-1][e[1]-1] = e[2]
                adj[e[1]-1][e[0]-1] = e[2]
            for i, row in enumerate(adj):
                if i == 0:
                    print('\n\n    ' + '   '.join((str(j) for j in range(1, size+1))))
                    print("  " + "+---" * len(row) + "+")
                print(f"{i+1} | {' | '.join((str(c) for c in row))} |")
                print("  " + "+---" * len(row) + "+")
            pprint.pp(model)
        elif not line.startswith("Answer"):
            print(line, end="")

except KeyboardInterrupt:
    print()