# uage like:
# { genmol to-factbase -f C5H11NO2; cat naive-SBASS.lp | sed '1,20d'; } | python naive-SBASS.py

# Dynamically modify the following rules to account for all present element symbols in the absence of #sum aggregates:
# - `size(S)` in line 23
# - `type(I, E)` in lines 39-41
# - `VALENCE_SUM` in line 78


from sys import stdin
from re import match

try:
    elements = []
    atomic_numbers = dict()
    type_done = False
    for line in stdin:
        if g := match(r'^molecular_formula\("([A-Za-z]+)",\s*([0-9]+)\).$', line):
            if g.group(1) != 'H':
                elements.append(g.group(1))
            print(line, end='')
        if g := match(r'^element\("([A-Za-z]+)",\s*([0-9]+),\s*([0-9]+)\).$', line):
            if g.group(1) != 'H':
                atomic_numbers[g.group(1)] = int(g.group(2))
            print(line, end='')
        elif line.startswith('size'):
            print(f'size({'+'.join(elements)}) :- {', '.join([f'molecular_formula("{e}", {e})' for e in elements])}.')
        elif line.startswith('type'):
            if not type_done:
                type_done = True
                sorted_elements = sorted(elements, key=lambda e: atomic_numbers[e])
                for k in range(len(sorted_elements)):
                    print(f'type(I, "{sorted_elements[k]}") :- node(I){', ' if len(sorted_elements)!=1 else ''}{', '.join([f'molecular_formula("{e}", {e})' for e in sorted_elements[:k+(1 if k!=len(sorted_elements)-1 else 0)]])}{f', {'+'.join(sorted_elements[:k])} < I' if k!=0 else ''}{f', I <= {'+'.join(sorted_elements[:k+1])}' if k!=len(sorted_elements)-1 else ''}.')
        elif g := match(r'^(\s*)VALENCE_SUM', line):
            print(g.group(1) + f'VALENCE_SUM = {' + '.join([f'{e}*V{e}' for e in elements])}, {', '.join([f'molecular_formula("{e}", {e}), element("{e}", _, V{e})' for e in elements])},')
        else:
            print(line, end='')

except KeyboardInterrupt:
    print()


"""
$ { genmol to-factbase -f C4H9NO2; cat naive-SBASS.lp | sed '1,20d'; } | python naive-SBASS.py | gringo -o smodels | ./sbass | clasp 0 --project=show -q
clasp version 3.3.4
Reading from stdin
Solving...
SATISFIABLE

Models       : 69324
Calls        : 1
Time         : 7.783s (Solving: 7.78s 1st Model: 0.01s Unsat: 0.20s)
CPU Time     : 7.780s

$ clingo 0 naive.lp <(genmol to-factbase -f C4H9NO2) -q
clingo version 5.4.0
Reading from naive.lp ...
Solving...
SATISFIABLE

Models       : 69324
Calls        : 1
Time         : 2.390s (Solving: 2.38s 1st Model: 0.00s Unsat: 0.00s)
CPU Time     : 2.389s
"""
