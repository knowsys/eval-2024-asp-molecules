import sys
import re

if len(sys.argv) > 1 and sys.argv[1] == "-h":
    exit("Usage: cat chemdata.csv | python chemdata-sort.py > chemdata-sort.csv")

try:
    lines = []
    SUB = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
    for i, line in enumerate(sys.stdin):
        if i == 0:
            print(line, end="") # csv header
        else:
            sumformula = line.split(",")[0].translate(SUB)
            atom_count = sum((int(sub) if sub != "" else 1 for sub in [x.group(1) for x in re.finditer( r'(?:[A-GI-Z][a-z]*|H[a-z]+)(\d*)', sumformula)]))
            #print(sumformula, atom_count, file=sys.stderr)
            lines.append((atom_count, line))
    lines.sort(key=lambda x: x[0])
    for atom_count, line in lines:
        print(line, end="")
except KeyboardInterrupt:
    print()
