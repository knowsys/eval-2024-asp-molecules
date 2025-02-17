from subprocess import run, TimeoutExpired, DEVNULL, call
from typing import List, Dict, Union, Tuple

from dataclasses import dataclass, field
from dataclasses_json import dataclass_json
from os.path import isfile, exists
from os import makedirs, environ
from sys import stderr

import argparse

from pprint import PrettyPrinter
pp = PrettyPrinter(indent=2)

import matplotlib.pyplot as plt

import re

path = "sb-results.json"
diagram_path = "diagrams"

PROG_NAIVE = environ.get('PROG_NAIVE', "alternative/naive.lp")
PROG_CANONICAL = environ.get('PROG_CANONICAL', "alternative/symmetry.lp")
PROG_SMILES = environ.get('PROG_SMILES', "smiles_min.lp")
PROG_NAIVE_SBASS = environ.get('PROG_NAIVE_SBASS', "naive-SBASS.lp")

CHEMDATA = environ.get('CHEMDATA', "chemdata-sort.csv")
GENMOL = environ.get('GENMOL', "LD_LIBRARY_PATH=/nix/store/a3zlvnswi1p8cg7i9w4lpnvaankc7dxx-gcc-12.3.0-lib/lib/:/nix/store/2n7vw6bfc8lbfhzlcahdzacp23si4rxc-openssl-1.1.1q/lib/:$LD_LIBRARY_PATH ../genmol")

if not isfile(PROG_SMILES):
    call(['bash', './prepare-asp-programs.sh'])

num_threads = 2
timeout = 60 # in seconds

record_count = 2500 # total number of datapoints to collect
min_num_models = 0

_DEBUG: bool = False


@dataclass_json
@dataclass
class NOMResultSet:
    formula: List[str] = field(default_factory=lambda: [])
    molgen_num_models: List[Union[int, None]] = field(default_factory=lambda: [])
    smiles_num_models: List[Union[int, None]] = field(default_factory=lambda: [])
    canonical_num_models: List[Union[int, None]] = field(default_factory=lambda: [])
    sbass_num_models: List[Union[int, None]] = field(default_factory=lambda: [])
    breakid_num_models: List[Union[int, None]] = field(default_factory=lambda: [])
    naive_num_models: List[Union[int, None]] = field(default_factory=lambda: [])


def to_num(txt: List[str], idx: int, f) -> Union[int, float, None]:
    try:
        return f(txt[idx])
    except (ValueError, IndexError) as error:
        return None


BLUE = '\33[34m'
VIOLET = '\33[35m'
C_END = '\33[0m'


def print_cmd(cmd: str, numbers1: Union[List[str], None] = None, numbers2: Union[List[str], None] = None):
    global _DEBUG

    if _DEBUG:
        print(BLUE + ' '.join([tok if ' ' not in tok else f'"{tok}"' for tok in [' '.join(c.replace('\n',"\"$'\\n'\" ").split()) for c in cmd]]) + C_END, file=stderr)
        if numbers1 is not None:
            print(f"{VIOLET}numbers1 = {numbers1}{C_END}", file=stderr)
        if numbers2 is not None:
            print(f"{VIOLET}numbers2 = {numbers2}{C_END}", file=stderr)


def check_formula_smiles_support(formula: str) -> bool:
    if 0 != run(["bash", "-c", f"{GENMOL} to-factbase -f {formula} || exit 1"], stdout=DEVNULL, stderr=DEVNULL).returncode:
        # have invalid sumformula or unsupported elements
        print("Reject", file=stderr)
        return False
    print("Pass", file=stderr)
    return True


def measure_num_models_smiles(formula: str) -> Union[int, None]:
    try:
        numbers = run((cmd := ["bash", "-c", f"clingo 0 --quiet=2,0,2 -t {num_threads} {PROG_SMILES} <({GENMOL} to-factbase -f {formula}) \
                    | grep -oP ':.*|^\d+$' \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout).stdout.splitlines()
        print_cmd(cmd, numbers)

        num_models = to_num(numbers, 0, int)
        return num_models
    except TimeoutExpired:
        return None


def measure_num_models_molgen(formula: str) -> Union[int, None]:
    try:
        numbers = run((cmd := ["bash", "-O", "expand_aliases", "-c", f"[ -f .bash_aliases ] && source .bash_aliases\n \
                    molgen {formula} -v 2>&1 \
                    | tail -1 \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout).stdout.splitlines()
        print_cmd(cmd, numbers)

        num_models = to_num(numbers, 0, int)
        return num_models

    except TimeoutExpired:
        return None


def measure_num_models_naive(formula: str) -> Union[int, None]:
    try:
        numbers = run((cmd := ["bash", "-c", f"clingo 0 --quiet=2,0,2 -t {num_threads} {PROG_NAIVE} <({GENMOL} to-factbase -f {formula}) \
                    | grep -oP ':.*|^\d+$' \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout).stdout.splitlines()
        print_cmd(cmd, numbers)

        num_models = to_num(numbers, 0, int)
        return num_models
    except TimeoutExpired:
        return None


def measure_num_models_canonical(formula: str) -> Union[int, None]:
    try:
        numbers = run((cmd := ["bash", "-c", f"clingo 0 --quiet=2,0,2 -t {num_threads} {PROG_NAIVE} {PROG_CANONICAL} <({GENMOL} to-factbase -f {formula}) \
                    | grep -oP ':.*|^\d+$' \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout).stdout.splitlines()
        print_cmd(cmd, numbers)

        num_models = to_num(numbers, 0, int)
        return num_models
    except TimeoutExpired:
        return None


def measure_num_models_sbass(formula: str) -> Union[int, None]:
    # With the (static) aggregate-free encoding, suitable only for CHNO-formulae:
    """
    m = re.search("^(?:C(\d*))?(?:H(\d*))?(?:N(\d*))?(?:O(\d*))?$", formula)
    if m is None:
        return -1
    c = int(m.group(1)) if m.group(1) not in ['', None] else 1
    h = int(m.group(2)) if m.group(2) not in ['', None] else 1
    n = int(m.group(3)) if m.group(3) not in ['', None] else 1
    o = int(m.group(4)) if m.group(4) not in ['', None] else 1
    try:
        output = run((cmd := ["bash", "-O", "expand_aliases", "-c", f"[ -f .bash_aliases ] && source .bash_aliases\n \
                    gringo {PROG_NAIVE_SBASS} --const c={c} --const h={h} --const n={n} --const o={o} -o smodels \
                    | sbass \
                    | clasp 0 --project=show --quiet=2,0,2 -t {num_threads} \
                    | grep -oP ':.*|^\d+$' \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout)
    """
    # With the dynamically generated aggregate-free encoding, not restricting element symbols (for comparability):
    try:
        output = run((cmd := ["bash", "-O", "expand_aliases", "-c", f"[ -f .bash_aliases ] && source .bash_aliases\n \
                    {{ {GENMOL} to-factbase -f {formula}; cat {PROG_NAIVE_SBASS} | sed '1,20d'; }} | python naive-SBASS.py \
                    | gringo -o smodels \
                    | sbass \
                    | clasp 0 --project=show --quiet=2,0,2 -t {num_threads} \
                    | grep -oP ':.*|^\d+$' \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout)
        numbers = output.stdout.splitlines()
        print_cmd(cmd, numbers)

        num_models = to_num(numbers, 0, int)
        return num_models
    except TimeoutExpired:
        return None


def measure_num_models_breakid(formula: str) -> Union[int, None]:
    try:
        output = run((cmd := ["bash", "-O", "expand_aliases", "-c", f"[ -f .bash_aliases ] && source .bash_aliases\n \
                    gringo {PROG_NAIVE} <({GENMOL} to-factbase -f {formula}) -o smodels \
                    | breakID -asp \
                    | tail -n +2 \
                    | cat - <(echo '0') \
                    | clasp 0 --project=show --quiet=2,0,2 -t {num_threads} \
                    | grep -oP ':.*|^\d+$' \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout)
        numbers = output.stdout.splitlines()
        print_cmd(cmd, numbers)

        num_models = to_num(numbers, 0, int)
        return num_models
    except TimeoutExpired:
        return None


def evaluate_num_models(records: int, result_set: NOMResultSet = NOMResultSet()) -> NOMResultSet:
    with open(CHEMDATA, "r") as chemdata:
        SUB = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
        recovery = result_set.formula[-1] if len(result_set.formula) > 0 else None
        i = 0
        if recovery is not None:
            i = len(list(filter(lambda x: x[0] > min_num_models and x[1] > min_num_models, zip(result_set.molgen_num_models, result_set.smiles_num_models))))
        for line in chemdata.readlines():
            if i == 0:
                i += 1
                continue
            if i >= records+1:
                break
            sumformula = line.split(",")[0].translate(SUB)
            if recovery is not None:
                if sumformula == recovery:
                    recovery = None
                continue
            print(f'{i}: {sumformula} ? ', end="", file=stderr)
            if check_formula_smiles_support(sumformula):
                molgen_num_models = measure_num_models_molgen(sumformula)
                if molgen_num_models is not None:
                    smiles_num_models = measure_num_models_smiles(sumformula)
                    if smiles_num_models is not None:
                        canonical_num_models = measure_num_models_canonical(sumformula)
                        breakid_num_models = measure_num_models_breakid(sumformula)
                        sbass_num_models = measure_num_models_sbass(sumformula)
                        naive_num_models = measure_num_models_naive(sumformula)
                        result_set.formula.append(sumformula)
                        result_set.molgen_num_models.append(molgen_num_models)
                        result_set.smiles_num_models.append(smiles_num_models)
                        result_set.canonical_num_models.append(canonical_num_models)
                        result_set.sbass_num_models.append(sbass_num_models)
                        result_set.breakid_num_models.append(breakid_num_models)
                        result_set.naive_num_models.append(naive_num_models)
                        print(f"{i},{sumformula},{molgen_num_models},{smiles_num_models},{canonical_num_models},{sbass_num_models},{breakid_num_models},{naive_num_models}")
                        if molgen_num_models > min_num_models and smiles_num_models > min_num_models:
                            i += 1

    return result_set


def diagram(filename: str, label: str, data: List[Tuple[int, int, int, int, int]], legend: Tuple[str, str, str, str, str]):
    plt.rcParams["figure.figsize"] = [8.00, 3.50]
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams['font.size'] = 12
    prop_cycle_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax = plt.subplots()

    plt.scatter(range(len(data)), list(map(lambda x: x[2], data)), s=5, color=prop_cycle_colors[1], marker='s')
    plt.scatter(range(len(data)), list(map(lambda x: x[1], data)), s=5, color=prop_cycle_colors[0], marker='o')
    sbass_dat = list(filter(lambda x: x[1] != -1, zip(range(len(data)), list(map(lambda x: x[4], data)))))
    plt.scatter(list(map(lambda x: x[0], sbass_dat)), list(map(lambda x: x[1], sbass_dat)), s=5, color=prop_cycle_colors[6], marker='<') # SBASS can only deal with CHNO formulae (-1 is placeholder for unsupported formula)
    plt.scatter(range(len(data)), list(map(lambda x: x[3], data)), s=5, color=prop_cycle_colors[4], marker='v')
    plt.scatter(range(len(data)), list(map(lambda x: x[5], data)), s=5, color=prop_cycle_colors[2], marker='_')
    plt.scatter(range(len(data)), list(map(lambda x: x[0], data)), s=5, color=prop_cycle_colors[3], marker='_')

    plt.scatter([], [], s=70, label=legend[0], color=prop_cycle_colors[3], marker='_')
    plt.scatter([], [], s=70, label=legend[1], color=prop_cycle_colors[0], marker='o')
    plt.scatter([], [], s=70, label=legend[2], color=prop_cycle_colors[1], marker='s')
    plt.scatter([], [], s=70, label=legend[3], color=prop_cycle_colors[4], marker='v')
    plt.scatter([], [], s=70, label=legend[4], color=prop_cycle_colors[6], marker='<')
    plt.scatter([], [], s=70, label=legend[5], color=prop_cycle_colors[2], marker='_')

    ax.set_yscale('log')
    ax.set_ylim(auto=True)
    ax.set_ylabel(label)
    ax.get_xaxis().set_visible(False)
    plt.yticks(fontsize=10)
    ax.legend(loc='upper left', ncols=len(legend)/2, labelspacing=0.5)

    # save the figure in PDF format and close it
    if not exists(diagram_path):
        makedirs(diagram_path)

    plt.savefig(f"{diagram_path}/diagram_{filename}.pdf")
    plt.savefig(f"{diagram_path}/diagram_{filename}.svg")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare number of models of ASP-programs for chemical structure exploration against Molgen.',
                                     epilog='Make sure to pipe stdout into a file, so that recovery is possible (e.g.: python eval-sb.py -d -t sb-data.csv 1>> sb-data.csv).')
    parser.add_argument('--debug', '-d', action='store_true',
                        help='show debug output with commands and resulting numbers (default: no)')
    parser.add_argument('--print-results', '-p', action='store_true',
                        help='pretty-print the results to the console (default: no)')
    parser.add_argument('--render-diagrams', '-r', action='store_true',
                        help='create diagrams from the results (default: no)')
    parser.add_argument('--try-recover', '-t', metavar='FILE',
                        help='Try to recover from output of previous run')
    eval_group = parser.add_argument_group('Evaluation parameters')
    eval_group.add_argument('--num_threads', default=num_threads, type=int,
                            help=f'Number of solver threads to use via Clingo option, (default: -t {num_threads})')
    eval_group.add_argument('--timeout', default=timeout, type=int,
                            help=f'Timeout for Clingo / Molgen invocations (default: {timeout}sec)')
    eval_group.add_argument('--record_count', default=record_count, type=int,
                            help=f'Number of records to collect (default: {record_count})')
    eval_group.add_argument('--min_num_models', default=min_num_models, type=int,
                            help=f'Minimum number of models / structures required for a measured record to count and be displayed in the diagram (default: {min_num_models})')
    args = parser.parse_args()
    if args.debug:
        _DEBUG = True
    num_threads = args.num_threads
    timeout = args.timeout
    record_count = args.record_count
    min_num_models = args.min_num_models

    if isfile(path):
        with open(path, "r") as fp:
            result_set = NOMResultSet.from_json(fp.read())
    else:
        if args.try_recover is not None:
            with open(args.try_recover, "r") as fp:
                print("Recovery...", file=stderr)
                result_set = NOMResultSet()
                i = 1
                canonical_done = []
                sbass_done = []
                breakid_done = []
                naive_done = []
                for line in fp.readlines():
                    parts = line.split(",")
                    if len(parts) >= 4:
                        if int(parts[0]) != i:
                            print(f"Warning: record {i} is missing...", file=stderr)
                        sumformula = parts[1]
                        result_set.formula.append(sumformula)
                        molgen = int(parts[2])
                        smiles = int(parts[3])
                        canonical = int(parts[4]) if len(parts) > 4 else None
                        breakid = int(parts[5]) if len(parts) > 5 else None
                        naive = int(parts[6]) if len(parts) > 6 else None
                        sbass = int(parts[7]) if len(parts) > 7 else None
                        result_set.molgen_num_models.append(molgen)
                        result_set.smiles_num_models.append(smiles)
                        result_set.canonical_num_models.append(canonical)
                        result_set.sbass_num_models.append(sbass)
                        result_set.breakid_num_models.append(breakid)
                        result_set.naive_num_models.append(naive)
                        if molgen > min_num_models and smiles > min_num_models:
                            i += 1
                    elif len(parts) == 3 and parts[0] == "canonical":
                        canonical_done.append(sumformula)
                        sumformula = parts[1]
                        canonical = int(parts[2]) if not parts[2].startswith("None") else None
                        for idx, s in enumerate(result_set.formula):
                            if s == sumformula:
                                result_set.canonical_num_models[idx] = canonical
                    elif len(parts) == 3 and parts[0] == "sbass":
                        sbass_done.append(sumformula)
                        sumformula = parts[1]
                        sbass = int(parts[2]) if not parts[2].startswith("None") else None
                        for idx, s in enumerate(result_set.formula):
                            if s == sumformula:
                                result_set.sbass_num_models[idx] = sbass
                    elif len(parts) == 3 and parts[0] == "breakid":
                        breakid_done.append(sumformula)
                        sumformula = parts[1]
                        breakid = int(parts[2]) if not parts[2].startswith("None") else None
                        for idx, s in enumerate(result_set.formula):
                            if s == sumformula:
                                result_set.breakid_num_models[idx] = breakid
                    elif len(parts) == 3 and parts[0] == "naive":
                        naive_done.append(sumformula)
                        sumformula = parts[1]
                        naive = int(parts[2]) if not parts[2].startswith("None") else None
                        for idx, s in enumerate(result_set.formula):
                            if s == sumformula:
                                result_set.naive_num_models[idx] = naive
                i = 1
                missing = []
                for j in range(len(result_set.formula)):
                    sumformula = result_set.formula[j]
                    molgen = result_set.molgen_num_models[j]
                    smiles = result_set.smiles_num_models[j]
                    canonical = result_set.canonical_num_models[j]
                    sbass = result_set.sbass_num_models[j]
                    breakid = result_set.breakid_num_models[j]
                    if molgen > min_num_models and smiles > min_num_models:
                        if ((canonical is None and sumformula not in canonical_done) or (sbass is None and sumformula not in sbass_done) or (breakid is None and sumformula not in breakid_done) or (naive is None and sumformula not in naive_done)) and i < record_count+1:
                            missing.append((molgen,smiles,canonical,sbass,breakid,j,sumformula))
                        i += 1
                for molgen, smiles, canonical, sbass, breakid, j, sumformula in sorted(missing):
                    if canonical is None and sumformula not in canonical_done:
                        canonical = measure_num_models_canonical(sumformula)
                        print(f"canonical,{sumformula},{canonical}")
                        for idx, s in enumerate(result_set.formula):
                            if s == sumformula:
                                result_set.canonical_num_models[idx] = canonical
                    if sbass is None and sumformula not in sbass_done:
                        sbass = measure_num_models_sbass(sumformula)
                        print(f"sbass,{sumformula},{sbass}")
                        for idx, s in enumerate(result_set.formula):
                            if s == sumformula:
                                result_set.sbass_num_models[idx] = sbass
                    if breakid is None and sumformula not in breakid_done:
                        breakid = measure_num_models_breakid(sumformula)
                        print(f"breakid,{sumformula},{breakid}")
                        for idx, s in enumerate(result_set.formula):
                            if s == sumformula:
                                result_set.breakid_num_models[idx] = breakid
                    if naive is None and sumformula not in naive_done:
                        naive = measure_num_models_naive(sumformula)
                        naive_done.append(sumformula)
                        print(f"naive,{sumformula},{naive}")
                        for idx, s in enumerate(result_set.formula):
                            if s == sumformula:
                                result_set.naive_num_models[idx] = naive
                j = 0
                i = 1
                for molgen, smiles, sumformula in zip(result_set.molgen_num_models, result_set.smiles_num_models, result_set.formula):
                    j += 1
                    if molgen > min_num_models and smiles > min_num_models:
                        i += 1
                        if i > record_count:
                            break
                del result_set.formula[j:]
                del result_set.molgen_num_models[j:]
                del result_set.smiles_num_models[j:]
                del result_set.canonical_num_models[j:]
                del result_set.breakid_num_models[j:]
                del result_set.naive_num_models[j:]
                #print(pp.pformat(result_set), file=stderr)
            result_set = evaluate_num_models(records=record_count, result_set=result_set)
        else:
            result_set = evaluate_num_models(records=record_count)
        with open(path, "w") as fp:
            fp.write(result_set.to_json())

    if args.print_results:
        print(pp.pformat(sorted(zip(result_set.formula, result_set.molgen_num_models, result_set.smiles_num_models), key=lambda x: x[1]-x[2])), file=stderr)

    if args.render_diagrams:
        # lexicographic sort --> 1st sort by molgen, 2nd sort by smiles
        data = sorted(zip(result_set.molgen_num_models, result_set.smiles_num_models, result_set.canonical_num_models, result_set.breakid_num_models, result_set.sbass_num_models, result_set.naive_num_models), \
            key=lambda x: (x[0], x[1], x[2] if x[2] is not None else 0, x[3] if x[3] is not None else 0, x[4] if x[4] is not None else 0, x[5] if x[5] is not None else 0))
        # check how many are equal at the start
        equal_up_to = 0
        last_equal_up_to = 0
        equal_count = 0
        for molgen, smiles, canonical, sbass, breakid, naive in data:
            if molgen != smiles:
                break
            equal_count += 1
            if last_equal_up_to < molgen:
                equal_up_to = last_equal_up_to
            last_equal_up_to = molgen
        print(f"Equal up to {equal_up_to} ({equal_count} records) !!!", file=stderr)
        print(f"Have {len(data)} data points...", file=stderr)
        data = [(m, s, c, b, sb, n) for (m, s, c, b, sb, n) in data if m > min_num_models and s > min_num_models]
        print(f"Filter for records with at least {min_num_models} models... Have {len(data)} data points...", file=stderr)
        print(data, file=stderr)
        diagram("number_of_models-comparison", "Number of models", data, ("Molgen", "Our encoding", "Canonical", "BreakID", "SBASS", "Naive"))

if False: #__name__ == "__main__":
    f_symmetry_breaking = open(f"results/symmetry_breaking.csv", "r")
    lines_symmetry_breaking = [l.strip().split(',') for l in f_symmetry_breaking.readlines()]
    lines_symmetry_breaking[0].append("sbass")

    for line in lines_symmetry_breaking[1:]:

        _DEBUG = True

        result = measure_num_models_sbass(line[0])
        if result is not None:
            print(f"{line[0]};{result}")
            line.append(result)
        else:
            print(f"{line[0]};None")
            line.append(None)

    open("results/new_symmetry_breaking.csv", "w").write('\n'.join([','.join([str(v) for v in l]) for l in lines_symmetry_breaking]))

if False: #__name__ == "__main__":
    f_symmetry_breaking = open(f"results/symmetry_breaking.csv", "r")
    lines_symmetry_breaking = [l.strip().split(',') for l in f_symmetry_breaking.readlines()]

    result_set = NOMResultSet()
    for line in lines_symmetry_breaking[1:]:
        for i in range(1, len(line)):
            v = line[i]
            line[i] = int(v) if v != 'None' else None
        # sumformula,molgen,smiles,canonical,breakid,naive,sbass
        result_set.formula.append(line[0])
        result_set.molgen_num_models.append(line[1])
        result_set.smiles_num_models.append(line[2])
        result_set.canonical_num_models.append(line[3])
        result_set.sbass_num_models.append(line[6])
        result_set.breakid_num_models.append(line[4])
        result_set.naive_num_models.append(line[5])

    # lexicographic sort --> 1st sort by molgen, 2nd sort by smiles
    data = sorted(zip(result_set.molgen_num_models, result_set.smiles_num_models, result_set.canonical_num_models, result_set.breakid_num_models, result_set.sbass_num_models, result_set.naive_num_models), \
        key=lambda x: (x[0], x[1], x[2] if x[2] is not None else 0, x[3] if x[3] is not None else 0, x[4] if x[4] is not None else 0, x[5] if x[5] is not None else 0))
    # check how many are equal at the start
    equal_up_to = 0
    last_equal_up_to = 0
    equal_count = 0
    for molgen, smiles, canonical, sbass, breakid, naive in data:
        if molgen != smiles:
            break
        equal_count += 1
        if last_equal_up_to < molgen:
            equal_up_to = last_equal_up_to
        last_equal_up_to = molgen
    print(f"Equal up to {equal_up_to} ({equal_count} records) !!!", file=stderr)
    print(f"Have {len(data)} data points...", file=stderr)
    data = [(m, s, c, b, sb, n) for (m, s, c, b, sb, n) in data if m > min_num_models and s > min_num_models]
    print(f"Filter for records with at least {min_num_models} models... Have {len(data)} data points...", file=stderr)
    print(data, file=stderr)
    diagram("new_number_of_models-comparison", "Number of models", data, ("Molgen", "Our encoding", "Canonical", "BreakID", "SBASS", "Naive"))

if False: #__name__ == "__main__":
    f_symmetry_breaking = open(f"results/symmetry_breaking.csv", "r")
    lines_symmetry_breaking = [l.strip().split(',') for l in f_symmetry_breaking.readlines()]

    #factors_genmol = [0 for _ in range(10000)]
    factors_sbass = [0 for _ in range(10000)]
    for line in lines_symmetry_breaking[1:]:
        for i in range(1, len(line)):
            v = line[i]
            line[i] = int(v) if v != 'None' else None

        #if line[2] is not None:
        #    factor = line[2] // line[1]
        #    for f in range(factor, 10000):
        #        factors_genmol[f] += 1

        if line[6] is not None and line[6] != -1:
            factor = line[6] // line[1]
            for f in range(factor, 10000):
                factors_sbass[f] += 1

    #print(factors_genmol)
    #print(factors_sbass)

    #ratios_genmol = []
    ratios_sbass = []
    for i in range(1, 10000):
        #ratios_genmol.append(100.0 * factors_genmol[i] / len(lines_symmetry_breaking[1:]))
        ratios_sbass.append(100.0 * factors_sbass[i] / len(lines_symmetry_breaking[1:]))

    from math import log10

    #print(ratios_genmol)

    # calculate x scale from breakid line in the SVG (line 843)
    #breakid = "m 589.15956,157.42023 h 2.55118 2.35276 2.15433 l 3.85512,-3.34488 4.13858,-2.12599 1.50236,-0.51023 2.0693,-1.89922 0.0283,-8.10709 10.91339,-6.94488 7.73858,-8.10709 6.00945,-6.51967 9.04252,-13.2378 9.60945,-5.89606 7.05827,-6.321265 11.5937,-10.28977 15.05197,-14.8252 9.60945,-11.67874 18.65197,-6.29291 18.65198,-6.63307 24.68976,-3.99685 18.62363,-1.70079 18.65197,-1.04882 24.66142,-0.62362 h 10.91338 l 7.73859,-0.22677 24.66142,-0.70867"
    #l = breakid.split()
    #start_pos = [float(v) for v in l[1].split(',')]
    #l = l[2:]
    #end_pos = [v for v in start_pos]
    #for v in l:
    #    if v == 'h':
    #        h_mode = True
    #    elif v == 'l':
    #        h_mode = False
    #    elif h_mode:
    #        end_pos[0] += float(v)
    #    else:
    #        dxdy = v.split(',')
    #        end_pos[0] += float(dxdy[0])
    #        end_pos[1] += float(dxdy[1])
    #print(end_pos[0] - start_pos[0])
    x_scale = 272.52283 / log10(10000)

    # the following line defines the outer border of the diagram in the SVG (line 128) --> both offsets and y_scale
    # 734.7603,189.13992 H 589.15956 V 5.029665 H 880.3344 V 189.13992 Z
    x_offset = 589.15956
    y_offset = 189.13992
    y_scale = - (189.13992 / 100.0)

    def to_svg(ratios):
        last_x = 0.0
        last_y = 0.0
        is_first = True
        for x,y in map(lambda v: (0.0 if v[0]==0 else log10(v[0]), v[1]), enumerate(ratios)):
            #print(last_x, x)
            #print(last_y, y)
            if is_first or last_x + log10(1000.0)/30 <= x:
                #print(f"{x},{y} ", end="")
                print(f"{'m' if is_first else ''} {x_offset + x_scale * x if is_first else x_scale * (x-last_x):.5f},{y_offset + y_scale * y if is_first else y_scale * (y-last_y):.5f}", end='')
                last_x = x
                last_y = y
                is_first = False
        print()

    #to_svg(ratios_genmol)
    to_svg(ratios_sbass)
    print("y position of end of line:", y_offset + y_scale * ratios_sbass[-1] - (50.762137+8.944272/2))
