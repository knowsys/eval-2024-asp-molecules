from subprocess import run, TimeoutExpired, call
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
import numpy as np
from statistics import median

PROG_NAIVE = environ.get('PROG_NAIVE', "alternative/naive.lp")
PROG_CANONICAL = environ.get('PROG_CANONICAL', "alternative/symmetry.lp")
PROG_SMILES = environ.get('PROG_SMILES_EVAL', "smiles_eval.lp")
PROG_NAIVE_SBASS = environ.get('PROG_NAIVE_SBASS', "naive-SBASS.lp")

if not isfile(PROG_SMILES):
    call(['bash', './prepare-asp-programs.sh'])

num_threads = 2
timeout = 60 # in seconds

series = [(0,0),(0,1),(1,1),(2,0),(2,2)] # (cycles, oxygens)
carbons_limit = 30

repetition_count = 5

path = "results.json"
diagram_path = "diagrams"

_DEBUG: bool = False

GREEN = '\033[32m'
BLUE = '\33[34m'
VIOLET = '\33[35m'
C_END = '\33[0m'

@dataclass_json
@dataclass
class Result:
    ground_prog_size: Union[int, None] = None
    num_models: Union[int, None] = None
    runtime: Union[float, None] = None
    solving_time: Union[float, None] = None
    first_model_time: Union[float, None] = None
    unsat_time: Union[float, None] = None
    grounding_time: Union[float, None] = None
    total_runtime: Union[float, None] = None
    per_model_runtime: Union[float, None] = field(init=False)

    breakid_runtime: Union[float, None] = None
    breakid_ground_prog_size: Union[int, None] = None
    breakid_num_symmetry_generators: Union[int, None] = None
    breakid_auxiliary_variables: Union[int, None] = None
    breakid_symmetry_breaking_clauses: Union[int, None] = None

    def __post_init__(self):
        if self.num_models is not None and self.total_runtime is not None:
            if self.num_models > 0:
                self.per_model_runtime = self.total_runtime * 1000.0 / self.num_models
            else:
                self.per_model_runtime = 0
        else:
            self.per_model_runtime = None


def to_num(txt: List[str], idx: int, f) -> Union[int, float, None]:
    try:
        return f(txt[idx])
    except (ValueError, IndexError) as error:
        return None


def print_cmd(cmd: str, numbers1: Union[List[str], None] = None, numbers2: Union[List[str], None] = None):
    global _DEBUG

    if _DEBUG:
        print(BLUE + ' '.join([tok if ' ' not in tok else f'"{tok}"' for tok in [' '.join(c.replace('\n',"\"$'\\n'\" ").split()) for c in cmd]]) + C_END, file=stderr)
        if numbers1 is not None:
            print(f"{VIOLET}numbers1 = {numbers1}{C_END}", file=stderr)
        if numbers2 is not None:
            print(f"{VIOLET}numbers2 = {numbers2}{C_END}", file=stderr)


def get_params_str(carbons: int, hydrogens: int, oxygens: int, nitrogens: int) -> str:
    from math import ceil, log
    atom_number = carbons + oxygens + nitrogens
    min_main_chain_len = min(2*ceil(log((atom_number-1)/2+1,3))+1, 2*ceil(log(atom_number+1,3)))
    return f"--const c={carbons} --const h={hydrogens} --const o={oxygens} --const n={nitrogens} --const m={min_main_chain_len}"


def measure_gringo(progs: str, params_str: str) -> Union[int, None]:
    try:
        numbers = run((cmd := ["bash", "-c", f"gringo {progs} {params_str} \
                    | wc -l"]), capture_output=True, text=True, timeout=timeout).stdout.splitlines()
        print_cmd(cmd, numbers)

        ground_prog_size = to_num(numbers, 0, int)
        return ground_prog_size
    except TimeoutExpired:
        return None


def measure_clingo(progs: str, params_str: str) -> Union[Result, None]:
    try:
        numbers = run((cmd := ["bash", "-c", f"clingo 0 --quiet=2,0,2 -t {num_threads} {progs} {params_str} \
                    | grep -oP ':.*|^\d+$' \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout).stdout.splitlines()
        print_cmd(cmd, numbers)

        num_models = to_num(numbers, 0, int)
        runtime = to_num(numbers, 2, float)
        solving_time = to_num(numbers, 3, float)
        first_model_time = to_num(numbers, 5, float)
        unsat_time = to_num(numbers, 6, float)
        total_cpu_time = to_num(numbers, 7, float)
        grounding_time = total_cpu_time - solving_time - first_model_time - unsat_time

        return Result(None, num_models, runtime, solving_time, first_model_time, unsat_time, grounding_time, total_cpu_time)
    except TimeoutExpired:
        return None


measure_naive = lambda c, h, o, n: measure_clingo(PROG_NAIVE, get_params_str(c, h, o, n))
measure_canonical = lambda c, h, o, n: measure_clingo(f"{PROG_NAIVE} {PROG_CANONICAL}", get_params_str(c, h, o, n))
measure_smiles = lambda c, h, o, n: measure_clingo(f"{PROG_SMILES}", get_params_str(c, h, o, n))

measure_ground_prog_size_naive = lambda c, h, o, n: measure_gringo(PROG_NAIVE, get_params_str(c, h, o, n))
measure_ground_prog_size_canonical = lambda c, h, o, n: measure_gringo(f"{PROG_NAIVE} {PROG_CANONICAL}", get_params_str(c, h, o, n))
measure_ground_prog_size_smiles = lambda c, h, o, n: measure_gringo(f"{PROG_SMILES}", get_params_str(c, h, o, n))


def measure_sbass(carbons: int, hydrogens: int, oxygens: int, nitrogens: int) -> Union[Result, None]:
    params_str = get_params_str(carbons, hydrogens, oxygens, nitrogens)

    try:
        output = run((cmd := ["bash", "-O", "expand_aliases", "-c", f"[ -f .bash_aliases ] && source .bash_aliases\n \
                    gringo {PROG_NAIVE_SBASS} {params_str} -o smodels \
                    | tee 1> >(wc -l) \
                        >(sbass --stats 2> \
                                >(grep -oP '=.*|^[\d.]+$' \
                                    | grep -oP '[0-9]+(\.[0-9]*)?' >&2) \
                                | tee 1> >(wc -l) \
                                    >(clasp 0 --project=show --quiet=2,0,2 -t {num_threads}) \
                        ) 2>&1 \
                    | grep -oP ':.*|^\d+$' 2>&1 \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout)
        numbers1 = output.stdout.splitlines()
        numbers2 = output.stderr.splitlines()
        print_cmd(cmd, numbers1, numbers2)

        ground_prog_size = to_num(numbers1, 0, int)
        num_models = to_num(numbers1, 2, int)
        runtime = to_num(numbers1, 4, float)
        solving_time = to_num(numbers1, 5, float)
        first_model_time = to_num(numbers1, 7, float)
        unsat_time = to_num(numbers1, 8, float)
        total_cpu_time = to_num(numbers1, 9, float)
        grounding_time = total_cpu_time - solving_time - first_model_time - unsat_time

        breakid_runtime = to_num(numbers2, -1, float)
        breakid_ground_prog_size = to_num(numbers1, 1, int) - ground_prog_size
        breakid_num_symmetry_generators = to_num(numbers2, 2, int)
        breakid_symmetry_breaking_clauses = to_num(numbers2, -3, int)

        return Result(ground_prog_size, num_models, runtime, solving_time, first_model_time, unsat_time, grounding_time, total_cpu_time,
            breakid_runtime, breakid_ground_prog_size, breakid_num_symmetry_generators, None, breakid_symmetry_breaking_clauses)

    except TimeoutExpired:
        return None


def measure_breakid(carbons: int, hydrogens: int, oxygens: int, nitrogens: int) -> Union[Result, None]:
    params_str = get_params_str(carbons, hydrogens, oxygens, nitrogens)

    try:
        output = run((cmd := ["bash", "-O", "expand_aliases", "-c", f"[ -f .bash_aliases ] && source .bash_aliases\n \
                    gringo {PROG_NAIVE} {params_str} -o smodels \
                    | tee 1> >(wc -l) \
                        >({{ TIMEFORMAT=%R; time breakID -asp; }} 2> \
                                >(grep -oP ':.*|^[\d.]+$' \
                                    | grep -oP '[0-9]+(\.[0-9]*)?' >&2) \
                                | tee 1> >(wc -l) \
                                    >(tail -n +2 \
                                        | cat - <(echo '0') \
                                        | clasp 0 --project=show --quiet=2,0,2 -t {num_threads}) \
                        ) 2>&1 \
                    | grep -oP ':.*|^\d+$' 2>&1 \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout)
        numbers1 = output.stdout.splitlines()
        numbers2 = output.stderr.splitlines()
        print_cmd(cmd, numbers1, numbers2)

        ground_prog_size = to_num(numbers1, 0, int)
        num_models = to_num(numbers1, 2, int)
        runtime = to_num(numbers1, 4, float)
        solving_time = to_num(numbers1, 5, float)
        first_model_time = to_num(numbers1, 7, float)
        unsat_time = to_num(numbers1, 8, float)
        total_cpu_time = to_num(numbers1, 9, float)
        grounding_time = total_cpu_time - solving_time - first_model_time - unsat_time

        breakid_runtime = to_num(numbers2, -1, float)
        breakid_ground_prog_size = to_num(numbers1, 1, int) - ground_prog_size
        breakid_num_symmetry_generators = to_num(numbers2, 0, int)
        breakid_auxiliary_variables = to_num(numbers2, -4, int)
        breakid_symmetry_breaking_clauses = to_num(numbers2, -3, int)

        return Result(ground_prog_size, num_models, runtime, solving_time, first_model_time, unsat_time, grounding_time, total_cpu_time,
            breakid_runtime, breakid_ground_prog_size, breakid_num_symmetry_generators, breakid_auxiliary_variables, breakid_symmetry_breaking_clauses)

    except TimeoutExpired:
        return None


def measure_molgen(carbons: int, hydrogens: int, oxygens: int, nitrogens: int) -> Union[Result, None]:
    params_str = f"C{carbons}H{hydrogens}{f'O{oxygens}' if oxygens != 0 else ''}{f'N{nitrogens}' if nitrogens != 0 else ''}"

    try:
        numbers = run((cmd := ["bash", "-O", "expand_aliases", "-c", f"[ -f .bash_aliases ] && source .bash_aliases\n \
                    {{ TIMEFORMAT=%R; time molgen {params_str} -v; }} 2>&1 \
                    | tail -2 \
                    | grep -oP '[0-9]+(\.[0-9]*)?'"]), capture_output=True, text=True, timeout=timeout).stdout.splitlines()
        print_cmd(cmd, numbers)

        num_models = to_num(numbers, 0, int)
        total_cpu_time = to_num(numbers, 1, float)
        return Result(None, num_models, None, None, None, None, total_cpu_time, total_cpu_time)

    except TimeoutExpired:
        return None


RESULTS = Dict[str, List[Result]]


def evaluate(func, cycles: int, oxygens: int, max_carbons: int, repetitions: int, results: RESULTS = dict()) -> RESULTS:
    for carbons in range(1, max_carbons+1):
        print(f"--- carbons={carbons} ---", file=stderr)
        hydrogens = 2*carbons+2-2*cycles
        key = str((carbons, hydrogens, oxygens))
        result_list = results[key] if key in results else []
        for k in range(repetitions - len(result_list)):
            result = func(carbons, hydrogens, oxygens, 0)
            if result is not None:
                print(f"{key};{result.to_json()}")
                result_list.append(result)
            else:
                (f"{key};None")
                break
        else:
            results[key] = result_list
            continue
        break
    return results


def evaluate_ground_prog_size(func, cycles: int, oxygens: int, max_carbons: int, repetitions: int, results: RESULTS) -> RESULTS:
    for carbons in range(1, max_carbons+1):
        print(f"--- ground_prog_size for carbons={carbons} ---", file=stderr)
        hydrogens = 2*carbons+2-2*cycles
        key = str((carbons, hydrogens, oxygens))
        ground_prog_size = results[key][0].ground_prog_size if key in results and len(results[key]) > 0 else None
        if ground_prog_size is None:
            ground_prog_size = func(carbons, hydrogens, oxygens, 0)
            if ground_prog_size is not None:
                print(f"{key};ground_prog_size;{ground_prog_size}")
            else:
                print(f"{key};ground_prog_size;None")
                break
        if key in results:
            for r in results[key]:
                r.ground_prog_size = ground_prog_size
        else:
            results[key] = []
            for i in range(repetitions):
                results[key].append(Result(ground_prog_size))
    return results


@dataclass_json
@dataclass
class ResultSet:
    data: Dict[str, Dict[str, RESULTS]]


def diagram(filename: str, title: str, label: str, data: Dict[str, List[Tuple[List[Union[int, float]], List[List[Union[int, float]]]]]], stack_labels: List[Tuple[str,str]]):
    plt.rcParams["figure.figsize"] = [8.00, 3.50]
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams['font.size'] = 17
    prop_cycle_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax = plt.subplots()

    # bar chart
    width = 0.9 / len(data.keys())
    stack_color_association = []
    for i, (attribute, stack) in enumerate(data.items()):
        offset = width * i
        last_measurement = None
        facecolor = None
        for j, (measurement, errorbars) in enumerate(stack):
            if j > 0:
                args = {
                    'bottom': last_measurement,
                    'fill': False,
                    'hatch': stack_labels[j][1],
                    'edgecolor': facecolor,
                }
            else:
                args = {
                    'label': attribute,
                    'color': prop_cycle_colors[i % len(prop_cycle_colors)],
                    'edgecolor': prop_cycle_colors[i % len(prop_cycle_colors)],
                }

            rects = ax.bar(
                np.arange(len(measurement)) + offset,
                measurement,
                width,
                **args,
                linewidth=1,
                yerr=np.transpose(errorbars) if not all(map(lambda x: x[0]==0 and x[1]==0, errorbars)) else None,
                capsize=2,
                ecolor='grey')
            if j == 0:
                last_measurement = measurement
                facecolor = rects[0].get_facecolor()
            else:
                last_measurement = [sum(m) for m in zip(last_measurement, measurement)]
            if len(stack_color_association) <= j:
                stack_color_association.append([])
            if not all(map(lambda x: x==0, measurement)):
                stack_color_association[j].append(facecolor)

    box = ax.get_position()
    ax.set_position([box.x0 - box.width * 0.05, box.y0 + box.height * 0.3,
                    box.width * 1.1, box.height * 0.75])

    # legend
    legend1 = ax.legend(loc='upper left', ncols=len(data.keys()), bbox_to_anchor=(-0.1, -0.1))

    num_stack_labels = 0
    for j, stack_label in enumerate(stack_labels):
        if stack_label[0] != "":
            num_stack_labels += 1
            label_color = stack_color_association[j][0] if len(stack_color_association[j]) == 1 else 'black'
            ax.bar(0, 0, 0,
                label=stack_label[0],
                fill=(j == 0),
                edgecolor=label_color,
                linewidth=1,
                color=label_color,
                hatch=stack_labels[j][1])

    if num_stack_labels > 0:
        handles, labels = plt.gca().get_legend_handles_labels()
        ax.legend(handles[-num_stack_labels:], labels[-num_stack_labels:],
            loc='upper left', ncols=num_stack_labels, bbox_to_anchor=(-0.1,-0.35))
    plt.gca().add_artist(legend1)

    # y axis logarithmic
    ax.set_yscale('log')
    ax.set_ylim(auto=True)
    ax.set_ylabel(label)
    #ax.yaxis.grid(True)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    # x axis ticks
    max_x = max((len(l[0][0]) for l in data.values()))
    ax.set_xticks(np.arange(max_x) + ((len(data.keys())-1)/2)*width, np.arange(1, max_x+1))

    # plot title
    plt.title(title)

    # save the figure in PDF format and close it
    if not exists(diagram_path):
        makedirs(diagram_path)

    plt.savefig(f"{diagram_path}/diagram_{'-'.join(f'{filename}_{label}'.split()).lower()}_log.pdf")
    plt.savefig(f"{diagram_path}/diagram_{'-'.join(f'{filename}_{label}'.split()).lower()}_log.svg")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate ASP-programs for chemical structure exploration.',
                                     epilog='Make sure to pipe stdout into a file, so that recovery is possible (e.g.: python eval.py -d -t data.csv 1>> data.csv).')
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
    eval_group.add_argument('--series', default=len(series), type=int,
                            help=f'Number of series to measure from {series} (default: all)')
    eval_group.add_argument('--carbons_limit', default=carbons_limit, type=int,
                            help=f'Maximum number of carbon atoms to consider (default: {carbons_limit})')
    eval_group.add_argument('--repetition_count', default=repetition_count, type=int,
                            help=f'Number of repetitions per experiment (default: {repetition_count})')
    args = parser.parse_args()
    if args.debug:
        _DEBUG = True
    num_threads = args.num_threads
    timeout = args.timeout
    del series[args.series:]
    carbons_limit = args.carbons_limit
    repetition_count = args.repetition_count

    if isfile(path):
        with open(path, "r") as fp:
            result_set = ResultSet.from_json(fp.read())
    else:
        result_set = ResultSet(dict())
        if args.try_recover is not None:
            with open(args.try_recover, "r") as fp:
                print("Recovery...", file=stderr)
                series_key = None
                name = None
                for line in fp.readlines():
                    parts = line.split(";")
                    if len(parts) == 2 and parts[0] == "===":
                        series_key = parts[1].strip()
                        if series_key not in result_set.data:
                            result_set.data[series_key] = dict()
                    elif series_key is not None and len(parts) == 2 and parts[0] == "---":
                        name = parts[1].strip()
                        if name not in result_set.data[series_key]:
                            result_set.data[series_key][name] = dict()
                    elif series_key is not None and name is not None and len(parts) > 1:
                        key = parts[0]
                        if key not in result_set.data[series_key][name]:
                            result_set.data[series_key][name][key] = []
                        if len(parts) == 3 and parts[1] == "ground_prog_size":
                            ground_prog_size = int(parts[2])
                            if len(result_set.data[series_key][name][key]) != 0:
                                for r in result_set.data[series_key][name][key]:
                                    r.ground_prog_size = ground_prog_size
                            else:
                                for i in range(repetitions):
                                    result_set.data[series_key][name][key].append(Result(ground_prog_size))
                        elif len(parts) == 2:
                            result_set.data[series_key][name][key].append(Result.from_json(parts[1]))
                print(pp.pformat(result_set), file=stderr)
        for cycles, oxygens in series:
            print(f"{GREEN}=== cycles={cycles}, oxygens={oxygens} ==={C_END}", file=stderr)
            series_key = str((cycles, oxygens))
            print(f"===;{series_key}")
            run_data = result_set.data[series_key] if series_key in result_set.data else dict()
            for name, (func, func_ground_prog_size) in zip(["Naive", "Canonical", "BreakID", "Smiles", "Molgen", "SBASS"], \
                                                           [(measure_naive, measure_ground_prog_size_naive), (measure_canonical, measure_ground_prog_size_canonical), (measure_breakid, None), (measure_smiles, measure_ground_prog_size_smiles), (measure_molgen, None), (measure_sbass, None)]):
                print(f"{GREEN}~~~ {name} ~~~{C_END}", file=stderr)
                print(f"---;{name}")
                run_data[name] = evaluate(func, cycles, oxygens, max_carbons=carbons_limit, repetitions=repetition_count, results=run_data[name] if name in run_data else dict())
                if func_ground_prog_size is not None:
                    evaluate_ground_prog_size(func_ground_prog_size, cycles, oxygens, max_carbons=carbons_limit, repetitions=repetition_count, results=run_data[name])
            result_set.data[series_key] = run_data
        with open(path, "w") as fp:
            fp.write(result_set.to_json())

    if args.print_results:
        #tmp = dict()
        #for n, d in result_set.data[str((0,0))].items():
        #    tmp[n] = [(v[0].grounding_time, v[0].solving_time) for k, v in sorted(d.items(), key=lambda x: eval(x[0]))]
        #print(pp.pformat(tmp), file=stderr)
        print(pp.pformat(result_set), file=stderr)

    if args.render_diagrams:
        for series, series_data in result_set.data.items():
            (cycles, oxygens) = eval(series)

            series_name = f"cycles={cycles} oxygens={oxygens}"

            formula_cycles = ' + 2' if cycles == 0 else ('' if cycles == 1 else f' - {cycles*2}')
            formula_oxygens = '' if oxygens == 0 else ('O' if oxygens == 1 else f'O_{oxygens}')
            title = f"$C_xH_{{2 \cdot x{formula_cycles}}}{formula_oxygens}$"

            for dependent_name, dependents, stack_labels, in zip(["Number of models", "Total runtime", "Ground program size"], \
                                                                #[["num_models"], ["grounding_time", "first_model_time", "solving_time", "unsat_time", "breakid_runtime"], ["ground_prog_size", "breakid_ground_prog_size"]], \
                                                                #[[], ["Grounding", "First model", "Solving", "Unsat", "BreakID run"], ["", "BreakID SBCs"]]):
                                                                [["num_models"], ["grounding_time", "solving_time", "breakid_runtime"], ["ground_prog_size", "breakid_ground_prog_size"]], \
                                                                [[], [("Solving", ''), ("Grounding", 'XXXXX'), ("BreakID run", '.....')], [("", ''), ("SBASS / BreakID SBCs", '.....')]]):
                def process(lst):
                    median_list = [median(l) for l in lst if None not in l]
                    if len(median_list) == 0:
                        return None
                    error_list = [[med - min(l), max(l) - med] for med, l in zip(median_list, lst[:len(median_list)])]
                    return (median_list, error_list)

                key_map = lambda x: x if x != "Smiles" else "Our encoding"
                values = { key_map(k): p for k, v in series_data.items() \
                    if (p := [l for d in dependents \
                        if (l := process([[getattr(vvv, d) for vvv in vv] for _, vv in sorted(v.items(), key=lambda x: eval(x[0]))])) is not None]) != [] }

                if "Our encoding" in values.keys():
                    max_len = len(values["Our encoding"][0][0])
                    for k in values.keys():
                        for l in values[k]:
                            if len(l[0]) > max_len:
                                del l[0][max_len:]
                            if len(l[1]) > max_len:
                                del l[1][max_len:]

                if args.print_results:
                    print(dependent_name, file=stderr)
                    print(pp.pformat(values), file=stderr)

                diagram(series_name, title, dependent_name, values, stack_labels)

if False: #__name__ == "__main__":
    for cycles in [0,1]:
        oxygens = 1
        print(f"===;{str((cycles, oxygens))}")
        f_num_models = open(f"results/performance_cycles={cycles}_oxygens={oxygens}_num_models.csv", "r")
        f_runtime = open(f"results/performance_cycles={cycles}_oxygens={oxygens}_runtime.csv", "r")

        lines_num_models = [l[:-1].split(',') for l in f_num_models.readlines()]
        lines_runtime = [l[:-1].split(',') for l in f_runtime.readlines()]
        max_carbons = int(lines_num_models[-1][0])
        repetitions = repetition_count

        lines_num_models[0] += [f'sbass_{i}' for i in range(1,repetitions+1)]
        lines_runtime[0] += [f'sbass_{i}' for i in range(1,repetitions+1)]

        _DEBUG = True

        done = False
        for carbons in range(1, max_carbons+1):
            print(f"--- carbons={carbons} ---", file=stderr)
            hydrogens = 2*carbons+2-2*cycles
            key = str((carbons, hydrogens, oxygens))
            for k in range(repetitions):
                result = None if done else measure_sbass(carbons, hydrogens, oxygens, 0)
                if result is not None:
                    print(f"{key};{result.to_json()}")
                    lines_num_models[carbons].append(result.num_models)
                    lines_runtime[carbons].append(result.runtime)
                else:
                    if not done:
                        print(f"{key};None")
                        done = True
                    lines_num_models[carbons].append(None)
                    lines_runtime[carbons].append(None)

        open(f"results/new_perperformance_cycles={cycles}_oxygens={oxygens}_num_models.csv", "w").write('\n'.join([','.join([str(v) for v in l]) for l in lines_num_models]))
        open(f"results/new_perperformance_cycles={cycles}_oxygens={oxygens}_runtime.csv", "w").write('\n'.join([','.join([str(v) for v in l]) for l in lines_runtime]))
