# eval-2024-asp-molecules

Evaluation files for an ASP-based tool that enumerates molecule shapes for a given sum fomula

## Usage

To enumerate all molecule structures matching a sum formula of the form $C_qH_rN_sO_t$, use:

```bash
clingo 0 smiles.lp --const c=q --const h=r --const n=s -const o=t
```

To convert the models in a more readable SMILES representation, pipe the output into the enclosed visualization script:

```bash
clingo 0 smiles.lp --const c=q --const h=r --const n=s -const o=t | python smiles-vis.py
```

See `python smiles-vis.py --help` for further options.

### Examples

```bash
# C7H16, acyclic
clingo 0 smiles.lp --const c=7 --const h=16 | python smiles-vis.py -c
# C7H14, 1 cycle
clingo 0 smiles.lp --const c=7 --const h=14 | python smiles-vis.py -c
```

To check which of the emmitted structures are isomorphic, use e.g.:

```bash
clingo 0 smiles.lp smiles-to-edge.lp --const c=7 --const h=14 | python smiles-vis.py -c
```

(Note that this is only feasible for relatively small model counts.)

## Preparation

For the correctness validation and the symmetry-breaking evaluation,
the Rust-based interface program is needed.
You can either build it locally or use the Docker image.

To build locally - ommitting the Yew frontend as it is not needed for the experiments, run the following commands:

```bash
git submodule update --init
carbo build -r
GENMOL=./tool//target/x86_64-unknown-linux-gnu/release/genmol
```

(Refer to the tool's [README](todo) for further details.)

To use the pre-built Docker image:

```bash
GENMOL_IMAGE="todo"
docker pull $GENMOL_IMAGE
GENMOL="docker run $GENMOL_IMAGE"
```

You can also use the tool to translate a given sumformula into a factbase.
This enables enumeration of molecules containing any main-group element, e.g.:

```bash
clingo 0 smiles.lp <($GENMOL to-factbase -f Si2C5H14) | python smiles-vis.py -c
```

See `$GENMOL --help` for a full list of supported options.

## Correctness validation

To assess the correctness of our encoding, the validation performs esperiments
to see, whether relevant molecules collected from Wikidata can be found
with `smiles.lp` given their sum formula.
It also reports the number of isomorphic models per compound.

This will write results to a `validation.csv` data file.
Each solver invocation is performed with a timeout of 40 seconds.
Note that you can terminate the program at any time with Ctrl+C,
which will cause it to report the number of satisfiable and unsatsifiable
cases encountered, as well as the average number isomorphic models.
Additionally, it reports the number of compounds with subgroup elements
and non-standard valences, which are skipped.
You can continue the validation from where you terminated it, by simply
re-running the same command. In case you increase the timeout value,
compounds which exhausted the previous timeout will be re-checked automatically.

```bash
$GENMOL selfcheck -t 40 -o validation.csv -- -t 2
```

See `$GENMOL selfcheck --help` for further options.

To regenerate the dataset of relevant chemical compounds, you can run:

```bash
curl https://query.wikidata.org/sparql \\
    --header 'Accept: text/csv' \\
    --data-urlencode query='SELECT ?sumformula ?smiles ?inchi ?name ?qid
WHERE
{{
  ?chemical wdt:P233 ?smiles;
            wdt:P234 ?inchi;
            wdt:P274 ?sumformula.
  ?article schema:about ?chemical;
             schema:isPartOf <https://en.wikipedia.org/>.
  BIND (REPLACE(STR(?chemical), \"^.*/([^/]*)$\", \"$1\") as ?qid).
  BIND (REPLACE(STR(?article), \"^.*/([^/]*)$\", \"$1\") as ?name).
}}' | python chemdata-sort.py > chemdata-sort.csv
```

This dataset consists of all chemical compounds on Wikidata, which are associated with
a SMILES and InChi, as well as a sum formula. To restrict to - in some sense - relevant
compounds, we only consider compounds having a matching article in english Wikipedia.
The `chemdata-sort.py` script sorts the dataset by atom count.

## Performance evaluation

We evaluate the performance of our encoding against a naive graph-based encoding,
symmetry-breaking via caninical graphs, automated symmetry-breaking using __BreakID__,
and the pre-existing commercial tool __Molgen__,
w.r.t. ground program size, total runtime, and number of models.

To this purpose, sum formulae of the form $C_xH_{2 \cdot x + 2 - 2 \cdot c}O_o$ are considered
for increasing carbon atom counts.
Here, $c$ represents the number of cycles or multi-bonds and $o$ the oxygen count in molecules matching the formula.
Series are collected for $\langle c,o \rangle \in \{0,1,2\}^2$.

This will produce a `results.json` data file as well as several diagrams for each series:

* `diagrams/diagram_cycles=<c>-oxygens=<o>_ground-program-size_log`,
* `diagrams/diagram_cycles=<c>-oxygens=<o>_number-of-models_log.pdf`, and
* `diagrams/diagram_cycles=<c>-oxygens=<o>_total-runtime_log.pdf`

You can terminate it with Ctrl+C at any time and continue by re-running the command.
This picks up your intermediate results from the recovery file `sb-data.csv`.

```bash
python eval.py -d -t data.csv -r 1>> data.csv
```

See ``python eval.py --help` for further options.

## Symmetry-breaking evaluation

We compare the number of models our encoding finds agaisnst the number of structures reported by __Molgen__,
for real-world sum formulae collected from Wikidata (see `chemdata-sort.csv`).

This will produce a `sb-results.json` data file as well as `diagrams/diagram_number_of_models-comparison.pdf`.
You can terminate it with Ctrl+C at any time and continue by re-running the command.
This picks up your intermediate results from the recovery file `sb-data.csv`.

```bash
python eval-sb.py -d -t sb-data.csv -r 1>> sb-data.csv
```

See ``python eval-sb.py --help` for further options.
