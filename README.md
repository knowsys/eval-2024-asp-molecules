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

## Correctness validation

```bash
$GENMOL selfcheck -t 20 -o validation.csv -- -t 2
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

## Performance evaluation

```bash
python eval.py -d -t data.csv 1>> data.csv
```

See ``python eval.py --help` for further options.

## Symmetry-breaking evaluation

```bash
python eval-sb.py -d -t sb-data.csv 1>> sb-data.csv
```

See ``python eval-sb.py --help` for further options.
