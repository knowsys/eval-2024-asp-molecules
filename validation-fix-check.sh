#!/usr/bin/env bash

# Show the solver invocations.
set -x

# Allow killing script with Ctrl-C, instead of individual command.
trap "exit" INT

PROG_SMILES="smiles.lp"

[ -f "$PROG_SMILES" ] || ./prepare-asp-programs.sh

# Iterate all unsatisfiable factbases.
for f in validation/unsat_*.lp; do
    # Exit on first UNSATISFIABLE factbase.
    # Project, so that different choices of atom_map(_, _, _) do not register, and the correct count of isomorphic structures is reported.
    clingo 0 $PROG_SMILES smiles-to-edge_min.lp smiles-check_min.lp <(cat $f) --project=show \
        | python smiles-vis.py -c
    [[ "${PIPESTATUS[0]}" == "20" ]] && exit 1
done

echo "ALL OK !!!"
