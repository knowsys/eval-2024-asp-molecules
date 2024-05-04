#!/usr/bin/env bash

# smiles.lp
sed '1,18d; 60,72d; s/[^%].*@min_main_chain_len.*/1{ main_chain_len(MIN_LEN..ATOM_COUNT) }1 :- non_hydrogen_atom_count(ATOM_COUNT), min_main_chain_len(MIN_LEN)./; /^[ \t]*%/d; s/[ \t]*%.*$//; /^[ \t]*# /d' smiles.lp | cat -s > smiles_min.lp
sed '60,72d; s/[^%].*@min_main_chain_len.*/1{ main_chain_len(MIN_LEN..ATOM_COUNT) }1 :- non_hydrogen_atom_count(ATOM_COUNT), min_main_chain_len(MIN_LEN)./; /^[ \t]*%/d; s/[ \t]*%.*$//; /^[ \t]*# /d; 7i #const m=0. min_main_chain_len(m) :- m != 0.' smiles.lp | cat -s > smiles_eval.lp
# smiles-to-edge.lp
sed '/#show/d; /^[ \t]*%/d; s/[ \t]*%.*$//' smiles-to-edge.lp | cat -s > smiles-to-edge_min.lp
# smiles-check.lp
sed '/^[ \t]*%/d; s/[ \t]*%.*$//' smiles-check.lp | cat -s > smiles-check_min.lp
