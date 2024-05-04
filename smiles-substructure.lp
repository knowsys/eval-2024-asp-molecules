% A version of this file without comments was produced with the following command:
% $> sed '/^[ \t]*%/d; s/[ \t]*%.*$//' smiles-substructure.lp | cat -s > smiles-substructure_min.lp

% SUBSTRUCTURE

% Project the fragment atoms.
2{ fragment_atom(FRAGMENT_IDX, IDX_1, ELEMENT_1); fragment_atom(FRAGMENT_IDX, IDX_2, ELEMENT_2) }2 :- fragment_bond(FRAGMENT_IDX, IDX_1, ELEMENT_1, IDX_2, ELEMENT_2, MULTIPLICITY).

% Nondeterministically decide for an atom assignment per fragment.
1{ fragment_atom_map(FRAGMENT_IDX, I, IDX, ELEMENT) : symbol(I, ELEMENT) }1 :- fragment_atom(FRAGMENT_IDX, IDX, ELEMENT).
:- fragment_atom_map(FRAGMENT_IDX, I_MAP, IDX_1, ELEMENT), fragment_atom_map(FRAGMENT_IDX, I_MAP, IDX_2, ELEMENT), IDX_1 < IDX_2.

% Ensure that the generated bonds match the fragment bonds under the atom assignment.
:- not edge(I1, I2, MULTIPLICITY),
   not edge(I2, I1, MULTIPLICITY),
   fragment_bond(FRAGMENT_IDX, IDX_1, ELEMENT_1, IDX_2, ELEMENT_2, MULTIPLICITY),
   fragment_atom_map(FRAGMENT_IDX, I1, IDX_1, ELEMENT_1), fragment_atom_map(FRAGMENT_IDX, I2, IDX_2, ELEMENT_2).
:- not weak_subgraph(FRAGMENT_IDX),
   edge(I1, I2, MULTIPLICITY),
   not fragment_bond(FRAGMENT_IDX, IDX_1, ELEMENT_1, IDX_2, ELEMENT_2, MULTIPLICITY),
   not fragment_bond(FRAGMENT_IDX, IDX_2, ELEMENT_2, IDX_1, ELEMENT_1, MULTIPLICITY),
   fragment_atom_map(FRAGMENT_IDX, I1, IDX_1, ELEMENT_1), fragment_atom_map(FRAGMENT_IDX, I2, IDX_2, ELEMENT_2).
