% A version of this file without comments was produced with the following command:
% $> sed '/^[ \t]*%/d; s/[ \t]*%.*$//' smiles-check.lp | cat -s > smiles-check_min.lp

% CHECK

% Nondeterministically decide for an atom assignment.
1{ atom_map(I, 1..COUNT, ELEMENT) : molecular_formula(ELEMENT, COUNT), symbol(I, ELEMENT) }1 :- atom(I).
:- atom_map(I1, IDX_MAP, ELEMENT), atom_map(I2, IDX_MAP, ELEMENT), I1 < I2.

% Ensure that the generated bonds match the required bonds under the atom assignment.
:- not edge(I1, I2, MULTIPLICITY),
   not edge(I2, I1, MULTIPLICITY),
   required_bond(IDX_1, ELEMENT_1, IDX_2, ELEMENT_2, MULTIPLICITY),
   atom_map(I1, IDX_1, ELEMENT_1), atom_map(I2, IDX_2, ELEMENT_2).
:- edge(I1, I2, MULTIPLICITY),
   not required_bond(IDX_1, ELEMENT_1, IDX_2, ELEMENT_2, MULTIPLICITY),
   not required_bond(IDX_2, ELEMENT_2, IDX_1, ELEMENT_1, MULTIPLICITY),
   atom_map(I1, IDX_1, ELEMENT_1), atom_map(I2, IDX_2, ELEMENT_2).
