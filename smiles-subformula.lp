% A version of this file without comments was produced with the following command:
% $> sed '/^[ \t]*%/d; s/[ \t]*%.*$//' smiles-subformula.lp | cat -s > smiles-subformula_min.lp

% SUBSTFORMULA

% Nondeterministically select atoms belonging to the fragment.
COUNT_FRAG{ fragment_atom(FRAGMENT_IDX, I) : symbol(I, ELEMENT) }COUNT_FRAG :- fragment_molecular_formula(FRAGMENT_IDX, ELEMENT, COUNT_FRAG), ELEMENT != "H".

% Connectivity of fragment atom via edge relation.
fragment_connected(FRAGMENT_IDX, I1, I2) :- fragment_atom(FRAGMENT_IDX, I1),
                                            fragment_atom(FRAGMENT_IDX, I2),
                                            edge(I1, I2, _).
fragment_connected(FRAGMENT_IDX, I1, I2) :- fragment_connected(FRAGMENT_IDX, I1, I),
                                            fragment_atom(FRAGMENT_IDX, I2),
                                            edge(I, I2, _).
fragment_connected(FRAGMENT_IDX, I1, I2) :- fragment_connected(FRAGMENT_IDX, I1, I),
                                            fragment_atom(FRAGMENT_IDX, I2),
                                            edge(I2, I, _).

% Ensure that all fragment atoms are connected.
:- fragment_atom(FRAGMENT_IDX, I1), fragment_atom(FRAGMENT_IDX, I2), I1 < I2,
   not fragment_connected(FRAGMENT_IDX, I1, I2).

% Ensure that the fragment has the required numner of (implicit) hydrogen atoms.
:- fragment_molecular_formula(FRAGMENT_IDX, "H", X),
   X != #sum{ FREE_BONDS, I : fragment_atom(FRAGMENT_IDX, I), free_bonds(I, FREE_BONDS) }.
