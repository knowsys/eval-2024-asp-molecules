% A version of this file without comments was produced with the following command:
% $> sed '/^[ \t]*%/d; s/[ \t]*%.*$//' smiles-check-standard-valence-test.lp | cat -s > smiles-check-standard-valence-test_min.lp

% STANDARD VALENCE TEST
% Tests, whether the structure can be attained with standard valence values.
% (Standalone; don't combine with smiles.lp)

% [see also smiles.lp]
% Compute, how many atoms have which valnce value.
valence(1;3;4).
valence_count(VALENCE, X) :- X = #sum{ COUNT, ELEMENT : molecular_formula(ELEMENT, COUNT), element(ELEMENT, _, VALENCE) }, valence(VALENCE), not element_valence(_, VALENCE).
element_valence(ELEMENT, VALENCE) :- element(ELEMENT, _, VALENCE), element_valence(ELEMENT, _).
1{ element_valence_count(ELEMENT, VALENCE, 0..COUNT) }1 :- element_valence(ELEMENT, VALENCE), molecular_formula(ELEMENT, COUNT).
:- molecular_formula(ELEMENT, COUNT), COUNT != #sum{ C, VALENCE : element_valence_count(ELEMENT, VALENCE, C) }, element_valence(ELEMENT, _).
valence_count(VALENCE, X) :- X = #sum{ COUNT, ELEMENT : element_valence_count(ELEMENT, VALENCE, COUNT) }, element_valence(_, VALENCE).
% Compute, how many cycles / multi-bonds there should be, based on the hydrogen count in the sum formula.
num_cycles_or_multi_bonds((2*C+2 + N - H) / 2) :- valence_count(4, C), valence_count(3, N), valence_count(1, H).
% SENIOR rules (adaptation)
:- valence_count(4, C), valence_count(3, N), valence_count(1, H), 2*C+2 + N < H. % --> cannot have too many atoms with valence 1
:- valence_count(4, C), valence_count(3, N), valence_count(1, H), (2*C+2 + N - H) \ 2 != 0. % --> cannot have "half a cycle"

% Sum up the multi-bonds.
num_multi_bonds(X) :- X = #sum{ MULTIPLICITY-1, IDX_1, ELEMENT_1, IDX_2, ELEMENT_2 : required_bond(IDX_1, ELEMENT_1, IDX_2, ELEMENT_2, MULTIPLICITY) }.

% Sort the required bonds lexicographically.
required_bond_count(X) :- X = #count{ (I1, E1, I2, E2) : required_bond(I1, E1, I2, E2, M) }.
ord(1, T) :- T = #min{ edge((I1, E1), (I2, E2)) : required_bond(I1, E1, I2, E2, M), (I1, E1) < (I2, E2);
                       edge((I1, E1), (I2, E2)) : required_bond(I2, E2, I1, E1, M), (I1, E1) < (I2, E2) }.
ord(O+1, T) :- ord(O, T'),
               T = #min{ edge((I1, E1), (I2, E2)) : required_bond(I1, E1, I2, E2, M), T' < edge((I1, E1), (I2, E2)), (I1, E1) < (I2, E2);
                         edge((I1, E1), (I2, E2)) : required_bond(I2, E2, I1, E1, M), T' < edge((I1, E1), (I2, E2)), (I1, E1) < (I2, E2) },
               O < O_MAX, required_bond_count(O_MAX).

% Put an (arbitrary) tree over the required bond relation.
tree(1, A, B) :- ord(1, edge(A, B)).

tree(O, A, B) :- tree(O', A, C), O' < O, A < B,
                 ord(O, edge(C, B)),
                 { tree(O'', A, B) : O'' < O }0.
tree(O, A, B) :- tree(O', A, C), O' < O, A < B,
                 ord(O, edge(B, C)),
                 { tree(O'', A, B) : O'' < O }0.

tree(O, A, B) :- ord(O, edge(A, B)), A < B,
                 { tree(O', _, A) : O' < O }0,
                 { tree(O'', _, B) : O'' < O }0.

tree(O, A, B) :- tree(O', A, C), tree(O'', C, B), A < B,
                 O = #max{ O1 : O1=O'; O2 : O2=O'' },
                 { tree(O''', A, B) : O''' < O }0, { tree(O, A', B) : A' < A }0.
tree(O, A, B) :- tree(O', A, C), tree(O'', B, C), A < B,
                 O = #max{ O1 : O1=O'; O2 : O2=O'' },
                 { tree(O''', A, B) : O''' < O }0, { tree(O, A', B) : A' < A }0.

cycle(C, B) :- tree(O', A, C), O' < O, A < B,
               ord(O, edge(C, B)),
               tree(O'', A, B), O'' < O.
cycle(B, C) :- tree(O', A, C), O' < O, A < B,
               ord(O, edge(B, C)),
               tree(O'', A, B), O'' < O.

num_cycles(X) :- X = #count{ A, B : cycle(A, B) }.

:- num_multi_bonds(X1), num_cycles(X2), not num_cycles_or_multi_bonds(X1 + X2).

% Additionally, check that each atom's degree is at most its valence.
:- required_bond(I, E, _, _, _), element(E, _, VALENCE), VALENCE < IN_DEG + OUT_DEG,
   IN_DEG = #sum{ M, I', E' : required_bond(I', E', I, E, M) },
   OUT_DEG = #sum{ M, I', E' : required_bond(I, E, I', E', M) }.
:- required_bond(_, _, I, E, _), element(E, _, VALENCE), VALENCE < IN_DEG + OUT_DEG,
   IN_DEG = #sum{ M, I', E' : required_bond(I', E', I, E, M) },
   OUT_DEG = #sum{ M, I', E' : required_bond(I, E, I', E', M) }.

#show num_cycles/1.
#show num_multi_bonds/1.
