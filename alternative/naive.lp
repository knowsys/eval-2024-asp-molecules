% INPUT: C_{c}H_{h},O_{o}N_{n}

#const c=0.
#const h=0.
#const n=0.
#const o=0.

molecular_formula("C", c) :- c != 0.
molecular_formula("H", h) :- h != 0.
molecular_formula("N", n) :- n != 0.
molecular_formula("O", o) :- o != 0.

% The element predicate contains the atomic number and the valence for each element
element("C", 6, 4) :- molecular_formula("C", C), C > 0.
element("H", 1, 1) :- molecular_formula("H", H), H > 0.
element("N", 7, 3) :- molecular_formula("N", N), N > 0.
element("O", 8, 2) :- molecular_formula("O", O), O > 0.

% ---------------------------------------------

% Construct an undirected graph,
% which represents the molecule

node(1..V) :- V = #sum { C : molecular_formula(E, C), E != "H" }.

{ edge(X, Y) : node(X), node(Y), X < Y }.
edge(Y, X) :- edge(X, Y).

% Assign an element to each node in the graph
%
% We order elements by their atomic number and ignore hydrogen,
% meaning that the first c atoms will be carbon, the next n nitrogen,
% and the final o atoms will be assigned to oxygen.
type(START+(1..C), E) :- molecular_formula(E, C), element(E, ATOMIC_NR, _), E != "H",
                         START = #sum { C', E' : molecular_formula(E', C'),
                                                 element(E', A, _), A < ATOMIC_NR, E' != "H" }.

% Make sure the constructed graph is connected by checking if every node is reachable from the fist node

reachable(1).
reachable(Y) :- reachable(X), edge(X, Y).

:- not reachable(X), node(X).

% For each edge, guess a (unique) bond-multiplicity

1{ edge(X, Y, 1..3) }1 :- edge(X, Y), X < Y.
edge(Y, X, M) :- edge(X, Y, M).

% Compute the degree of each node
degree(N, D) :- node(N), D = #sum { C, X : edge(N, X, C) }.

% No atom may have a larger degree than the valance of its assigned element
:- node(N), type(N, E), degree(N, D), element(E, _, V), D > V.

% Ensure that the required amount of bonds has been placed
:- EDGE_COUNT = #sum { M, X, Y : edge(X, Y, M), X < Y },
   VALENCE_SUM = #sum { V*C, E : molecular_formula(E, C), element(E, _, V), E != "H" },
   molecular_formula("H", H_COUNT),
   EDGE_COUNT != (VALENCE_SUM - H_COUNT)/2.

%#show type/2.
#show edge/3.
