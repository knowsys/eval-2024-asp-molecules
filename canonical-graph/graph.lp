% INPUT: C_{c}H_{h},O_{o}

#const c=0.
#const h=0.
#const n=0.
#const o=0.

molecular_formula("C", c) :- c != 0.
molecular_formula("H", h) :- h != 0.
molecular_formula("N", n) :- n != 0.
molecular_formula("O", o) :- o != 0.

element("C", 6, 4) :- molecular_formula("C", C), C > 0.
element("H", 1, 1) :- molecular_formula("H", H), H > 0.
element("N", 7, 3) :- molecular_formula("N", N), N > 0.
element("O", 8, 2) :- molecular_formula("O", O), O > 0.

% ---------------------------------------------

% Multigraph
node(1..V) :- V = #sum { C : molecular_formula(E, C), E != "H" }.
type(START+(1..C), E) :- molecular_formula(E, C), element(E, ATOMIC_NR, _), E != "H",
                         START = #sum { C', E' : molecular_formula(E', C'),
                                                 element(E', A, _), A < ATOMIC_NR, E' != "H" }.

{ edge(X, Y) : node(X), node(Y), X < Y }.
edge(Y, X) :- edge(X, Y).

% Connected
reachable(1).
reachable(Y) :- reachable(X), edge(X, Y).

:- not reachable(X), node(X).

% Edge labels
1{ edge(X, Y, 1..3) }1 :- edge(X, Y), X < Y.
edge(Y, X, M) :- edge(X, Y, M).

% Degrees
degree(N, D) :- node(N), D = #sum { C, X : edge(N, X, C) }.
:- node(N), type(N, E), degree(N, D), element(E, _, V), D > V.
:- EDGE_COUNT = #sum { M, X, Y : edge(X, Y, M), X < Y },
   VALENCE_SUM = #sum { V*C, E : molecular_formula(E, C), element(E, _, V), E != "H" },
   molecular_formula("H", H_COUNT),
   EDGE_COUNT != (VALENCE_SUM - H_COUNT)/2.

%#show type/2.
#show edge/3.
