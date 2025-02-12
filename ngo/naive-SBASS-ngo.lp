#program base.
% INPUT: C_{c}H_{h},O_{o}N_{n}
#const c = 0.
#const h = 0.
#const n = 0.
#const o = 0.
% (do not check for zero, since we rely on all four predicates being present to substitute the #sum aggregates)
molecular_formula("C",c).
molecular_formula("H",h).
molecular_formula("N",n).
molecular_formula("O",o).
% The element predicate contains the atomic number and the valence for each element
element("C",6,4).
element("N",7,3).
element("O",8,2).
% ---------------------------------------------
%%%node(1..V) :- V = #sum { C, E : molecular_formula(E, C), E != "H" }.
% (rewritten without #sum aggregate; only for CHNO formulas)
size(((C+N)+O)) :- molecular_formula("C",C); molecular_formula("N",N); molecular_formula("O",O).
node(S) :- size(S).
node((C-1)) :- node(C); C > 1.
{ edge(X,Y): node(X), node(Y), X < Y }.
edge(Y,X) :- edge(X,Y).
% Assign an element to each node in the graph
%
% We order elements by their atomic number and ignore hydrogen,
% meaning that the first c atoms will be carbon, the next n nitrogen,
% and the final o atoms will be assigned to oxygen.
%%%type(START+(1..C), E) :- molecular_formula(E, C), element(E, ATOMIC_NR, _), E != "H",
%%%                         START = #sum { C', E' : molecular_formula(E', C'),
%%%                                                 element(E', A, _), A < ATOMIC_NR, E' != "H" }.
% (rewritten without #sum aggregate; only for CHNO formulas)
type(I,"C") :- node(I); molecular_formula("C",C); I <= C.
type(I,"N") :- node(I); molecular_formula("C",C); molecular_formula("N",N); C < I; I <= (C+N).
type(I,"O") :- node(I); molecular_formula("C",C); molecular_formula("N",N); (C+N) < I.
% Make sure the constructed graph is connected by checking if every node is reachable from the fist node
reachable(1).
reachable(Y) :- reachable(X); edge(X,Y).
#false :- not reachable(X); node(X).
% For each edge, guess a (unique) bond-multiplicity
1 <= { edge(X,Y,(1..3)) } <= 1 :- edge(X,Y); X < Y.
edge(Y,X,M) :- edge(X,Y,M).
% Compute the degree of each node
%%%degree(N, D) :- node(N), D = #sum { C, X : edge(N, X, C) }.
% (rewritten without #sum aggregate; only for CHNO formulas)
tmp_degree(0,N,0) :- node(N).
tmp_degree(I,N,(PRE+C)) :- tmp_degree((I-1),N,PRE); edge(N,I,C).
tmp_degree(I,N,(PRE+0)) :- tmp_degree((I-1),N,PRE); not edge(N,I); size(S); I <= S.
degree(N,D) :- tmp_degree(S,N,D); size(S).
% No atom may have a larger degree than the valance of its assigned element
#false :- type(N,E); degree(N,D); element(E,_,V); D > V.
% Ensure that the required amount of bonds has been placed
%%%:- EDGE_COUNT = #sum { M, X, Y : edge(X, Y, M), X < Y },
%%%   VALENCE_SUM = #sum { V*C, E : molecular_formula(E, C), element(E, _, V), E != "H" },
%%%   molecular_formula("H", H_COUNT),
%%%   EDGE_COUNT != (VALENCE_SUM - H_COUNT)/2.
% (rewritten without #sum aggregate; only for CHNO formulas)
tmp_edge_count(0,0,0).
__aux_1(C,I1,I2) :- edge(I1,I2,C); I1 < I2; size(S); I2 <= S.
tmp_edge_count(I1,I2,(PRE+C)) :- tmp_edge_count(I1,(I2-1),PRE); __aux_1(C,I1,I2).
tmp_edge_count(I1,I2,(PRE+0)) :- tmp_edge_count(I1,(I2-1),PRE); not edge(I1,I2,_); I1 < I2; size(S); I2 <= S.
tmp_edge_count(I1,I1,PRE) :- tmp_edge_count((I1-1),S,PRE); size(S); I1 < S.
edge_count(EDGE_COUNT) :- tmp_edge_count((S-1),S,EDGE_COUNT); size(S).
#false :- edge_count(EDGE_COUNT); molecular_formula("C",C); element("C",_,VC); molecular_formula("N",N); element("N",_,VN); molecular_formula("O",O); element("O",_,VO); molecular_formula("H",H_COUNT); EDGE_COUNT != (((((C*VC)+(N*VN))+(O*VO))-H_COUNT)/2).
