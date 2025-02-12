%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra Predicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The meaning of `same(J,I,J,K)` is that
% -- J and K are of the same color and
% -- all L<I have an edge to K if and only if they have an edge (of equal multiplicity) to J.
tmp_same(0,I,J,J,K) :- node(I), type(J,E), type(K,E), I < J, J < K.
same(I,J,J,K) :- tmp_same(I-1,I,J,J,K).

% The meaning of `same(J,I,I,K)` is that
% -- I and K are of the same color and
% -- all L<J,L!=I,L!=K have an edge to K if and only if they have an edge (of equal multiplicity) to I.
tmp_same(0,I,J,I,K) :- type(I,E), node(J), type(K,E), I < J, I < K, J != K.
tmp_same(I,I,J,I,K) :- tmp_same(I-1,I,J,I,K).
tmp_same(K,I,J,I,K) :- tmp_same(K-1,I,J,I,K).
same(I,J,I,K) :- tmp_same(J-1,I,J,I,K).

% Since the graph is undirected, we only look at 'forward' edges.
tmp_same(L,I,J,T,K) :- tmp_same(L-1,I,J,T,K), edge(L,T, M), edge(L,K, M).
tmp_same(L,I,J,T,K) :- tmp_same(L-1,I,J,T,K), not edge(L,T), not edge(L,K), node(L).

%%#show same/4.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SBCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We would expect a hypothesis like:

%%:- not edge(I,J), edge(I,K), same(I,J,J,K).
% Intuitively, this means that J and K can be swapped to produce a lecicographically smaller isomorphic graph.

%%:- not edge(I,J), edge(J,K), same(I,J,I,K).
% Intuitively, this means that I and K can be swapped to produce a lecicographically smaller isomorphic graph.

%%#show prune/3.
