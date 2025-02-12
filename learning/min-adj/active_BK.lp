lt(I,J) :- node(I), node(J), I < J.
eq(I,I) :- node(I).
lte(M,N) :- edge(_,_,M), edge(_,_,N), M <= N.
