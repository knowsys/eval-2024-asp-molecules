% Symmetry breaking based on 
%
% Codish, M., Miller, A., Prosser, P., & Stuckey, P. J. 
% Constraints for symmetry breaking in graph representation (2014). 

sat(I, K, J) :- 
    type(I, T), type(J, T), type(K, T), type(L, T), J > I, J - I != 2,
    edge(I, K), edge(J, L), L < K, L != I.
sat(I, K, J) :- 
    type(I, T), type(J, T), type(K, T), J > I, J - I != 2,
    edge(I, K, N), edge(J, K, M), N <= M.

:- type(I, T), type(J, T), type(K, T), 
   edge(I, K), J > I, not sat(I, K, J), J - I != 2, K != J.
