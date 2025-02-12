% Additional active background knowledge:

sat(I, K, J) :- 
    type(I, T), type(J, T), type(K, T), type(L, T), J > I, J - I != 2,
    edge(I, K), edge(J, L), L < K, L != I.
sat(I, K, J) :- 
    type(I, T), type(J, T), type(K, T), J > I, J - I != 2,
    edge(I, K, N), edge(J, K, M), N <= M.

% In theory, this would permit the constraint of the graph symmetry breking of Codish et al. (see `lex.lp`) as a hypothesis.

% Run with the following ILASP options:
% --max-rule-length=7 -ml=7 -na

#modeb(1, sat(var(node),var(node),var(node))).
#modeb(3, type(var(node),var(type)), (positive)).
#modeb(1, edge(var(node),var(node)), (anti_reflexive,symmetric)).
#modeb(1, lt(var(node),var(node)), (positive)).
#modeb(1, eq(var(node),var(node)), (symmetric,negative)).

#maxv(4).

#bias("
:- body(lt(A,B)), body(naf(eq(A,B))).
:- body(lt(A,B)), body(naf(eq(B,A))).

:- #count{ B : body(B) } < 7.

% To see whether the Codish conditions are contained, exclude all other rules.
%     ilasp --max-rule-length=7 -ml=7 -na min-adj/ILASP_iterative_step_1_SMALL.lp -st | clingo 0
%%rule3 :- #count{ B : body(B) } == 7, body(type(I,T)), body(type(J,T)), body(type(K,T)), body(edge(I,K)), body(lt(I,J)), body(naf(sat(I,K,J))), body(naf(eq(K,J))).
%%:- not rule3.
").

% Print the first three rules in the hypothesis space for testing.
%     ilasp --max-rule-length=7 -ml=7 -na min-adj/ILASP_iterative_step_1_SMALL.lp
%%#ilasp_script
%%ilasp.cdilp.initialise()
%%print(ilasp.hypothesis_to_string(list(range(1))).replace(';', ','), end='')
%%#end.

% Check the size of the hypothesis space.
%     ilasp --max-rule-length=7 -ml=7 -na min-adj/ILASP_iterative_step_1_SMALL.lp -st | clingo 0 -q | grep Models | grep -oP '\d+'
