% In theory, this would permit the graph symmetry breking rules from Codish et al. (see `lex.lp`) as a hypothesis.

% Run with the ILASP option:
% --override-default-sm
% (BIAS will generate the same hypothesis space as MODE+BIAS, but take only a fraction of the time.)

#bias("
{ head(H); body(H); body(naf(H)) } = 1 :- H = sat(var__(1),var__(2),var__(3)).

{ body(B) } :- possible_body(B).
:- body(naf(B)), body(B).

body(type(var__(1),var__(t))).
body(type(var__(2),var__(t))).
body(type(var__(3),var__(t))).
possible_body(type(var__(4),var__(t))).

possible_body(edge(A,B)) :- body(type(A,_)), body(type(B,_)), A!=B.
var_multiplicity(var__(m1)). var_multiplicity(var__(m2)).
possible_body(edge(A,B,M)) :- possible_body(edge(A,B)), var_multiplicity(M).
:- body(edge(A,B)), body(edge(C,D,M)).
:- body(edge(A,B)), body(edge(B,A)).
:- body(edge(A,B,M)), body(edge(B,A,M)).
:- body(edge(A,B,M)), body(edge(C,D,N)), (A,B)<(C,D), M>N.

body_edge_cnt(CNT) :- CNT = #count{ A,B : body(edge(A,B)); A,B : body(edge(A,B,M)) }.
:- body_edge_cnt(CNT), CNT < 1.
:- body_edge_cnt(CNT), CNT > 2.

possible_body(lte(A,B)) :- body(edge(_,_,A)), body(edge(_,_,B)), A!=B.
:- body(lte(A,B)), body(lte(B,A)).

possible_body(lt(A,B)) :- body(type(A,_)), body(type(B,_)), A!=B.
possible_body(naf(eq(A,B))) :- body(type(A,_)), body(type(B,_)), A<B.
:- body(lt(A,B)), body(naf(eq(A,B))).
:- body(lt(A,B)), body(naf(eq(B,A))).
:- body(lt(A,B)), body(lt(B,A)).

body_comp_cnt(CNT) :- CNT = #count{ A,B : body(lt(A,B)); A,B : body(lte(A,B)); A,B : body(eq(A,B)); A,B : body(naf(eq(A,B))) }.
:- body_comp_cnt(CNT), CNT < 2.
:- body_comp_cnt(CNT), CNT > 3.

:- #count{ B : body(B) } < 7.
:- #count{ B : body(B) } > 9.

#show body/1.
#show head/1.

% To see whether the Codish conditions are contained, exclude all other rules.
%     ilasp --override-default-sm min-adj/ILASP_iterative_step_1_BIAS.lp -st | clingo 0
%%rule1 :- head(sat(I,K,J)), #count{ B : body(B) } == 9, body(type(I,T)), body(type(J,T)), body(type(K,T)), body(type(L,T)), body(lt(I,J)), body(edge(I,K)), body(edge(J,L)), body(lt(L,K)), body(naf(eq(I,L))).
%%rule2 :- head(sat(I,K,J)), #count{ B : body(B) } == 7, body(type(I,T)), body(type(J,T)), body(type(K,T)), body(lt(I,J)), body(edge(I,K,N)), body(edge(J,K,M)), body(lte(N,M)).
%%rule3 :- #count{ B : body(B) } == 7, body(type(I,T)), body(type(J,T)), body(type(K,T)), body(edge(I,K)), body(lt(I,J)), body(naf(sat(I,K,J))), body(naf(eq(K,J))).
%%:- not rule1, not rule2, not rule3.
").

% Print the first three rules in the hypothesis space for testing.
%     ilasp --override-default-sm min-adj/ILASP_iterative_step_1_BIAS.lp
%%#ilasp_script
%%ilasp.cdilp.initialise()
%%print(ilasp.hypothesis_to_string(list(range(3))).replace(';', ','), end='')
%%#end.

% Check the size of the hypothesis space.
%     ilasp --override-default-sm min-adj/ILASP_iterative_step_1_BIAS.lp -st | clingo 0 -q | grep Models | grep -oP '\d+'
