% In theory, this would permit the graph symmetry breking rules from Codish et al. (see `lex.lp`) as a hypothesis.

% Run with the following ILASP options:
% --max-rule-length=10 -ml=9 -na

#modeh(1, sat(var(node),var(node),var(node))).
#modeb(1, sat(var(node),var(node),var(node))).
#modeb(3, type(var(node),var(type)), (positive)).
#modeb(2, edge(var(node),var(node)), (anti_reflexive,symmetric)).
#modeb(2, edge(var(node),var(node),var(multiplicity))).
#modeb(2, lt(var(node),var(node)), (positive)).
#modeb(1, eq(var(node),var(node)), (symmetric,negative)).
#modeb(1, lte(var(multiplicity),var(multiplicity)), (positive)).

#maxv(5).
#maxhl(1).
