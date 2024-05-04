% A version of this file without comments was produced with the following command:
% $> sed '/#show/d; /^[ \t]*%/d; s/[ \t]*%.*$//' smiles-to-edge.lp | cat -s > smiles-to-edge_min.lp

% TO EDGE
% (required by CHECK, SUBFORMULA, SUBSTRUCTURE)

% Derive a labeled edge relation from the SMILES model.
% - bonds due to tree structure
edge(I1, I2, 1) :- parent(I1, _, I2), not multi_bond(I2, _), I1 < I2.
edge(I1, I2, MULTIPLICITY) :- parent(I1, _, I2), multi_bond(I2, MULTIPLICITY), I1 < I2.
% - bonds due to cycle markers
edge(I1, I2, 3) :- I1 < I2,
                   cycle_start(I1, CYCLE_1), cycle_end(I2, CYCLE_1),
                   cycle_start(I1, CYCLE_2), cycle_end(I2, CYCLE_2), CYCLE_1 < CYCLE_2,
                   cycle_start(I1, CYCLE_3), cycle_end(I2, CYCLE_3), CYCLE_2 < CYCLE_3.
edge(I1, I2, 2) :- not edge(I1, I2, 3), I1 < I2,
                   cycle_start(I1, CYCLE_1), cycle_end(I2, CYCLE_1),
                   cycle_start(I1, CYCLE_2), cycle_end(I2, CYCLE_2), CYCLE_1 < CYCLE_2.
edge(I1, I2, 1) :- not edge(I1, I2, 3), not edge(I1, I2, 2), I1 < I2,
                   cycle_start(I1, CYCLE), cycle_end(I2, CYCLE).

#show edge/3.
