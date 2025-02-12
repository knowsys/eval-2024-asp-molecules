# Non-ground optimizer (NGO)

We applied the non-ground ad-hoc optimizer `ngo` to all our evaluated ASP encodings, to see whether they could be improved by simple rewriting.
Documentation regarding `ngo` can be found at <https://potassco.org/ngo/ngo.html> and <https://github.com/potassco/ngo>.

The following commands remove `#show` statements from the ASP encodings and pass them to `ngo` (using its default configuration) respectively.
The rewritten logic program is written to `<encoding>-ngo.lp` and the info printout to `<encoding>-ngo-info.txt`:

```bash
cat ../naive.lp | grep -v -E '#show' | \
    ngo --input-predicate="molecular_formula/2, element/3" --output-predicate="edge/3" \
    > naive-ngo.lp 2> naive-ngo-info.txt
cat ../naive.lp ../lex.lp | grep -v -E '#show' | \
    ngo --input-predicate="molecular_formula/2, element/3" --output-predicate="edge/3" \
    > lex-ngo.lp 2> lex-ngo-info.txt
cat ../naive-SBASS.lp | grep -v -E '#show' | \
    ngo --input-predicate="molecular_formula/2, element/3" --output-predicate="edge/3" \
    > naive-SBASS-ngo.lp 2> naive-SBASS-ngo-info.txt
cat ../smiles.lp | grep -v -E '#show' | \
    ngo --input-predicate="molecular_formula/2, element/3" --output-predicate="branching/2, symbol/2, multi_bond/2, cycle_start/2, cycle_end/2" \
    > smiles-ngo.lp 2> smiles-ngo-info.txt
```

With all the above encodings, `ngo` only performs milde (mostly small syntactic) changes, if any.
To see whether the rewritten ASP programs perform better than the original ones, test them with a medium-sized problem instance ($C_7H_{14}$, which has one cycle or multi-bond):

```bash
$ clingo 0 ../naive.lp --const c=7 --const h=14 -q | grep -P '^(Models|Time)'
Models       : 155130
Time         : 11.096s (Solving: 11.09s 1st Model: 0.00s Unsat: 0.00s)
$ clingo 0 naive-ngo.lp --const c=7 --const h=14 -q | grep -P '^(Models|Time)'
Models       : 155130
Time         : 11.249s (Solving: 11.23s 1st Model: 0.00s Unsat: 0.00s)

$ clingo 0 ../naive.lp ../lex.lp --const c=10 --const h=20 -q | grep -P '^(Models|Time)'
Models       : 52345
Time         : 9.349s (Solving: 9.31s 1st Model: 0.00s Unsat: 0.08s)
$ clingo 0 lex-ngo.lp --const c=10 --const h=20 -q | grep -P '^(Models|Time)'
Models       : 52345
Time         : 10.663s (Solving: 10.62s 1st Model: 0.00s Unsat: 0.01s)

$ clingo 0 ../smiles.lp --const c=13 --const h=26 -q | grep -P '^(Models|Time)'
Models       : 16615
Time         : 5.695s (Solving: 4.90s 1st Model: 0.03s Unsat: 0.00s)
$ clingo 0 smiles-ngo.lp --const c=13 --const h=26 -q | grep -P '^(Models|Time)'
# Produces 15 warnings about global variable in tuple of aggregate element !!! (omitted here)
Models       : 16615
Time         : 6.482s (Solving: 5.57s 1st Model: 0.10s Unsat: 0.13s)
```

You can see, that in all cases the version rewritten by `ngo` actually performs *slightly worse* than the original.
This suggests, that our encodings are reasonably well written.
