#ilasp_script
ilasp.cdilp.initialise()
solve_result = ilasp.cdilp.solve()

c_egs = None
if solve_result is not None:
  c_egs = ilasp.find_all_counterexamples(solve_result)

conflict_analysis_strategy = {
  'positive-strategy': 'all-ufs',
  'negative-strategy': 'single-as',
  'brave-strategy':    'all-ufs',
  'cautious-strategy': 'single-as-pair'
}

while c_egs and solve_result is not None:
  ce = ilasp.get_example(c_egs[0]['id'])
  constraint = ilasp.cdilp.analyse_conflict(solve_result['hypothesis'], ce['id'], conflict_analysis_strategy)

  # An example with recorded penalty of 0 is in reality an example with an
  # infinite penalty, meaning that it must be covered. Constraint propagation is,
  # therefore, unnecessary.
  if not ce['penalty'] == 0:
    c_eg_ids = list(map(lambda x: x['id'], c_egs))
    prop_egs = []
    if ce['type'] == 'positive':
      prop_egs = ilasp.cdilp.propagate_constraint(constraint, c_eg_ids, {'select-examples': ['positive'], 'strategy': 'cdpi-implies-constraint'})
    elif ce['type'] == 'negative':
      prop_egs = ilasp.cdilp.propagate_constraint(constraint, c_eg_ids, {'select-examples': ['negative'], 'strategy': 'neg-constraint-implies-cdpi'})
    elif ce['type'] == 'brave-order':
      prop_egs = ilasp.cdilp.propagate_constraint(constraint, c_eg_ids, {'select-examples': ['brave-order'],    'strategy': 'cdoe-implies-constraint'})
    else:
      prop_egs = [ce['id']]

    ilasp.cdilp.add_coverage_constraint(constraint, prop_egs)

  else:
    ilasp.cdilp.add_coverage_constraint(constraint, [ce['id']])

  solve_result = ilasp.cdilp.solve()

  if solve_result is not None:
    c_egs = ilasp.find_all_counterexamples(solve_result)


if solve_result:
  print(ilasp.hypothesis_to_string(solve_result['hypothesis']))
else:
  print('UNSATISFIABLE')

ilasp.stats.print_timings()

#end.
