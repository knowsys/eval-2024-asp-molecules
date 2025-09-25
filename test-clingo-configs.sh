#!/usr/bin/env bash

# scp naive.lp lex.lp smiles_eval.lp min_main_chain_len.py test-clingo-configs.sh wille:/home/nilsk/genmol_clingo_config
# nohup nix develop --command ./test-clingo-configs.sh

DIR='results_clingo_config'
mkdir -p $DIR

cat naive.lp lex.lp > naive_lex.lp

MAX_C=15

TIMEOUT='1m'

for encoding in "naive.lp" "naive_lex.lp" "smiles_eval.lp"
do
    for configuration in "auto" "frumpy" "jumpy" "tweety" "handy" "crafty" "trendy" "many"
    do
        for c in $(seq 1 $MAX_C)
        do
            ((h=2*c))
            echo "${encoding} ${configuration} C_${c} H_${h} O"
            output_file="${DIR}/${encoding}_${configuration}.csv"
            echo -n "$c," >> ${output_file}
            ((natoms=c+1))
            m=$(python min_main_chain_len.py $natoms)
            timeout --foreground $TIMEOUT clingo 0 ${encoding} --const c=${c} --const h=${h} --const o=1 --const m=${m} -q --configuration=${configuration} > ${DIR}/tmp
            if [ $? -eq 124 ]; then
                echo "TIMEOUT=${TIMEOUT}" | tee -a "${output_file}"
                break
            else
                cat ${DIR}/tmp | grep -P '^Time' | grep -oP '[0-9]+\.[0-9]+' | head -1 >> ${output_file}
            fi
        done
    done
done
