#!/bin/bash

[ ! "$#" -eq 1 ] && echo "Usage: $0 {lex-leader,min-adj}" && exit 1
CRITERIUM="$1"

function gen {
    mkdir -p $1
    rm -f $1/*

    for c in $(seq 1 $2); do
        for n in $(seq 0 $3); do
            for o in $(seq 0 $4); do
                for cycles in $(seq 0 $5); do
                    local h=$((2*c+2+n-2*cycles))
                    if (( $h >= 0 )); then
                        local FILENAME="$1/C${c}H${h}N${n}O${o}.lp"
                        #echo "$cycles --> $FILENAME"
                        echo "molecular_formula(\"C\", ${c})." >> $FILENAME
                        echo "molecular_formula(\"H\", ${h})." >> $FILENAME
                        echo "molecular_formula(\"N\", ${n})." >> $FILENAME
                        echo "molecular_formula(\"O\", ${o})." >> $FILENAME
                    fi
                done
            done
        done
    done

    echo "$(($(ls -l $1 | wc -l)-2)) instances in $1"
}

#   dir   C N O cycles
gen $CRITERIUM/Gen  12 2 2      3
gen $CRITERIUM/S     7 2 2      2
