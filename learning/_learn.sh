#!/usr/bin/env bash

[ "$#" -lt 1 ] && echo "Usage: $0 {lex-leader,min-adj} [$(ls */ILASP_iterative_step_1*.lp | grep -oP '\d+_(\K[^.]+)' | tr '\n' ',')]" && exit 1
CRITERIUM="$1"
[ ! -d "$CRITERIUM" ] && echo "Invalid criterium \"$CRITERIUM\", choose from {lex-leader,min-adj}." && exit 1
HYP="$([ "$#" -eq 1 ] && echo "" || echo "_$2")"
[ ! -f "$CRITERIUM/ILASP_iterative_step_1$HYP.lp" ] && echo "Hypothesis space \"$CRITERIUM/ILASP_iterative_step_1$HYP.lp\" does not exist, choose 2nd argument from {$(ls $CRITERIUM/ILASP_iterative_step_1*.lp | grep -oP '\d+_(\K[^.]+)' | tr '\n' ',')}" && exit 1

trap 'echo "[Killed]"; ps -a | grep ILASP; kill -9 $(ps -a | grep ILASP | awk "{ print \$1 }"); trap - INT; kill -s INT "$$"' INT;

SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)
cd $SCRIPT_DIR

MIN_C=3
MAX_C=4
N=1
O=1
CYCLES=2

TIMEOUT_EXAMPLE_GEN="30s"
#TIMEOUT_ILASP="40m"
TIMEOUT_ILASP="0" # (0 disables the timeout)

Symmetry_Breaking_with_ILP="$([ -z "$Symmetry_Breaking_with_ILP" ] && echo "/home/nilsk/Symmetry_Breaking_with_ILP" || echo "$Symmetry_Breaking_with_ILP")"
ILASP="$Symmetry_Breaking_with_ILP/src/ILASP4.3/ILASP"

# The options specific to the hypothesis space definitions are contained in a commented-out line respectively.
ILASP_OPTIONS="$(cat $CRITERIUM/ILASP_iterative_step_1$HYP.lp | grep -oP -- '% \K-[^ ]+( -[^ ]+)*')"

#PyLASP_SCRIPT="$Symmetry_Breaking_with_ILP/src/ilasp4.3.py"
PyLASP_SCRIPT="pylasp_script_mod.py"
# (`pylasp_script_mod.py` was produced by asking ILASP to print the version 4 PyLASP script with the command below and adding the same print statements as in the 'Symmetry_Breaking_with_ILP' repo)
#$ILASP <(echo "") --version=4 -p > pylasp_script.py 

function check_age {
    local target="$1"
    shift
    if [ ! -e "$target" ]; then
        return 0 # == true
    else
        for file in "$@"; do
            if [ "$target" -ot "$file" ]; then
                return 0 # == true
            fi
        done
    fi
    return 1 # == false
}

if check_age $CRITERIUM/examples.lasp "${BASH_SOURCE[0]}" input.lp group.py; then

    echo $'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%% Examples '"C=${C}, N=${N}, O=${O}, CYCLES=${CYCLES}"$'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' | tee $CRITERIUM/examples.lasp

    START=$(date +%s)
    for c in $(seq $MIN_C $MAX_C); do
        for n in $(seq 0 $N); do
            for o in $(seq 0 $O); do
                for cycles in $(seq 0 $CYCLES); do
                    h=$((2*c+2+n-2*cycles))
                    if (( $h >= 0 )); then
                        echo $'\n'"% C${c}H${h}N${n}O${o}" | tee -a $CRITERIUM/examples.lasp
                        timeout $TIMEOUT_EXAMPLE_GEN bash -c "clingo 0 input.lp --const 'c=${c}' --const 'h=${h}' --const 'n=${n}' --const 'o=${o}' | python group.py $CRITERIUM 2> >(tee -a $CRITERIUM/examples.lasp)"
                    fi
                done
            done
        done
    done
    END=$(date +%s)
    DIFF=$(echo "$END - $START" | bc)
    echo "% example generation took $DIFF seconds" >> $CRITERIUM/examples.lasp

fi

if check_age $CRITERIUM/hypothesis$HYP.lp "${BASH_SOURCE[0]}" ILASP_BK.lp $CRITERIUM/active_BK.lp $CRITERIUM/ILASP_iterative_step_1$HYP.lp $CRITERIUM/examples.lasp $PyLASP_SCRIPT; then

    cat ILASP_BK.lp $CRITERIUM/active_BK.lp $CRITERIUM/ILASP_iterative_step_1$HYP.lp $CRITERIUM/examples.lasp > $CRITERIUM/ILASP_input$HYP.lp
    # (fix for debug_print out-of-order when piping to tee)
    cat $PyLASP_SCRIPT | sed -e '/ilasp_script/a import sys' -e 's/^\(\s*\)\(.*print_new_iteration.*\)$/\1sys.stdout.flush(); \2; sys.stdout.flush()/g' >> $CRITERIUM/ILASP_input$HYP.lp

    echo "Starting ILASP..."
    {
        # (show number of examples)
        echo "number of positive examples: $(cat $CRITERIUM/examples.lasp | grep '#pos' | wc -l), number of negative examples: $(cat $CRITERIUM/examples.lasp | grep '#neg' | wc -l)"
        # (show size of the hypothesis space)
        #   (Note that the following only works for mode declarations: $ILASP $ILASP_OPTIONS -s $CRITERIUM/ILASP_iterative_step_1.lp | wc -l)
        echo "size of the hypothesis space: $($ILASP $ILASP_OPTIONS -st $CRITERIUM/ILASP_iterative_step_1$HYP.lp | clingo 0 --project=show -q | grep Models | grep -oP '\d+')"
        echo "ilasp options: $ILASP_OPTIONS"
        # (search for "best" hypothesis --> cover all positive examples and as many negative examples as possible)
        timeout $TIMEOUT_ILASP bash -c "$ILASP $ILASP_OPTIONS --cache-path=$CRITERIUM/ilasp-cache$HYP.bin $CRITERIUM/ILASP_input$HYP.lp | tee >(grep '^\s*:-' | tr ';' ',' > $CRITERIUM/hypothesis$HYP.lp)"
        if [ "$?" -eq "124" ]; then echo "[Timeout]"; fi
    } | tee $CRITERIUM/ilasp-run$HYP.log &

    # (monitor ILASP's memory usage)
    while ! cat $CRITERIUM/ilasp-run$HYP.log | grep "ilasp options" >/dev/null; do
        # (wait for learning to start, so that we can get the PID)
        sleep 10
    done
    # (append a line every 10 seconds)
    # TIME, PID, COMMAND, VIRT, RSS, SHR, %CPU, %MEM
    top -b -d 10 -p $(ps -u$USER | grep ILASP | awk '{ print $1 }') | awk -v OFS="," '$1=="top"{ time=$3 } $1+0>0 { print time,$1,$NF,$5,$6,$7,$9,$10; fflush() }' > $CRITERIUM/ilasp-mem$HYP.log &

    wait

fi
