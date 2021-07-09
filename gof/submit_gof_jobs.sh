#!/bin/bash
# Short script to write arguments.txt file for gof jobs and submit them.

MODE=$1

source gof/build_masks.sh

# Remove arguments file if it already exists as we are only appending to the file.
[[ -f gof/arguments.txt ]] && rm gof/arguments.txt

case "$MODE" in
    category)
        # Write arguments file containing all categories to be considered.
        echo "[INFO] Writing out arguments used in the submitted jobs."
        for ERA in 2016 2017 2018; do
            for CH in et mt tt em; do
                catlist=($(build_categories $CH mod-indep))
                for CAT in ${catlist[@]}; do
                    echo $ERA $CH $CAT $PWD >> gof/arguments.txt
                done
            done
        done
        ;;

    channel)
        # Write arguments file containing all categories to be considered.
        echo "[INFO] Writing out arguments used in the submitted jobs."
        for ERA in 2016 2017 2018; do
            for CH in et mt tt em; do
                echo $ERA $CH all $PWD >> gof/arguments.txt
            done
        done
        ;;

    *)
        echo "[ERROR] Given mode $MODE not known. Aborting."
        exit 1
        ;;
esac

echo "[INFO] Submit gof jobs."
condor_submit gof/submit_gof_jobs.jdl
