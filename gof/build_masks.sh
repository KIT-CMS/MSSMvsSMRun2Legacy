#!/bin/bash

build_categories () {
    channel=$1
    mode=$2

    if [[ "$mode" == "mod-indep" ]]; then
        if [[ "$channel" == "et" || "$channel" == "mt" ]]; then
            catlist=("32" "33" "35" "36")
        elif [[ "$channel" == "tt" ]]; then
            catlist=("32" "35")
        elif [[ "$channel" == "em" ]]; then
            catlist=("2" "32" "33" "34" "35" "36" "37")
        fi
    else
        echo "[ERROR] Mode currently not supported. Aborting."
        exit 1
    fi
    echo ${catlist[@]}
}

build_masks () {
    era=$1
    channel=$2
    category=$3
    mode=$4

    tempmask=""
    finalmask=""
    # Loop over all possible masks.
    for ERA in "2016" "2017" "2018"; do
        for CH in "et" "mt" "tt" "em"; do
            catlist=($(build_categories $CH $mode))
            for CAT in "${catlist[@]}"; do
                # Get decision for this mask.
                to_mask=1
                if [[ "$era" == "combined" ]]; then
                    # In case we want to perform an inclusive test unmask all categories.
                    to_mask=0
                elif [[ "$era" == "$ERA" ]]; then
                    if [[ "$channel" == "cmb" ]]; then
                        # In case we run on a specific era, we want to unmask al channels for a combined test.
                        to_mask=0
                    elif [[ "$CAT" == "2" ]]; then
                        to_mask=0
                    elif [[ "$channel" == "$CH" ]]; then
                        # In case we want to run on a specific channel, we want to include all categories of this channel
                        if [[ "$category" == "all" ]]; then
                            to_mask=0
                        elif [[ "$category" == "$CAT" ]]; then
                            to_mask=0
                        fi
                    fi
                fi
                tempmask+="mask_htt_${CH}_${CAT}_${ERA}=${to_mask},"
            done
        done
    done
    finalmask=${tempmask:0:-1}
    echo $finalmask
}

build_masks_evaluation () {
    era=$1
    channel=$2
    category=$3
    mode=$4

    tempmask=""
    finalmask=""
    # Loop over all possible masks.
    for ERA in "2016" "2017" "2018"; do
        for CH in "et" "mt" "tt" "em"; do
            catlist=($(build_categories $CH $mode))
            for CAT in "${catlist[@]}"; do
                # Get decision for this mask.
                to_mask=1
                if [[ "$era" == "combined" ]]; then
                    # In case we want to perform an inclusive test unmask all categories.
                    to_mask=0
                elif [[ "$era" == "$ERA" ]]; then
                    if [[ "$channel" == "cmb" ]]; then
                        # In case we run on a specific era, we want to unmask al channels for a combined test.
                        to_mask=0
                    elif [[ "$channel" == "$CH" ]]; then
                        # In case we want to run on a specific channel, we want to include all categories of this channel
                        if [[ "$category" == "all" ]]; then
                            to_mask=0
                        elif [[ "$category" == "$CAT" ]]; then
                            to_mask=0
                        fi
                    fi
                fi
                tempmask+="mask_htt_${CH}_${CAT}_${ERA}=${to_mask},"
            done
        done
    done
    finalmask=${tempmask:0:-1}
    echo $finalmask
}
