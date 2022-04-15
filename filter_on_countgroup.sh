#!/bin/bash

# Select some libs based on a match in lib_taxa. Could be species name or typically countgroup.
# Go into the directory were you want to select things
# arg1 = keyword, arg2 = lib_taxa file (path)
# Results go in a tmp directory created in working dir

mkdir -p TMP_select_$1

grep "$1" $2 |while IFS=$'\t', read -r a b c d e; do 
cp -r *"$a"* TMP_select_$1; 
done
