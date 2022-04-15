#!/bin/bash

# To be executed in a folder containing a collection of library folders
# in which the Novoplasty.sh script has already been run.
# It will reda an assembly using Novoplasty results as seeds.
# -p number proc required; -m memory required; -n name of the node.

njobs=16
proc=4
while getopts ":n:m:p:" opt; do
  case $opt in
    j) njobs="$OPTARG"
    ;;
    p) proc="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


echo "Will process by array of $njobs jobs with $proc processors each"
wdir=$(pwd)
echo $wdir

for d in $(find "$wdir"/ -mindepth 1 -maxdepth 1 -type d -name '[A-Za-z]*[0-9]' );
    do
        libid=$(basename $d)
	echo $libid
	arr=($(grep "$libid" $wdir/lib_taxa.tsv))
        countgroup=${arr[3]}
        lineage=${dict_lineages[$countgroup]}
	echo moving into "$d/Mitogenome"
	cd $d
        mkdir -p Mitogenome
	cd Mitogenome
	
        

        sem -j $njobs get_organelle_from_reads.py -1 ${d}/${libid}.autotrim.filter1.paired_1.fq.gz -2 ${d}/${libid}.autotrim.filter1.paired_2.fq.gz -s C*.fasta -R 10 -k 105 -F animal_mt -t $proc -o ./Get_organelle
      
    done
    sem --wait

cd $wdir

