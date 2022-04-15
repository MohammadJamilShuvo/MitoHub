#!/bin/bash

# To be executed in a folder containing a collection of library folders.
# A taxonomic assignation file of the libraries must be present (lib_taxa.tsv).
# Will crawl in the folder and run Novoplasty for each libraries using a set
# of pre-selected seed.
# -p number proc required; -m memory required; -n name of the node.

njobs=64
proc=1
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

declare -A dict_lineages=( ["Collembola"]="Collembola" ["Acari"]="Acari" ["Enchytraeidae"]="Enchytraeidae" ["Nematoda"]="Nematoda" ["Aves"]="aves" ["Lumbricina"]="Enchytraeidae" ["Myriapoda"]="Myriapoda" ["Oribatida"]="Acari"  ["Gamasina"]="Acari"  ["Chilopoda"]="Chilopoda" ["Isopoda"]="Collembola" ["Diplopoda"]="Myriapoda" ["Diplura"]="Collembola" ["Lumbricidae"]="Enchytraeidae" ["Symphyla"]="Myriapoda" ["Pauropoda"]="Myriapoda" )


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
	
        
	cat << EOF > NOVO_Config.txt
Project:
-----------------------
Project name          = ${libid}_${countgroup}
Type                  = mito
    Genome Range          = 11000-25000
    K-mer                 = 32
    Max memory            = 50
    Extended log          = 0
    Save assembled reads  = no
    Seed Input            = /home/cschneider/shared_data/Configs/Seeds_Mitogenome/${lineage}.fa
    Extend seed directly  = no
    Reference sequence    = 
    Variance detection    = no
    Heteroplasmy          = 
    HP exclude list       =
    Chloroplast sequence  = 
    
    Dataset 1:
    -----------------------
    Read Length           = 151
    Insert size           = 300
    Platform              = illumina
    Single/Paired         = PE
    Combined reads        =
    Forward reads         = ${d}/${libid}.autotrim.filter1.paired_1.fq.gz
    Reverse reads         = ${d}/${libid}.autotrim.filter1.paired_2.fq.gz
    
    Optional:
    -----------------------
    Insert size auto      = yes
    Insert Range          = 1.8
    Insert Range strict   = 1.3
    Use Quality Scores    = no    
EOF

           sem -j $proc NOVOPlasty4.3.pl -c NOVO_Config.txt
      
    done
    sem --wait

cd $wdir

