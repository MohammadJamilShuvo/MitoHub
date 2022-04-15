mkdir -p mitogenomes/01_Assemblers_output

for d in $(find ./ -maxdepth 1 -type d -iname "[A-Za-z]*[0-9]"); do 

    cp $d/Mitogenome/C*_1_*.fasta ./mitogenomes/01_Assemblers_output/$(basename $d)__Novo.fasta;
    
    # Conveying the information weither the contig is circular or not
    for f in $(find $d/Mitogenome -maxdepth 1 -type f -iname "C*_1_*.fasta"); do
	    echo $( basename $f)
	    if [[ $( basename $f) == "Circular"* ]]; then
		 echo "CIRCULAIRE"
		 sed -E -i 's/>(.+)/>\1__circular/' ./mitogenomes/01_Assemblers_output/$(basename $d)__Novo.fasta
            fi
    done
    
    cp $d/Mitogenome/Get_organelle/animal_mt*graph1.1*.fasta ./mitogenomes/01_Assemblers_output/$(basename $d)__GO.fasta;
    cp $d/Mitogenome/Get_organelle/animal_mt*graph1*.gfa ./mitogenomes/01_Assemblers_output/$(basename $d)__GO.gfa;

    cp $d/Mitogenome/Get_organelle/animal_mt*graph2*.fasta ./mitogenomes/01_Assemblers_output/$(basename $d)__GO2.fasta;
    cp $d/Mitogenome/Get_organelle/animal_mt*graph2*.gfa ./mitogenomes/01_Assemblers_output/$(basename $d)__GO2.gfa;

done

cp ./lib_taxa.tsv ./mitogenomes/lib_taxa.tsv
