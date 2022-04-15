"""
Clément Schneider 2020-2021
"""
from Bio import SeqIO, Seq, motifs
import re
import sys
import csv as csv2
import difflib
import subprocess
import logging
import os
import argparse

parser=argparse.ArgumentParser()

parser.add_argument('--dirpath', 
            help='''project path''')

args=parser.parse_args()


path = args.dirpath.rstrip('/') + '/'

assemblies_path = path + '01_Assemblers_output/'
preproc_out_path = path + '02_Assemblies_Preprocessed/'

print("assemblies path=", assemblies_path)

logging.basicConfig(filename=path + 'Preprocessing.log',
                    level=logging.INFO)

logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

#ch = logging.StreamHandler()
#ch.setLevel(logging.INFO)
#logging.addHandler(ch)

def harmonize_strand(lib, path):
    """
    lib: Sequenced library unique identifier.
    path: path to the fasta files directory.
    """
    print(path+'{}__Novo.fasta'.format(lib))
    # Open the two resultings fasta, Novoplasty and Get_Organelle. 
    # If one is missing, go for manual check.
    missing = list()
    try:
        novo = list(SeqIO.parse(path+'{}__Novo.fasta'.format(lib), "fasta"))
    except Exception as e:
        print(e)
        missing.append("Novoplasty")
    try:
        get_org = list(SeqIO.parse(path+'{}__GO.fasta'.format(lib), "fasta"))
    except:
        missing.append("Get_Organelle")
    if len(missing) > 0:
        logging.warning("{}: result(s) from {} missing".format(lib, " and ".join(missing)))
        return(None, None)
    
    # If the files contain several contigs, go for manual check.
    if len(novo) >= 2 or len(get_org) >= 2:
        logging.warning("{}: more than one contig in fasta files !!! --> manual check".format(lib))
        return(None, None)
    
    # Make sure that results from Novoplasty and Get_Organelle are the same strand
    novo_seq = novo[0].seq
    get_org_seq = get_org[0].seq
    i = 0
    stop = len(novo) - 50
    check = False
    while not check:
        check = novo_seq[i:i+50] in get_org_seq
        i = i + 50
        if i >= stop:
            break
    if not check:
        i=0
        print(lib, 'Get_Organelle get reverse complemented')
        get_org_seq = get_org_seq.reverse_complement()
        while not check:
            check = novo_seq[i:i+50] in get_org_seq
            i = i + 50
            if i >= stop:
                break
    
    # If Novoplasty and Get_Organelle contigs cannot be matched at all, go for manual check
    # (should not be possible)
    if not check:
        logging.warning("{}: Novoplasty and GetOrganelle contigs do not match !!! --> manual check".format(lib))
        return None, None
    
    # Get the correct strand for Get_Organelle and return both results.
    get_org[0].seq = get_org_seq
    return novo[0], get_org[0]

def approx_align(seq1, seq2):
    """
    Shift the sequence seq1 until it is aligned with seq2 on a sub-sequence of 50 nucleotides
    picked arbitrarily. seq1 and seq2 are from a circular genome and are supposed to be
    nearly identical.
    """
    s = 0
    found = list()
    while len(found) == 0 and s <= len(seq2.seq):
        motif = str(seq2.seq[s:s+50])
        found = list(re.finditer(motif,str(seq1.seq)))
        s = s+50
    
    # If no motif found, then problem. Should not be possible in pipeline.
    if len(found) == 0:
        logging.warning("{}: quick alignment was not possible".format(lib))
        return None, None
    
    for i in re.finditer(motif,str(seq1.seq)):
        start  = i.start()
    seq1 = seq1[start:].seq + seq1[:start]
    
    return seq1, seq2

##### BEGIN SCRIPT #####

try:
    os.mkdir(preproc_out_path)
except:
    pass
    
#Get list of libraries in the batch (specified in the lib_taxa.tsv
with open(path+'lib_taxa.tsv') as f:
    sp_list = csv2.reader(f, delimiter='\t')
    libs = [lib[0] for lib in sp_list]


# For each libs, run the harmonize_strand. If contigs from Novo and GO could be matched
# then continue
for lib in libs:
    novo, get_org = harmonize_strand(lib, assemblies_path)
    if not novo:
        continue
        
    novo, get_org = approx_align(novo, get_org)
    if not novo:
        continue
    # If sequences are identical, then write only one copy in a new fasta file.
    # If sequences not identical, write the two in the new fasta file, and run a Muscle alignement to help
    # spot check the differences.
    # Topology is reported in the file name. I believe Get_Organelle to be more trustworthy to
    # deal with repeated patterns and avoid ambiguous closing of the sequence.
    # Therefore, if Get_Organelle say circular, I accept it. If Get_Organel says linear
    # agains Novoplasty, it need to be checked.
    # Sometimes, Novoplasty report circularity but missed repeats that Get_Organelle could solve !

    
    if novo.seq == get_org.seq:
        logging.warning('{}: Novoplasty and GetOrganelle sequences are identical'.format(lib))
        if "circular" in get_org.id:
            logging.warning('{}: at least GetOrganelle reported circular topology, it should be fine'.format(lib))
            topology = "circular"
        elif ("circular" in novo.id) and ("linear" in get_org.id):
            logging.warning('{}: Only Novoplasty reported circularity, it not usual. Maybe check GO GFA ?'.format(lib))
            topology = "ambiguous"
        else:
            logging.warning('{}: topology is linear'.format(lib))
            topology = "linear"
            
        get_org.id=lib
        get_org.description = "organism={} assembler=Get_Organelle len={} topology={}"\
                            .format("species", len(get_org.seq), topology)
            
        with open(preproc_out_path+"{}__{}_mtGenome.fasta".format(lib, topology), "w") as output_handle:
            SeqIO.write(get_org, output_handle, "fasta")
            
    else:
        logging.warning('{}: Novoplasty and GetOrganelle sequences are not identical, they will be aligned with Muscle for comparison and differences are will be reported'.format(lib))
        if "circular" in novo.id:
            topology="circular".format(len(novo.seq))
        else:
            topology="linear".format(len(novo.seq))
            
        novo.id=lib
        novo.description = "organism={} assembler=Novoplasty len={} topology={}"\
                            .format("species", len(get_org.seq), topology)

        if "circular" in get_org.id:
            topology="circular".format(len(get_org.seq))
        else:
            topology="linear".format(len(get_org.seq))
        get_org.id=lib
        get_org.description = "organism={} assembler=Get_Organelle len={} topology={}"\
                            .format("species", len(get_org.seq), topology)
        
        with open(preproc_out_path+"to_check__{}_mtGenome.fasta".format(lib), "w") as output_handle:
            SeqIO.write(get_org, output_handle, "fasta")
            SeqIO.write(novo, output_handle, "fasta")
        
        cmd = [path+"../bin/muscle", "-in", preproc_out_path + "to_check__{}_mtGenome.fasta".format(lib),
              "-out", preproc_out_path+"to_check__{}_mtGenome.fastaln".format(lib)]
        subprocess.run(cmd)
        
        os.remove(preproc_out_path + "to_check__{}_mtGenome.fasta".format(lib))
        os.rename(preproc_out_path+"to_check__{}_mtGenome.fastaln".format(lib), 
                  preproc_out_path+"to_check__{}_mtGenome.fasta".format(lib))
        
        # Report differences in aligned sequences.
        aligned = list(SeqIO.parse(preproc_out_path+"to_check__{}_mtGenome.fasta".format(lib), "fasta"))
        for i in range(len(aligned[0].seq)):
            if aligned[0].seq[i] != aligned[1].seq[i]:
                if aligned[0].seq[i] != '-' or aligned[1].seq[i] != '-':
                    logging.info(lib + " l"+ str(int(i / 60)+1) + " c"+ str(i - (int(i / 60) * 60) + 1) + " " + aligned[0].seq[i] + " " + aligned[1].seq[i])
                    
logging.warning('Preprocessing over. Find this log again in the "preprocessing.log" file.'.format(lib))
logging.warning('In the 02_Assemblies_Preprocessed, files with "circular" in name should be allright. The files with "to_check" in name means that one of the two sequences inside must be selected'.format(lib))
        
