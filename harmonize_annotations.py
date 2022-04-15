from Bio import SeqIO, Seq, motifs
from Bio.SeqRecord import SeqRecord
import re
import sys
import csv as csv2
import difflib
import subprocess
import os
import glob
import pandas as pd
import shutil
import numpy as np
import argparse
from io import StringIO
import logging
from pathlib import Path

gene_dict={
'atp6':'ATP6',
'atp8':'ATP8',
'cox1':'COX1',
'cox2':'COX2',
'cox3':'COX3',
'cob':'CYTB',
'nad1':'ND1',
'nad2':'ND2',
'nad3':'ND3',
'nad4':'ND4',
'nad4l':'ND4L',
'nad5':'ND5',
'nad6':'ND6',
'rrnS':'rRNA_12S',
'rrnL':'rRNA_16S'}   

tRNA_list=['Ala',
'Arg',
'Asn',
'Asp',
'Cys',
'Gln',
'Glu',
'Gly',
'His',
'Ile',
'Leu',
'Lys',
'Met',
'Phe',
'Pro',
'Ser',
'Thr',
'Trp',
'Tyr',
'Val']

def df2fas(df, fasta, seqName="seqName", sequence="sequence", write=True):
    """
    Dump all sequences as fasta for each markers
    """
    res = '>' + df[seqName] + '\n' + df[sequence]
    res = res.to_csv(index=False, header=False).replace('"', '')
     
    if write:
        with open(fasta, 'w') as f:
            f.write(res)
        f.close()
    return res

def fas2df(fasta, marker, source='', is_extract_tagged=False):
    """
    From a fasta file create a DataFrame [extractId, sequence] 
    Uses Biopython to parse the fasta
    """
    recGen = SeqIO.parse(fasta, "fasta")
    recId, recSeq = zip(*[(x.id , str(x.seq)) for x in recGen])
    #x.seq = class Bio.Seq.Seq, attribute str get the sequence string
    df = pd.DataFrame.from_dict(dict([('seqName', recId), ('sequence', recSeq)]))
    df['marker'] = marker
    if is_extract_tagged:
        df["extract"] = df["seqName"].str.split('__').str[1]
        df["extractSeqId"] = 1
        df["comments"] = None
        df['rawChromat'] = None
        df['cdhSeqName'] = df['seqName']
        df['id'] = None
        df['accession'] = None
    df['source'] = source
    return df


def dump_bed_n_fasta(df, path, prefix, name="all", bed_only=False, split=False):
    bedcols = ['Seq_id', 'Start', 'End', 
               'Gene_name', 'evalue', 'Direction',
               'thickStart', 'thickEnd', 'RGB']
    
    if prefix != "":
        prefix = prefix + "_"
    if split:
        try:
            os.mkdir(path+"/beds")
        except:
            pass
        for lib in df.Seq_id.unique():
            sdf = df[df.Seq_id == lib]
            sdf[bedcols]\
            .to_csv(path+"/beds/{}{}.bed".format(prefix,lib), sep='\t', index=False, header=False)

            if not bed_only:
                df2fas(sdf.drop_duplicates(subset=['Seq_id'], keep='first'),
                    path+"/genomes/{}{}.fasta".format(prefix,lib),
                    seqName="description",
                    sequence="sequence")
            
    print("writing new bed file")
    df[bedcols]\
    .to_csv(path+"/{}{}.bed".format(prefix, name), sep='\t', index=False, header=False)

    if not bed_only:
        print("writing new fasta")
        df.drop_duplicates(subset=['Seq_id'], keep='first')
        df2fas(df.drop_duplicates(subset=['Seq_id'], keep='first'),
            path +"/{}{}.fasta".format(prefix, name),
            seqName="description",
            sequence="sequence")
        
def cgi2bed_df(cgifile):
    with open(cgifile, 'r') as f:
        t = f.readlines()
        is_comment_line = False
        tsv = ("index\tGene_name\tDirection\tStart\t"+
            "End\tlength\tcodon\tSeq_id\n")
    for l in t:
        if is_comment_line:
            is_comment_line = False
            continue
        if l[0] == ">":
            Seq_id = re.sub(">(\S+).*", r"\1", l)
            is_comment_line = True
        else:
            is_comment_line = False
            line = l.strip() + '\t'+Seq_id+'\n'
            line = re.sub("c\[","-\t",line)
            line = re.sub("\[","+\t",line)
            line = re.sub("\]+","",line)
            line = re.sub(",","\t",line)
            line = re.sub("[ \t]+","\t",line)
            tsv = tsv + line
    
    arwen_df = pd.read_csv(StringIO(tsv), sep="\t", index_col=0)
    
    arwen_df['evalue'] = 0
    arwen_df = arwen_df[['Seq_id', 'Start', 'End', 'Gene_name', 
                         'evalue', 'Direction', 'codon']] 
    arwen_df.Gene_name = arwen_df.Gene_name + arwen_df.codon
    
    bedcol =  ['Seq_id', 'Start', 'End', 'Gene_name', 
               'evalue', 'Direction']
    
    return arwen_df[bedcol].dropna()
        
def gather_all_bed_file(refSeq63, refSeq89, arwen):
    bedcol =  ['Seq_id', 'Start', 'End', 'Gene_name', 'evalue', 'Direction']
    refSeq63 = pd.read_csv(refSeq63, sep="\t", header=None, names=bedcol)
    refSeq89 = pd.read_csv(refSeq89, sep="\t", header=None, names=bedcol)
    #mitoz = pd.read_csv('./Mitoz/mitoz.bed', sep="\t", header=None, names=bedcol)
    arwen = cgi2bed_df(arwen)    
    refSeq63['RGB'] = "0,255,0"
    refSeq89['RGB'] = "51,204,204"
    #mitoz['RGB'] = "255,102,153"
    arwen['RGB'] = "102,153,255"
    
    ref = pd.concat([refSeq63, refSeq89, arwen])
    ref['thickStart'] = ref['Start']
    ref['thickEnd'] = ref['End']
    ref = ref[['Seq_id', 'Start', 'End', 
               'Gene_name', 'evalue', 'Direction',
              'thickStart', 'thickEnd', 'RGB']]
    
    ref['Start'] = ref['Start'].astype('int')
    ref['End'] = ref['End'].astype('int')
    ref['thickStart'] = ref['thickStart'].astype('int')
    ref['thickEnd'] = ref['thickEnd'].astype('int')
    return ref

def consensus_63_89(bed_df):
    return bed_df.drop_duplicates(subset=['Seq_id', 'Gene_name'], keep='first')

def load_bed_n_fasta(bedfile, fastafile):
    df = pd.read_csv(bedfile, sep="\t",
                    header=None)
    df.columns = ['Seq_id', 'Start', 'End', 
                'Gene_name', 'evalue', 'Direction',
                'thickStart', 'thickEnd', 'RGB']
    
    genomes = [[g.id,g.description,g] for g in SeqIO.parse(fastafile, "fasta")]

    genomes = pd.DataFrame(genomes, columns=['Seq_id', 'description', 'genome'])
    print(df.columns)
    print(genomes.columns)
    df = df.merge(genomes, how='left', on ="Seq_id")
    
    print(df[df.genome.apply(lambda x: isinstance(x, float))])
    
    #print(df.head())
    df['length'] = df.genome.apply(lambda x: len(x.seq))
    return df

def shift_genome_and_annot_to_end_with_gene(sdf, gene_priority, lib):
    gene=""
    for g in gene_priority:
        if sdf[sdf.Gene_name == g].shape[0] != 0:
            print(lib, ": genome will be shifted to the start of", g)
            gene = g
            break
        else:
            print(lib, g, "not found")
    if gene == "":
        #print("Genome will not be shifted")
        return sdf
    else:
        
        gene_start = int(sdf[sdf.Gene_name == gene]['Start'].iloc[0])
        
        
        x = lambda r : r["genome"][gene_start:] + r["genome"][:gene_start]
        sdf['genome'] = sdf.apply(x, axis=1)
        sdf.Start = sdf.Start - gene_start
        sdf.End = sdf.End - gene_start

        sdf.Start = sdf.Start.where(sdf.Start >= 0, sdf.length + sdf.Start)
        sdf.End = sdf.End.where(sdf.End >= 0, sdf.length + sdf.End)
        
        
        print(sdf[sdf.Gene_name == gene]['Start'].iloc[0])
        print(sdf[sdf.Gene_name == gene]['End'].iloc[0])
        return sdf

def set_genome_origin(df, gene_priority_list, wdir, bedfile, fastafile):
    
    test = df.copy()
    lambda_rc = lambda seq : seq.reverse_complement()
    
    #First, get all genomes in the same direction
    
    directed_df = pd.DataFrame()
      
    for seq_id in df['Seq_id'].unique():
        #Get sub-dataframe for each genomes
        sdf = df[df['Seq_id'] == seq_id].copy()
        
        #Get the most common direction of all cox genes
        is_cox=1
        try:
            cox_dir = sdf[sdf['Gene_name'].str.contains('cox')]['Direction'].mode()[0]
    
            #If it is '-' = cox genes being on the reverse strand then we reverse everything
            if cox_dir == '-':
                print('{}: Reverse the genome to put COX genes on the forward strand'.format(seq_id))
                # 1 - Change direction marker
                sdf['OldDir'] = sdf['Direction']
                sdf['Direction'] = '+'
                sdf['Direction'] = sdf['Direction'].where(sdf['OldDir'] == '-', '-')
                sdf = sdf.drop('OldDir', axis=1)
                # 2 - reverse_complement the genome           
                sdf.genome = sdf.genome.apply(lambda_rc)
                # 3 - reverse the gene coordinate
                sdf['StartTMP'] = sdf['Start']
                sdf['thickStartTMP'] = sdf['thickStart']
                sdf['Start'] = sdf.genome.str.len() - sdf['End']
                sdf['End'] = sdf.genome.str.len() - sdf['StartTMP']
                sdf['thickStart'] = sdf.genome.str.len() - sdf['thickEnd']
                sdf['thickEnd'] = sdf.genome.str.len() - sdf['thickStartTMP']
                sdf = sdf.drop(['StartTMP', 'thickStart'], axis=1)
        except:
            print('COX genes not found in {}. The genome will not be directed'.format(seq_id))
            pass
        directed_df = pd.concat([directed_df, sdf])
    df = directed_df
    
    
    # Shift genome origin to start with the genes in gene_priority_list
    circular_df = pd.DataFrame()
    for seq_id in df[df.description.str.contains('circular')].Seq_id.unique():
        sdf = df[df['Seq_id'] == seq_id].copy()
        sdf = shift_genome_and_annot_to_end_with_gene(sdf, gene_priority_list, seq_id)
        circular_df = pd.concat([circular_df,sdf])
        
    df = pd.concat([circular_df, df[~df.description.str.contains('circular')]])

    df['sequence'] = df.genome.apply(lambda x: str(x.seq))
    
    
    dump_bed_n_fasta(df, wdir, "realigned")
    try:
        os.mkdir(wdir + "/backup/")
    except:
        pass
    
    os.rename(bedfile, wdir + "/backup/" + os.path.basename(bedfile))
    os.rename(fastafile, wdir + "/backup/" + os.path.basename(fastafile))
    return df
    
def update_annotation(gene_dir, df, bedfile, df2):
    # df2 is the annotated dataframe
    # df is the 63_89 consensus of df2
    
    udf = pd.DataFrame()
    for f in glob.glob(gene_dir+"/*.tsv"):
        gene = os.path.basename(f)[:-4]
        print(gene)
        with open(f, 'r') as ff:
            lines = ff.readlines()
            lines = [re.sub(r'[ !\-]', '', line) for line in lines]
            lines = [re.sub(r'\t+', r'\t', line) for line in lines]
            lines = [re.split(r'\t', line) for line in lines]
        sudf = pd.DataFrame.from_records(lines, columns=['Seq_id','seq_to_annotate'])        
        sudf.Seq_id = sudf.Seq_id.str.replace(">","")
        sudf.seq_to_annotate = sudf.seq_to_annotate.str.replace("[- !]","").str.strip()
        sudf['Gene_name'] = gene
        udf=pd.concat([udf, sudf])
    
    print(df.head())
    
    df = df.merge(udf, on=["Seq_id","Gene_name"], how='left')
    
    print(df[df.Gene_name=="12SrDNA_rrnS"])
    
    def get_seq_start_position(row):
        if row['seq_to_annotate'] == "":
            return -1
        
        sta = Seq.Seq(row['seq_to_annotate'])
              
        x = row['genome'].seq.find(sta)
        
        if x == -1:
            x = row['genome'].seq.find(sta.reverse_complement())
        
        if x == -1:
            print('Huho, annotated sequence for {}: {} not found in genome, keeping'.format(row['Seq_id'], row['Gene_name'])+
                  'original annotation')
        return x
    
    df['seq_to_annotate'] = df['seq_to_annotate'].fillna("")
    df['updated_start'] = df.apply(get_seq_start_position, axis=1)
    df['updated_end'] = df['seq_to_annotate'].str.len() + df['updated_start']
    
    df['Start'] = df['updated_start'].where(df['updated_start'] > -1, df['Start'])
    df['End'] = df['updated_end'].where(df['updated_end'] > 0, df['End'])
    
    df['thickStart'] = df['Start']
    df['thickEnd'] = df['End']
    df['RGB'] = df['RGB'].where(df['updated_start'] == -1, "175,255,40")
    
    df = df[df['RGB'] =="175,255,40"]
    
    bedcols = ['Seq_id', 'Start', 'End', 
               'Gene_name', 'evalue', 'Direction',
               'thickStart', 'thickEnd', 'RGB']
    
    df = pd.concat([df[bedcols], df2[bedcols]])
    
    
    dump_bed_n_fasta(df, wdir, "r_corrected", bed_only=True)
    os.rename(bedfile, wdir + "/backup/" + os.path.basename(bedfile))
    
    return df

def combine_annotations(wdir, cgifile, refSeq63, refSeq89):
    """
    Simply wrap the process of gathering all the annotations from MITOS2 and
    Arwen, and then dump everything in a bed file.
    """
    bed_tot = gather_all_bed_file(refSeq63, refSeq89, cgifile)
                        
    #beds = glob.glob(wdir+'/beds/*.bed')
    #for f in beds:
    #    os.rename(f, wdir + "/backup/beds/" + os.path.basename(f))
    #os.rename(bedfile, wdir + "/backup/" + os.path.basename(bedfile))
                        
    dump_bed_n_fasta(bed_tot, wdir, "", name="all", bed_only=True, split=True)

    return bed_tot


def blast_tRNAS(df, fastafile, wdir):
    Path(wdir + '/blast').mkdir(parents=True, exist_ok=True)
    
    fpath = Path(fastafile)
    nfpath = Path(wdir + '/blast/all.fasta')
    shutil.copy(fpath, nfpath)
    
    print(BLAST_BIN_DIR)
    
    makeblastdb_cmd = ("{} ".format(BLAST_BIN_DIR + "/makeblastdb") + 
                 "-in {} ".format(wdir + '/blast/all.fasta') +
                 "-out {} ".format(wdir + '/blast/all') +
                 "-dbtype nucl")
    
    subprocess.run(makeblastdb_cmd.split(' '))
    
    trna_df = pd.DataFrame()
    for trna in tRNA_list:
        tr = fas2df(wdir + '/../../RefSeq/refseq89_tRNA_{}.fasta'.format(trna), marker=trna, source='', is_extract_tagged=False)
        trna_df = pd.concat([trna_df, tr])
        
    trna_df['seqName'] = 'tRNA_' + trna_df['marker'] + '_' + trna_df['seqName']
    trna_df = trna_df[~trna_df['sequence'].str.contains('nformation')]
    
    df2fas(trna_df, wdir + '/blast/tRNA_query.fa')
    
    
    blast_cmd = [BLAST_BIN_DIR + '/blastn', 
                 '-db', wdir + '/blast/all', 
                 '-task', 'blastn', '-query', 
                 wdir + '/blast/tRNA_query.fa', '-outfmt', 
                 '6 qseqid sseqid sstart send sseq evalue', '-max_target_seqs', 
                 '100', '-max_hsps', '1', '-evalue', '1e-01', 
                 '-num_threads', '1', 
                 '-out', wdir + '/blast/blast_results.tsv']


    subprocess.run(blast_cmd)
    
    
    blast_res = pd.read_csv(wdir + '/blast/blast_results.tsv', sep='\t', header=None)
    blast_res.columns = ['Gene_name','Seq_id','Start','End','sseq', 'evalue']
    
    blast_res['Direction'] = "+"
    blast_res['Direction'] = blast_res['Direction'].where(blast_res['Start'] < blast_res['End'], '-')
    
    blast_res['thickStart'] = blast_res['Start'].where(blast_res['Direction'] == "+", blast_res['End'])
    blast_res['thickEnd'] = blast_res['End'].where(blast_res['Direction'] == "+", blast_res['Start'])
    blast_res['Start'] = blast_res['thickStart']
    blast_res['End'] = blast_res['thickEnd']
    
    blast_res['Start'] = blast_res['Start'] - 1
    blast_res['thickStart'] = blast_res['thickStart'] - 1
    blast_res['RGB'] = "76,0,153"
       
    df = pd.concat([df, blast_res])
    
    dump_bed_n_fasta(df, wdir, "realigned", name="all", bed_only=True, split=False)



def translate_CDS_and_control_RC(row):
    """
    Function that apply to a dataframe combining genomes and annotation.
    Will make sure that the CDS nucleotide sequences are in the correct direction,
    find the best reading frame (minimizing stop Codons), and translate it.
    Returns: Adjusted nucleotide sequence, AA sequence, reading frame from original annotation, and 
    if it is reverse complement of original annotation.
    """
    if row["Type"] == 'CDS':
        seq = row['nuc_sequence'].seq
        c = 10000
        frame = 0
        is_rc_record = False
        correct_rc = seq
        for is_rc in [False]:
            if not is_rc:
                rc = seq
            else:
                rc = seq.reverse_complement()
            for i in range(0,3):
                transl = rc[i:].translate(table = 5)
                stop = transl.count("*")
                if stop < c:
                    frame = i
                    c = stop
                    correct_rc = rc
                    is_rc_record = is_rc
        
        return SeqRecord(correct_rc[frame:]), correct_rc[frame:].translate(table = 5), frame, is_rc_record
    else:
        return row['nuc_sequence'],"" ,9, False
    
def elongate_till_stop(r):
    if r['aa_sequence'].endswith('*') or r['reading_frame'] == 9:
        return r['nuc_sequence']
    else:
        seq = r['nuc_sequence'][r['reading_frame']:]
        rest = len(seq.seq)%3
        elongation_start = r['EndF'] + 1 - rest 

        no_stop = True
        elongation = 3
        elongation_end = elongation_start + elongation
        while(no_stop):
            suffixe = r['genome'][elongation_start:elongation_end]
            if suffixe.seq.translate(table = 5).endswith('*'):

                return r['genome'][r['StartF']: elongation_start] + suffixe
            elif elongation > 120:
                print(r['Seq_id'], r['Gene_name'], "Stop not found")
                return r['nuc_sequence']
            
            elongation = elongation + 3
            elongation_end = elongation_start + elongation

def select_genes_sequences(df_annotation):
    """
    Function that apply to a dataframe combining genomes and annotation.
    Will get the sequences indicated by automatic annotations and do some sanity check.
    If stop codon absent, will look for the next one in the genome. Will also elongate in 5'
    to offer some manual verification options.
    """
    
    # 1 - reverse_complement the genome for gene in revese strand
    
    lambda_rc = lambda seq : seq.reverse_complement()
      
    df_annotation['genome'] = df_annotation['genome'].where(df_annotation['Direction'] == "+", 
                                        df_annotation.genome.apply(lambda_rc))
    
    df_annotation['StartF'] = df_annotation['Start'].where(df_annotation['Direction'] == "+", 
                                      df_annotation['genome'].str.len() - df_annotation['End'])
    df_annotation['EndF'] = df_annotation['End'].where(df_annotation['Direction'] == "+", 
                                      df_annotation['genome'].str.len() - df_annotation['Start'])
        
    x = lambda r : r["genome"][r["StartF"]:r["EndF"]]
    
    df_annotation['nuc_sequence'] = df_annotation.apply(x, axis=1)
    
    
    # TRANSLATE ORIGINAL NUC SEQUENCE AND RC IT IF NECESSARY
    
    PCG_names = 'atp8|atp6|cob|nad3|cox1|cox2|cox3|nad1|nad6|nad2|nad5|nad4|nad4l'
    df_annotation['Type'] = 'CDS'
    df_annotation['Type'] = df_annotation['Type'].where(df_annotation['Gene_name'].str.contains(PCG_names), 'nc')
    
   
    unpackdf=pd.DataFrame(df_annotation.apply(translate_CDS_and_control_RC, axis=1).tolist(),
                          columns=['nuc_sequence_2', 'aa_sequence','reading_frame', 'rev_co'], index=df_annotation.index)
    
    df_annotation=pd.concat([df_annotation,unpackdf],axis=1)
    
    df_annotation['diff_nuc_seq'] = df_annotation['nuc_sequence'].str.len() - df_annotation['nuc_sequence_2'].str.len()
    df_annotation['StartF'] = df_annotation['StartF'] + df_annotation['diff_nuc_seq']
    df_annotation.drop(labels=['nuc_sequence', 'diff_nuc_seq'], axis = 1, inplace = True)
    df_annotation.rename(columns={'nuc_sequence_2': 'nuc_sequence'}, inplace = True)
    
    
    #print(df_annotation[(df_annotation['Gene_name'] == 'cox2') & (df_annotation['Seq_id'].str.contains('a41'))]['nuc_sequence'].iloc[0].seq)
    
    
    # IF THE NUC WAS RC, THEN RC THE WHOLE GENOME FOR ELONGATION STEPS:
    df_annotation['genome'] = df_annotation['genome']\
                        .apply(lambda x : SeqRecord(x.seq.reverse_complement()))\
                        .where(df_annotation['rev_co'], df_annotation['genome'])
                    
    # I THINK I AM FORGETTING CHANGING THE GENES COORDINATE HERE X_X
    #################################################################
    #################################################################
    #################################################################
    #################################################################
    #################################################################

    
    # ELONGATE THE NUC SEQUENCE TILL A CODON STOP, IF NOT ENDS WITH STOP...
    # THIS WAS DONE BECAUSE MITOZ WAS CRAP. MAYBE USELESS WITH MITOS2
    df_annotation['is_elongated'] = df_annotation['aa_sequence'].apply(lambda x : not x.endswith('*'))
    df_annotation['nuc_sequence'] = df_annotation.apply(elongate_till_stop, axis = 1)
    
    unpackdf=pd.DataFrame(df_annotation.apply(translate_CDS_and_control_RC, axis=1).tolist(),
                          columns=['nuc_sequence', 'aa_sequence','reading_frame','rev_co'], index=df_annotation.index)
    
    df_annotation.drop(labels=['nuc_sequence', 'aa_sequence','reading_frame', 'rev_co'], axis = 1, inplace = True)
    df_annotation=pd.concat([df_annotation,unpackdf],axis=1)
    
    df_annotation['str_aa_sequence'] = df_annotation.aa_sequence.apply(lambda x: str(x))
    df_annotation['str_nuc_sequence'] = df_annotation.nuc_sequence.apply(lambda x: str(x.seq))
    
    return df_annotation

def align_genes(df_annotation, path_alignement):
    try:
        os.mkdir(path_alignement)
    except:
        pass
    
    df_annotation['aligned_nuc'] = np.NaN
    df_annotation['aligned_aa'] = np.NaN

    for gene in df_annotation.Gene_name.unique():
        
        gset = df_annotation[df_annotation.Gene_name == gene]
        
        gset['seq_length'] = (gset.Start - gset.End).abs()
        gset = gset.sort_values(by="seq_length", ascending = False).drop_duplicates(subset=['Seq_id'], keep="first")

        try:
            refset = fas2df(path_alignement + '../../../../RefSeq/refseq89_{}.fasta'.format(gene_dict[gene]), gene)[['seqName','sequence']]
            refset.columns = ['Seq_id','str_nuc_sequence']
            grset = pd.concat([gset[['Seq_id', 'str_nuc_sequence']], refset])
        except:
            grset = gset[['Seq_id', 'str_nuc_sequence']]
            pass
        
        valid_file_name = gene.replace('|','_').replace('(','_').replace(')','_')

        
        if gset.iloc[0]["Type"] != 'CDS':
            
        # If not a coding gene,
        # Write nucleotide sequence to fasta, align with muscle, 
        # then add aligned sequences in dataframe
        
            try:
                df2fas(grset, 
                    path_alignement + "{}.fasta".format(valid_file_name), 
                    seqName="Seq_id", sequence="str_nuc_sequence")
                cmd = [MUSCLE_BIN, "-in", path_alignement + "{}.fasta".format(valid_file_name),
                    "-out", path_alignement + "{}.fastaln".format(valid_file_name)]

                subprocess.run(cmd)
            except Exception as e:
                print(e)
                pass


        # If CDS align NUC along AA sequences

        if gset.iloc[0]["Type"] == 'CDS':
            
            df2fas(grset, 
                path_alignement + "{}.fasta".format(valid_file_name), 
                seqName="Seq_id", sequence="str_nuc_sequence")
            
            cmd = ["java", "-jar", MACSE_JAR, 
                "-prog", "alignSequences", "-gc_def", "5", 
                "-seq", path_alignement + "{}.fasta".format(valid_file_name),
                "-out_NT", path_alignement + "{}.fastaln".format(valid_file_name),
                "-out_AA", path_alignement + "{}_AA.fasta".format(valid_file_name)]
            
            subprocess.run(cmd)
            
        try:
            aligned_nuc = [[g.id,str(g.seq), gene] 
                for g in SeqIO.parse(path_alignement + "{}.fastaln".format(valid_file_name), "fasta")]
        except:
            aligned_nuc = [['','','']]

        df_annotation = df_annotation.merge(pd.DataFrame(aligned_nuc, 
                                                columns=['Seq_id', 'aligned_nuc_t', 'Gene_name']),
                                                how='left', on = ["Seq_id", "Gene_name"])
        
        df_annotation['aligned_nuc'] = df_annotation['aligned_nuc'].where(df_annotation['aligned_nuc'].notnull(),
                                                                         df_annotation['aligned_nuc_t'])
        df_annotation.drop(labels=['aligned_nuc_t'], axis=1, inplace=True)
        
            #df2fas(gset, 
            #    path_alignement + "{}_aa.fasta".format(gene), 
            #    seqName="Seq_id", sequence="str_aa_sequence")
        
            #cmd = ["muscle", "-in", path_alignement + "{}_aa.fasta".format(gene),
            #    "-out", path_alignement + "{}_aa.fastaln".format(gene)]
            #
            

            #aligned_aa = [[g.id,str(g.seq), gene] 
            #    for g in SeqIO.parse(path_alignement + "{}_aa.fastaln".format(gene), "fasta")]

            #df_annotation = df_annotation.merge(pd.DataFrame(aligned_aa, 
            #                                        columns=['Seq_id', 'aligned_aa_t', 'Gene_name']),
            #                                        how='left', on = ["Seq_id", "Gene_name"])
            
            #df_annotation['aligned_aa'] = df_annotation['aligned_aa'].where(df_annotation['aligned_aa'].notnull(),
            #                                                             df_annotation['aligned_aa_t'])
            #df_annotation.drop(labels=['aligned_aa_t'], axis=1, inplace=True)
        
    return df_annotation
    
def get_extra_nuc_where_gap_are_Ied(r, mol_type=""):
    
    #Check if we are working on the protein seq or the nucleotidic seq column
    if mol_type == 'missing_nuc_start':
        missing_n = r[mol_type]
    else:
        print('THIS PATH IS NOT POSSIBLE RIGHT')
        missing_n = r[mol_type] * 3
        sys.exit()

    if missing_n <= 0:
        #If no missing position, then send an empty sequence
        return SeqRecord(Seq.Seq(''))
    elif missing_n < r["StartF"]:
        # If number of missing positions is smaller than the position of the gene
        # then simply select the sequence that preceed the start of the gene
        print(r["Gene_name"], r["Seq_id"], r["StartF"], r["genome"][r["StartF"]:r["StartF"]+4].seq)
        return r["genome"][r["StartF"]-missing_n-1:r["StartF"]]
    
    else:

        # if not, then we have to loop in order to get a full sequence.
        if 'circular' in r["description"]:
            rest = missing_n - (r["StartF"])
            return r["genome"][-rest:] + r["genome"][:r["StartF"]]
        else:
            # if not circular, then its not possible to loop so we only take everything before the start of the gene
            return r["genome"][:r["StartF"]]    
        

def tabulate_fasta(path):
    with open(path, 'r') as f:
        t = f.read() 
    with open(path, 'w') as f:
        f.write(re.sub(r'(?<!^)>','\n>', re.sub(r'\n','', re.sub(r'([0-9a-z])\n',r'\1\t\t\t\t\t\t', t))))
    return None

def investigate_starting_gaps(df_annotation, path_alignement):

    df_annotation["missing_nuc_start"] = df_annotation["aligned_nuc"].str.extract(pat = '(^[\-!]+)')[0].str.len()
    df_annotation["missing_nuc_start"] = df_annotation["missing_nuc_start"].fillna(0).astype('int')
    
    #df_annotation["missing_aa_start"] = df_annotation["aligned_aa"].str.extract(pat = '(^-+)')[0].str.len()
    #df_annotation["missing_aa_start"] = df_annotation["missing_aa_start"].fillna(0).astype('int')      

    df_annotation["prefix_nuc"] = df_annotation.apply(get_extra_nuc_where_gap_are_Ied, mol_type="missing_nuc_start", axis=1)
    #df_annotation["prefix_aa"] = df_annotation.apply(get_extra_nuc_where_gap_are_Ied, mol_type="missing_aa_start", axis=1)
    
    
    def make_multiple_of_3(r):
        if len(r) >= 3:
            while len(r)%3 != 0:
                r = r[1:]
            return r
        else:
            return SeqRecord(Seq.Seq(''))
    
    df_annotation["prefix_nuc"] = df_annotation["prefix_nuc"].apply(make_multiple_of_3)
    
    #print(df_annotation[(df_annotation.Gene_name == "cox2") & (df_annotation.Seq_id.str.contains("a41"))]['prefix_nuc'].iloc[0].seq)
    #print(df_annotation[(df_annotation.Gene_name == "cox2") & (df_annotation.Seq_id.str.contains("a41"))]['str_nuc_sequence'].iloc[0])
    #print(df_annotation[(df_annotation.Gene_name == "cox2") & (df_annotation.Seq_id.str.contains("a41"))]['aligned_nuc'].iloc[0])
    #sys.exit()
    
    df_annotation["prefixed_nuc_sequence"] = df_annotation["prefix_nuc"]\
                                            .apply(lambda x : str(x.seq))\
                                            +'NNN'+ df_annotation["aligned_nuc"].str.strip('!-')
    

    #df_annotation["prefixed_aa_sequence"] = df_annotation["prefix_aa"]\
    #                                       .apply(lambda x : str(x.seq.translate(table = 5)))\
    #                                        + 'X' + df_annotation["str_aa_sequence"]
    
    for gene in df_annotation.Gene_name.unique():
        gset = df_annotation[df_annotation.Gene_name == gene].copy()
        #print(gset[gset.prefixed_nuc_sequence.str.len() > 10000].iloc[0])
        # Keep only one version of the annotation, the longest.

        gset['seq_length'] = (gset.Start - gset.End).abs()
        gset = gset.sort_values(by="seq_length", ascending = False).drop_duplicates(subset=['Seq_id'], keep="first")
        
        valid_file_name = gene.replace('|','_').replace('(','_').replace(')','_')
        
        try:
            refset = fas2df(path_alignement + '../../../../RefSeq/refseq89_{}.fasta'.format(gene_dict[gene]), gene)[['seqName','sequence']]
            refset.columns = ['Seq_id','prefixed_nuc_sequence']
            grset = pd.concat([gset[['Seq_id', 'prefixed_nuc_sequence']], refset])
        except: 
            grset = gset[['Seq_id', 'prefixed_nuc_sequence']]
            pass
        
        #After looking for potential missing nucleotides, realign (only for non coding sequences)
        if gset.iloc[0]["Type"] != 'CDS':
            grset["prefixed_nuc_sequence"] == grset["prefixed_nuc_sequence"].str.replace('-', '').str.replace('!', '')
            
            df2fas(grset, path_alignement + "{}_prefixed.fasta".format(valid_file_name), 
               seqName="Seq_id", sequence="prefixed_nuc_sequence")
        
            cmd = [MUSCLE_BIN, "-in", path_alignement + "{}_prefixed.fasta".format(valid_file_name),
              "-out", path_alignement + "{}_prefixed.tsv".format(valid_file_name)]

            subprocess.run(cmd)
            

        if gset.iloc[0]["Type"] == 'CDS':
            
            df2fas(grset, path_alignement + "{}_prefixed.fasta".format(valid_file_name), 
               seqName="Seq_id", sequence="prefixed_nuc_sequence")
            
            cmd = ["java", "-jar", MACSE_JAR, 
                "-prog", "alignSequences", "-gc_def", "5", 
                "-seq", path_alignement + "{}_prefixed.fasta".format(valid_file_name),
                "-out_NT", path_alignement + "{}_prefixed.tsv".format(valid_file_name),
                "-out_AA", path_alignement + "{}_prefixed_AA.tsv".format(valid_file_name)]
            
            subprocess.run(cmd)
            
        try:
            tabulate_fasta(path_alignement + "{}_prefixed.tsv".format(valid_file_name))
            tabulate_fasta(path_alignement + "{}_prefixed_AA.tsv".format(valid_file_name))
        except Exception as e:
            print(e)
            pass

            
    return df_annotation

    
if __name__ == "__main__":
    
    parser=argparse.ArgumentParser()

    parser.add_argument('--cmd', 
            help='''Command to execute''')

    parser.add_argument('--fasta', 
            help='''fasta path''')
    
    parser.add_argument('--bed', 
            help='''bed path''')
    
    parser.add_argument('--gene_dir', 
            help='''tsv gene dir''')
    
    parser.add_argument('--prefix_marker', 
            help='''Adding or not the artificial codon to indicated elongated sequences''')

    args=parser.parse_args()
    
    fastafile = args.fasta
    wdir = os.path.dirname(fastafile).rstrip('/')
    if wdir == "":
        wdir = "."
    
    MACSE_JAR= wdir + "/../../bin/macse_v2.05.jar"
    MUSCLE_BIN= wdir + "/../../bin/muscle"
    BLAST_BIN_DIR = wdir + "/../../bin/ncbi-blast-2.10.1+/bin"
    
    if args.cmd == "combine_annotations":
        bedfile = wdir + '/Mitos2/RefSeq63.bed'
        ref89bed = wdir + '/Mitos2/RefSeq89.bed'
        cgifile = wdir + '/Arwen/arwen.cgi'
        df = combine_annotations(wdir, cgifile, bedfile, ref89bed)  

    if args.cmd == "set_origin":
        """
        all.bed and all.fasta must be in working directory. This will get the genomes on the same strand and rotate them
        to the new point of origin.
        """
        df = load_bed_n_fasta(wdir + '/all.bed', fastafile)
        df = set_genome_origin(df, ['cox1', 'cox2', 'cox3'], wdir, wdir + '/all.bed', fastafile)
        
    if args.cmd == "blast_tRNAs":
        """
        """
        df = load_bed_n_fasta(wdir + '/r_corrected_all.bed', fastafile)
        blast_tRNAS(df, fastafile, wdir)
        
    if args.cmd == "extract_CDS":
        Path(wdir + '/workbench').mkdir(parents=True, exist_ok=True)
        df_annotation = load_bed_n_fasta(wdir + '/realigned_all.bed', fastafile)
        df_annotation = select_genes_sequences(df_annotation)
        df_annotation = align_genes(df_annotation, wdir + '/workbench/aligned/')
        
        #print("START INVESTIGATE STARTING GAPS")
        #input("Press Enter to continue...")
        
        df_annotation = investigate_starting_gaps(df_annotation, wdir + '/workbench/aligned/')

    if args.cmd == "update_annot":
        bedfile = args.bed
        df2 = load_bed_n_fasta(bedfile, fastafile)
        df = consensus_63_89(df2)
        df = update_annotation(args.gene_dir, df, bedfile, df2)
    
