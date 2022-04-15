import os
import pandas as pd
from Bio import SeqIO, Seq, motifs
from Bio.SeqRecord import SeqRecord
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
out_dir = "/Users/jamilshuvo/Documents/Thesis_Project/24_June/f"
fasta_set = [f"{out_dir}/{f}" for f in os.listdir(out_dir) if ('.fasta' in f)]
lib_taxa = pd.read_csv('./lib_taxa.tsv', sep="\t", header=None)
lib_taxa.columns = ['seqName', 'batch', 'countgroup', 'phylum', 'taxon', 'taxid']
for fasta in fasta_set:
    df = fas2df(fasta, '')
    df = df.merge(lib_taxa[['seqName', 'taxon']], how='left', on='seqName')
    df['seqName'] = df['seqName'] + '__' + df['taxon'].str.replace(' ', '_')
    df2fas(df, fasta + '_renamed')