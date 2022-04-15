"""
Clément Schneider 2020-2021
"""
from Bio import SeqIO, Seq, motifs
import pandas as pd
import os
import re
import argparse
from pathlib import Path


parser=argparse.ArgumentParser()
parser.add_argument('--dirpath', 
            help='''project path''')

args=parser.parse_args()

path = args.dirpath.rstrip('/') + '/'


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

def fas2df(fasta, marker, source=''):
    """
    From a fasta file create a DataFrame [extractId, sequence] 
    Uses Biopython to parse the fasta
    """
    recGen = SeqIO.parse(fasta, "fasta")
    recId, recSeq = zip(*[(x.id , str(x.seq)) for x in recGen])
    #x.seq = class Bio.Seq.Seq, attribute str get the sequence string
    df = pd.DataFrame.from_dict(dict([('seqName', recId), ('sequence', recSeq)]))
    df['marker'] = marker
    df['source'] = source
    return df

fastas = [f for f in os.listdir(path + '03_Assemblies_Validated') if re.match(r'.+\.fa(sta)?', f)]

fdf = pd.DataFrame()

for fa in fastas:
    
    with open(path + '03_Assemblies_Validated/'+ fa, "r") as f:
        data = f.read()
        data = data.replace('*', '')
        data = data.replace('-', '')

    with open(path + 'temp.fa', 'w') as f:
        f.write(data)
        
    fdf = pd.concat([fdf, fas2df(path + 'temp.fa', 'mito', source='')])
    
Path(path + '04_Annotations').mkdir(parents=True, exist_ok=True)

df2fas(fdf, path + '04_Annotations/all.fasta', seqName="seqName", sequence="sequence", write=True)

os.remove(path + 'temp.fa')

