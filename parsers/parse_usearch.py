#!/usr/bin/env python
import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import subprocess
from Bio import SeqIO
from subprocess import DEVNULL
from collections import defaultdict

def main(args):
    Rdb = parse_usearch_clustering(args.input)
    Rdb['scaffold'] = ["_".join(prot.split("_")[:-1]) for prot in Rdb['sequence']]
    Rdb['centroid_scaffold'] = ["_".join(prot.split("_")[:-1]) for prot in Rdb['centroid']]

    ## Get scaffold lengths
    s2l = {}
    seqs = {}
    for fn in os.listdir(args.fasta):
        if fn.endswith(".fa"):
            for record in SeqIO.parse(os.path.join(args.fasta,fn), "fasta"):
                s2l[record.id] = len(record.seq)
                seqs[record.id] = record.seq
    Rdb['length'] = Rdb['scaffold'].map(s2l)

    ## Get largest contig
    f_long = open(args.output + "rplF_contigs.fna", "w+")
    for cluster, db in Rdb.groupby('cluster'):
        cluster = str(cluster)
        leng = str(db.sort_values('length', ascending=False)['length'].tolist()[0])
        largest_name = str(db.sort_values('length', ascending=False)['scaffold'].tolist()[0])
        centroid = str(db.sort_values('length', ascending=False)['cluster'].tolist()[0])
        f_long.write(">" + cluster + "_" + largest_name + ":" + centroid + "." + str(leng) + "\n")
        f_long.write(str(seqs[largest_name]) + "\n")
    f_long.close()

    ## Print info
    Rdb = Rdb.rename(columns={'sequence':'gene', 'centroid':'centroid_gene', 'length':'scaffold_length'})
    Rdb.to_csv(args.output + "clustering.info.tsv", sep='\t', index=False)

def parse_usearch_clustering(loc):
    '''
    From the location of a .uc usearch file, return something like Cdb
    https://www.drive5.com/usearch/manual/cmd_calc_distmx.html
    https://www.drive5.com/usearch/manual/opt_uc.html
    '''
    dtypes = {0:'category', 1:'category', 2:np.int32, 8:'object'}
    ucols = [0,1,2,8]
    Rdb = pd.read_csv(loc, header=None, usecols=ucols,dtype=dtypes, sep='\t')
    table = defaultdict(list)

    # Find the centroids
    sdb  = Rdb[Rdb[0] == 'S']
    shdb = Rdb[Rdb[0].isin(['H', 'S'])]
    for centroid, cdb in sdb.groupby(1):
        cent = cdb[8].tolist()[0].split()[0]
        db = shdb[shdb[1] == centroid]

        for seq in db[8].tolist():
            table['cluster'].append(int(centroid))
            table['members'].append(len(db))
            table['sequence'].append(seq.split()[0])
            table['centroid'].append(cent)
    return pd.DataFrame(table)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= """Parsers UC file""")
    parser.add_argument('--input', help="path to uc file", action="store")
    parser.add_argument('--fasta', help="path to fasta files", action="store")
    parser.add_argument('--output', help="path to output files", action="store")
args = parser.parse_args()
main(args)
