#!/usr/bin/env python

import sys,pandas as pd, numpy as np
from argparse import ArgumentParser

def file_to_df(infile):
    r = []
    with open(infile, 'r') as f:
        for line in f:
            if line[0]=="#": continue
            line = line.rstrip()
            line = line.rsplit()
            items = [line[0],line[3],line[6],line[7],line[11],line[13],line[19],line[20]]
            r.append(items)
    df = pd.DataFrame(r)
    df.columns = ["target","accession","seq_eval","seq_score","dom_eval","dom_score","env_from","env_to"]
    df[["seq_eval","seq_score","dom_eval","dom_score","env_from","env_to"]] = df[["seq_eval","seq_score","dom_eval","dom_score","env_from","env_to"]].apply(pd.to_numeric)
    df.index = list(range(0,len(df)))
    return df

def checkoverlap(r,overlap_frac):
    if not overlap_frac and overlap_frac!=0: return [0]
    alis = []
    store = []
    for i in range(0,len(r)): 
        f = int(r.iloc[i,6])
        t = int(r.iloc[i,7])
        this_ali = range(f,t+1)
        this_ali_len = len(this_ali)
        ## Check overlap
        if len(alis)==0: 
            store.append(i)
            alis.append(this_ali)
            continue
        overlapping = False
        for ali in alis:
            if len(ali)<this_ali_len: shortest = len(ali)
            else: shortest = this_ali_len
            o = set(this_ali).intersection(set(ali))
            overlap_len = float(len(o))/shortest
            if overlap_len>overlap_frac:
                overlapping = True
                break
                
        if not overlapping:
            alis.append(this_ali)
            store.append(i)
    return store

def parse_trusted(df,t):
    trusted = pd.read_csv(t, sep="\t", usecols=["knum","threshold","score_type"],index_col=0, dtype={"threshold":np.float64})
    ## Intersect with trusted
    ## This means only the accessions in trusted are parsed
    ta = list(set(trusted.index).intersection(set(df.accession)))
    df = pd.merge(df.loc[df.accession.isin(ta)],trusted.loc[ta],left_on="accession",right_index=True)
    parsed = df[(df.score_type == "domain") & (df.dom_score >= df.threshold) | (df.score_type == "full") & (df.seq_score >= df.threshold)]
    return parsed

def main():
    parser = ArgumentParser(description='''Parses HMMER output (--tblout and --domtblout output) and reports 1 or 
    several non-overlapping hits matching thresholds per query''')

    parser.add_argument("-i", "--infile", required=True, type=str,
            help="HMMER results file")
    parser.add_argument("-e", "--evalue", type=float, default=1e-5,
            help="E-value to use as threshold. Defaults to 1e-5.")
    parser.add_argument("-t", "--trusted", type=str,
            help="Provide trusted score cutoffs for HMMs to parse with")
    parser.add_argument("--overlap", type=float,
            help="Allowed overlapping fraction for HMM hits on the same query. Default is 0 so only best hit is stored.")

    args = parser.parse_args()
    
    ## Read results
    df = file_to_df(args.infile)

    ## Parse by trusted cutoffs and by e-value
    df = parse_trusted(df,args.trusted)
    df = df.loc[df.seq_eval<args.evalue]

    ## Count hits
    c = pd.DataFrame(df.groupby("target").count().iloc[:,0])
    c.columns = ["hit_count"]

    ## Get queries with single hits
    single = list(c[c.hit_count==1].index)

    ## Get queries with multiple hits
    multi = list(set(c[c.hit_count>1].index))

    d = df.loc[df.target.isin(single)]
    for gene in multi:
        r = df.loc[df.target==gene]
        rs = r.sort_values("dom_score",ascending=False, inplace=False)
        store = checkoverlap(rs,args.overlap)
        d = pd.concat([d,rs.iloc[store]])
    
    ## Join multiple non-overlapping hits for each query
    d.index = d.target
    d.drop("target",axis=1, inplace=True)
    d.to_csv(sys.stdout, sep="\t")
    
if __name__ == '__main__':
    main()
