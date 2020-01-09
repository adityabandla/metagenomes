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
            items = [line[0],line[2],line[3],line[12],line[15],line[16],line[19],line[20]]
            r.append(items)
    df = pd.DataFrame(r)
    df.columns = ["target","target_len","query_seq","dom_eval","hmm_from","hmm_to","env_from","env_to"]
    df[["target_len","dom_eval","hmm_from","hmm_to","env_from","env_to"]] = df[["target_len","dom_eval","hmm_from","hmm_to","env_from","env_to"]].apply(pd.to_numeric)
    df["target"] = df["target"].str.replace(r'.hmm', '')
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

def main():
    parser = ArgumentParser(description='''Parses HMMER output (--tblout and --domtblout output) and reports 1 or 
    several non-overlapping hits matching thresholds per query''')

    parser.add_argument("-i", "--infile", required=True, type=str,
            help="HMMER results file")
    parser.add_argument("-e", "--evalue", type=float, default=1e-15,
            help="E-value to use as threshold. Defaults to 1e-15.")
    parser.add_argument("--overlap", type=float, default=0,
            help="Allowed overlapping fraction for HMM hits on the same query. Default is 0 so only best hit is stored.")
    parser.add_argument("--coverage", type=float, default=0.35,
            help="Minimum hmm coverage. Default is 0.35")

    args = parser.parse_args()
    
    ## Read results
    df = file_to_df(args.infile)

    ## Count hits
    c = pd.DataFrame(df.groupby("query_seq").count().iloc[:,0])
    c.columns = ["hit_count"]

    ## Get queries with single hits
    single = list(c[c.hit_count==1].index)

    ## Get queries with multiple hits
    multi = list(set(c[c.hit_count>1].index))

    d = df.loc[df.query_seq.isin(single)]
    for gene in multi:
        r = df.loc[df.query_seq==gene]
        rs = r.sort_values("dom_eval",ascending=True, inplace=False)
        store = checkoverlap(rs,args.overlap)
        d = pd.concat([d,rs.iloc[store]])
        
    ## Calculate coverage
    d["coverage"] = (d["hmm_to"]-d["hmm_from"])/d["target_len"]
    ## Parse by e-value & coverage
    d = d[(d["coverage"]>=args.coverage) & (d["dom_eval"]<=args.evalue)]
    d.index = d.target
    d.drop("target", axis=1, inplace=True)    
    d.to_csv(sys.stdout, sep="\t")
    
if __name__ == '__main__':
    main()
