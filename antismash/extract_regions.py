#!/usr/bin/env python
import sys, os
from argparse import ArgumentParser
from antismash.common.secmet import Record

def main():
	parser = ArgumentParser(description='''Parses ANTISMASH output and reports information about regions''')
	parser.add_argument("-i", "--infile", required=True, type=str, help="ANTISMASH results directory")
	parser.add_argument("-o", "--outfile", required=True, type=str, help="output directory")
	args = parser.parse_args()
    
	files = []
	dict_cluster = {}
	dict_cds = {}
	out_file = open(os.path.join(args.outfile, 'regions.txt'), "w")
    
	for i in os.listdir(args.infile):
		for j in (os.listdir(os.path.join(args.infile,i))):
			dirpath = os.path.join(args.infile,i)
			if ".region" in j:
				f = os.path.join(dirpath, j)
				files.append(f)

	for f in files:
		genome = f.split('/')[-2]
		gene = []
		for record in Record.from_genbank(f, taxon="bacteria"): 
			for region in record.get_regions():
				if record.id not in dict_cluster:
					candidates = region.candidate_clusters
					candidate_types = []
					for cand in candidates:
						candidate_types.append("{}".format(cand.kind))
						types = ",".join(candidate_types)
					dict_cluster[genome, record.id] = [str(len(record.seq)), region.get_product_string(), str(len(region.candidate_clusters)), str(types)]
				for cds in region.cds_children:
					gene.append(cds.get_name())
				if record.id not in dict_cds:
					dict_cds[genome, record.id] = gene
        
	for k, v in sorted(dict_cluster.items()):
		out_file.write(  '\t'.join(k) + '\t' +  "\t".join(v) + '\t' + ",".join(dict_cds[k])  + '\n')
	out_file.close()

if __name__ == '__main__':
	main()
