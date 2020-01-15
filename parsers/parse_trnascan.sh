#!/bin/bash

if [ "$#" -lt 1 ]
then
  echo "usage: parse_trnascan.sh <tRNAscan output folder>"
  echo "tRNAscan output folder should contain the output and stats files"
  echo "output files should be named with sample names and stats files as samplename_stats" 
exit 1
fi

# tRNAscan output folder
trna_out=$1

# Print headers
printf '%s\t%s\t%s\t%s\n' "Bin_ID" "standard_AA_count" "AA_count" "Unique"

# Loop through samples and compute stats
for sample in $(find ${trna_out} -mindepth 1 | grep -vE 'stats|fa'); do
  trna_stats=${sample}_stats
  # Extract the total tRNAs decoding standard AA from stats file
  std_AA_count=$(cat ${trna_stats} | grep 'tRNAs decoding Standard 20 AA:' | cut -d : -f2 | tr -d " ")
  # Make sure we are counting the AA's correctly
  AA_count=$(cat ${sample} | tail -n+4 | grep -vE 'pseudo|Undet|SeC|Sup' | cut -f5 | sort | wc -l)
  # Count instances of unique tRNA's
  unique_trna=$(cat ${sample} | tail -n+4 | grep -vE 'pseudo|Undet|fMet|SeC|Sup' | cut -f5 | sed 's/Ile2/Ile/g' | sed 's/iMet/Met/g' | sort | uniq | wc -l) 
  # Print stats
  printf '%s\t%s\t%s\t%s\n' $(basename ${sample}) ${std_AA_count} ${AA_count} ${unique_trna}
done
