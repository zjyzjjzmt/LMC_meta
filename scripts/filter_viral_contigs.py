import pandas as pd
from Bio import SeqIO
import argparse

def main():
    parser = argparse.ArgumentParser(description="Filter viral contigs based on Methods")
    parser.add_argument("--vs2", help="VirSorter2 score file")
    parser.add_argument("--vf", help="VirFinder score file")
    parser.add_argument("--genomad", help="geNomad summary file")
    parser.add_argument("--fasta", help="Input contigs fasta")
    parser.add_argument("--out", help="Output filtered fasta")
    args = parser.parse_args()

    # 1. Load VirSorter2 (viral_gene > 0 OR ...)
    vs2 = pd.read_csv(args.vs2, sep="\t")
    # Methods: viral_gene > 0 OR (viral_gene=0 AND host_gene=0 OR score >= 0.95 OR hallmark > 2)
    # Note: Adjust column names based on actual VS2 output
    keep_vs2 = set(vs2[
        (vs2['viral_genes'] > 0) | 
        ((vs2['viral_genes'] == 0) & (vs2['host_genes'] == 0)) |
        (vs2['max_score'] >= 0.95) |
        (vs2['hallmark_genes'] > 2)
    ]['seqname'])

    # 2. Load VirFinder (score >= 0.9 AND p <= 0.01)
    vf = pd.read_csv(args.vf, sep="\t")
    keep_vf = set(vf[(vf['score'] >= 0.9) & (vf['pvalue'] <= 0.01)]['name'])

    # 3. Load geNomad (Default high confidence, usually score > 0.7 or specific flag)
    genomad = pd.read_csv(args.genomad, sep="\t")
    keep_genomad = set(genomad['seq_name']) # Assuming input file is already filtered summary

    # Union of all valid IDs
    valid_ids = keep_vs2.union(keep_vf).union(keep_genomad)

    # Write output
    with open(args.out, "w") as f_out:
        for record in SeqIO.parse(args.fasta, "fasta"):
            if record.id in valid_ids:
                SeqIO.write(record, f_out, "fasta")

if __name__ == "__main__":
    main()