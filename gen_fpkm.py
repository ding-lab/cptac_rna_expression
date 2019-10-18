import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def readcount_to_fpkm(gene_info_df, rc_df):
    # FPKM and FPKM-UQ conversion formula from GDC
    # FPKM: <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#fpkm>
    # FPKM-UQ: <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#upper-quartile-fpkm>
    protein_coding_gene_info_df = gene_info_df.query("gene_type == 'protein_coding'")
    gene_length = gene_info_df['exon_length']
    rc_per_gene = rc_df.iloc[:, -1]
    rc_total_protein_coding = rc_per_gene[protein_coding_gene_info_df.index].sum()
    rc_quantile_75th = np.quantile(rc_per_gene[protein_coding_gene_info_df.index], 0.75)

    with np.errstate(divide='ignore'):
        shared_part = np.log(rc_per_gene) + np.log(10**9) - np.log(gene_length)
        fpkm = np.exp(shared_part - np.log(rc_total_protein_coding))
        fpkm_uq = np.exp(shared_part - np.log(rc_quantile_75th))

    fpkm_df = pd.concat([rc_per_gene, fpkm, fpkm_uq], axis=1, join='inner', sort=False)
    fpkm_df.columns = ['read_count', 'fpkm', 'fpkm_uq']
    fpkm_df.index.name = 'gene_id'
    return fpkm_df.reset_index()


def main(args):
    gene_info_df = pd.read_table(args.gene_info, index_col='gene_id')
    rc_df = pd.read_table(args.rc, comment='#', index_col='Geneid')
    fpkm_df = readcount_to_fpkm(gene_info_df, rc_df)
    fpkm_df.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Readcount to FPKM and FPKM-UQ"
    )
    parser.add_argument('gene_info', help="Path to gene annotation info")
    parser.add_argument('rc', help="Path to readcount tsv")
    parser.add_argument('out', help="Path to output tsv")
    args = parser.parse_args()
    main(args)
