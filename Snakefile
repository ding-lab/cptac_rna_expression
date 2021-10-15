import csv
from collections import namedtuple
from textwrap import dedent
from pathlib import Path
import logging

# List of sample names to run the pipeline
SAMPLE_LIST_PTH = config['sample_list']
# The mapping of sample name to other information.
BAM_MAP_PTH = config['bam_map']
GENE_GTF_PTH = config['gdc_gtf']  # gencode.v22.annotation.gtf
GENE_INFO_PTH = config['gdc_gene_info']  # gencode.gene.info.v22.tsv
WORKFLOW_ROOT = config['workflow_root']  # Path to this repository
# A [src_part, dst_part] to replace the file path of the BAM path
# Useful to change the BAM map in Docker
REPLACE_BAM_PTH = config.get('replace_bam_path', None)

_logger = logging.getLogger(__name__)

# Read all the cases to process
with open(SAMPLE_LIST_PTH) as f:
    SAMPLES = f.read().splitlines()
if len(SAMPLES) != len(set(SAMPLES)):
    _logger.error('There are duplicated samples in the sample list!')
SAMPLES = set(SAMPLES)


# Select all the available samples of the selected cases.
# SAMPLE_INFO is a mapping of:
#   sample_name -> SampleInfo(case, sample_type, disease, UUID, BAM_src_path)
SampleInfo = namedtuple('SampleInfo', 'case, disease, uuid, bam_pth')
SAMPLE_INFO = {}
with open(BAM_MAP_PTH) as f:
    reader = csv.DictReader(f, dialect='excel-tab')
    for row in reader:
        sample_name = row['# sample_name']
        if sample_name not in SAMPLES:
            continue
        is_genomic_rna_bam = (
            row['experimental_strategy'] == 'RNA-Seq' and
            row['data_format'] == 'BAM' and
            row['reference'] == 'hg38' and
            row['data_path'].endswith('.rna_seq.genomic.gdc_realn.bam')
        )
        if not is_genomic_rna_bam:
            _logger.warning(f'{sample_name} is not a genomic RNA-Seq BAM')

        bam_pth = Path(row['data_path'])
        SAMPLE_INFO[sample_name] = SampleInfo(
            row['case'], row['disease'], row['UUID'], bam_pth
        )


def find_sample_bam_path(wildcards, replace=REPLACE_BAM_PTH):
    bam_pth = str(SAMPLE_INFO[wildcards.sample].bam_pth)
    if replace is not None:
        src_part, dst_part = replace
        bam_pth = bam_pth.replace(src_part, dst_part)
    return {
        'bam': bam_pth,
        'bai': bam_pth + '.bai',
    }


rule featurecounts_stranded_readcount:
    """Readcount by featureCounts (stranded)."""
    output: count_tsv=temp('featurecounts_stranded_readcount/{sample}.tsv')
    input: unpack(find_sample_bam_path)
    log: 'logs/featurecounts_stranded/{sample}.log'
    params:
        gtf=GENE_GTF_PTH
    resources:
        io_heavy=1,
        mem_mb=lambda wildcards, attempt: 16000 + 16000 * (attempt - 1)
    threads: 16
    group: "featurecounts"
    shell:
        'featureCounts '
        '-g gene_id '  # feature id (-i in htseq)
        '-t exon '  # feature type (-t in htseq)
        '-T {threads} '
        '-Q 10 '  # htseq set this minimal mapping quality by default
        '-p '  # pair-end reads are considered one fragment; default HTSeq behavior
        '-B '  # both reads of a read pair need to be mapped
        '-s 2 '  # reversely stranded
        '-a {params.gtf} '
        '-o {output.count_tsv} {input.bam} 2> {log}'


rule compress_featurecounts:
    """Shrink and compress featureCounts output."""
    output: 'featurecounts_stranded_readcount/{sample}.tsv.gz'
    input: rules.featurecounts_stranded_readcount.output.count_tsv
    threads: 2
    group: "featurecounts"
    shell: 'python {WORKFLOW_ROOT}/shrink_featurecounts.py {input} | gzip -9 -c > {output}'


rule generate_fpkm:
    """Generate FPKM and FPKM-UQ from the readcount."""
    output: fpkm='readcount_and_fpkm/{sample}.tsv.gz'
    input: rc=rules.compress_featurecounts.output[0],
           gene_info=GENE_INFO_PTH
    shell: 'python {WORKFLOW_ROOT}/gen_fpkm.py {input.gene_info} {input.rc} {output.fpkm}'


rule all_featurecounts_stranded_readcount:
    input:
        counts=expand(rules.compress_featurecounts.output[0], sample=SAMPLES)


rule all_fpkms:
    input: fpkms=expand(rules.generate_fpkm.output.fpkm, sample=SAMPLES)


rule make_analysis_summary:
    input: rules.all_fpkms.input.fpkms
    output: analysis_summary='analysis_summary.dat'
    run:
        with open(output.analysis_summary, 'w') as f:
            writer = csv.writer(f, dialect='excel-tab', lineterminator='\n')
            # Write column header
            cols = ['# run_name', 'case', 'disease',
                    'data_path', 'file_format',
                    'sample_name', 'sample_uuid']
            writer.writerow(cols)

            for sample, info in SAMPLE_INFO.items():
                count_tsv_pth = Path(
                    rules.generate_fpkm.output.fpkm.format(sample=sample)
                ).resolve(strict=True)
                writer.writerow([
                    # sample name is unique to be run name
                    sample, info.case, info.disease,
                    str(count_tsv_pth), 'TSV',
                    sample, info.uuid
                ])
