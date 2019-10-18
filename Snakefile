import csv
from textwrap import dedent
from pathlib import Path

CASE_LIST_PTH = 'cases.list'
BAM_MAP_PTH = '/home/mwyczalk_test/Projects/CPTAC3/CPTAC3.catalog/katmai.BamMap.dat'
GENE_GTF_PTH = '/diskmnt/Datasets/Reference/GDC/gencode.v22.annotation.gtf'
GENE_INFO_PTH = '/diskmnt/Datasets/Reference/GDC/gencode.gene.info.v22.tsv'

# Read all the cases to process
CASES = set()
with open(CASE_LIST_PTH) as f:
    CASES = set(f.read().splitlines())


# Select all the available samples of the selected cases
# sample_name -> UUID
SAMPLES = {}
# sample_name -> (case, sample_type, disease)
SAMPLE_INFO = {}
with open(BAM_MAP_PTH) as f:
    reader = csv.DictReader(f, dialect='excel-tab')
    for row in reader:
        match_case = row['case'] in CASES
        is_genomic_rna_bam = (
            row['experimental_strategy'] == 'RNA-Seq' and
            row['data_format'] == 'BAM' and
            row['reference'] == 'hg38' and
            row['data_path'].endswith('.rna_seq.genomic.gdc_realn.bam')
        )
        sample_name = row['# sample_name']
        if match_case and is_genomic_rna_bam:
            SAMPLES[sample_name] = row['UUID']
            SAMPLE_INFO[sample_name] = (row['case'], row['sample_type'], row['disease'])


rule link_gdc_rna_bams:
    """Link all available GDC RNA-seq BAMs."""
    input: local_map=BAM_MAP_PTH
    output:
        all_bams=expand('external_data/gdc_bam/{sample}.bam', sample=SAMPLES),
        all_bais=expand('external_data/gdc_bam/{sample}.bam.bai', sample=SAMPLES)
    run:
        # Create a UUID to local file path map
        all_uuids = set(SAMPLES.values())
        with open(input['local_map']) as f:
            reader = csv.DictReader(f, dialect='excel-tab')
            for row in reader:
                if not row['UUID'] in all_uuids:
                    continue
                sample = row['# sample_name']
                # Link BAM
                src_pth = Path(row['data_path'])
                dst_pth = Path('external_data/gdc_bam/{sample}.bam'.format(sample=sample))
                dst_pth.symlink_to(src_pth)
                # Link BAI
                src_pth = src_pth.with_suffix('.bam.bai')
                dst_pth = dst_pth.with_suffix('.bam.bai')
                dst_pth.symlink_to(src_pth)


rule featurecounts_stranded_readcount:
    """Readcount by featureCounts (stranded)."""
    output: count_tsv=temp('processed_data/featurecounts_stranded_readcount/{sample}.tsv')
    input: bam='external_data/gdc_bam/{sample}.bam',
           bai='external_data/gdc_bam/{sample}.bam.bai',
    log: 'logs/featurecounts_stranded/{sample}.log'
    params:
        gtf=GENE_GTF_PTH
    resources:
        io_heavy=1
    threads: 8
    shell:
        'featureCounts '
        '-g gene_id '  # feature id (-i in htseq)
        '-t exon ' # feature type (-t in htseq)
        '-T {threads} '
        '-Q 10 ' # htseq set this minimal mapping quality by default
        '-p '  # pair-end reads are considered one fragment; default HTSeq behavior
        '-B '  # both reads of a read pair need to be mapped
        '-s 2 ' # reversely stranded
        '-a {params.gtf} '
        '-o {output.count_tsv} {input.bam} 2> {log}'


rule compress_featurecounts:
    """Shrink and compress featureCounts output."""
    output: 'processed_data/featurecounts_stranded_readcount/{sample}.tsv.gz'
    input: 'processed_data/featurecounts_stranded_readcount/{sample}.tsv'
    shell: 'python shrink_featurecounts.py {input} | gzip -9 -c > {output}'


rule all_featurecounts_stranded_readcount:
    input:
        counts=expand(rules.compress_featurecounts.output[0], sample=SAMPLES)


rule generate_fpkm:
    """Generate FPKM and FPKM-UQ from the readcount."""
    output: fpkm='processed_data/readcount_and_fpkm/{sample}.tsv.gz'
    input: rc=rules.compress_featurecounts.output[0],
           gene_info=GENE_INFO_PTH
    shell: 'python gen_fpkm.py {input.gene_info} {input.rc} {output.fpkm}'


rule all_fpkms:
    input: fpkms=expand(rules.generate_fpkm.output.fpkm, sample=SAMPLES)


rule make_analysis_summary:
    input: rules.all_fpkms.input.fpkms
    output: analysis_summary='processed_data/analysis_summary.dat'
    run:
        with open(output.analysis_summary, 'w') as f:
            writer = csv.writer(f, dialect='excel-tab', lineterminator='\n')
            # Write column header
            cols = ['# case', 'disease', 'file_path', 'file_format',
                    'sample_type', 'sample_name', 'sample_uuid']
            writer.writerow(cols)

            for sample, sample_uuid in SAMPLES.items():
                case, sample_type, disease = SAMPLE_INFO[sample]
                count_tsv_pth = Path(rules.generate_fpkm.output.fpkm
                                     .format(sample=sample)).resolve(strict=True)
                writer.writerow([case, disease, count_tsv_pth.as_posix(), 'TSV',
                                 sample_type, sample, sample_uuid])
