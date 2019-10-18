# CPTAC RNA Expression pipeline
The pipeline generates the readcount per gene from GDC aligned stranded RNA-Seq BAMs.

- Pipeline URL: <https://github.com/ding-lab/cptac_rna_expression>
- Author: Liang-Bo Wang <liang-bo.wang@wustl.edu>


## Release
- v1.0 (2019-10-18): Initial release


## Processing description
Gene annotation in use is identical to GDC: GDC.h38 GENCODE v22 GTF `gencode.v22.annotation.gtf.gz` (md5: `291330bdcff1094bc4d5645de35e0871`), which is available on the [GDC Reference Files] page.

[GDC Reference Files]: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files


## Installation
Create a conda environment for all the dependencies:

    conda create -n cptac_expression \
        python=3.7 \
        snakemake-minimal=5.7.0 \
        samtools=1.9 htslib=1.9 \
        subread=1.6.4


## Pipeline execution
First, create a `cases.list` to contain all the CPTAC case IDs to run the
pipeline.  All samples with GDC aligned genomic RNA-Seq BAMof a listed case
will be processed.

Once the case list is available, run the following steps.  Change `$NCPUS` to
the suitable number for maximal CPU usage.

    # Link all RNA-seq BAMs
    snakemake link_gdc_rna_bams

    # Run featureCounts
    snakemake -j $NCPUS --resources io_heavy=8 -- all_featurecounts_stranded_readcount

    # Generate summary dat
    snakemake make_analysis_summary
