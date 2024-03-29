# CPTAC RNA Expression pipeline
The pipeline generates the readcount per gene from GDC aligned stranded RNA-Seq BAMs.

- Pipeline URL: <https://github.com/ding-lab/cptac_rna_expression>
- Author: Liang-Bo Wang <liang-bo.wang@wustl.edu>


## Release
- v1.0 (2019-10-18): Initial release


## Processing description
The readcount is generated by featureCounts (subread v1.6.4) under stranded mode with parameters: `-g gene_id -t exon -Q 10 -p -B -s 2`. The readcount is later converted to FPKM and FPKM-UQ using [GDC's formula].

Gene annotation in use is identical to GDC (`gencode.v22.annotation.gtf.gz`; md5: `291330bdcff1094bc4d5645de35e0871`), which is available on the [GDC Reference Files] page. Note that although the gene annotation is identical to GDC, GDC readcount is done by HTSeq under unstranded mode, so the readcounts reported by GDC and this pipeline are not identical.


## Output format
Each sample gets its TSV in the exact same gene order. The output TSV file has the following columns:

|   Column   |             Description              |
| :--------- | :----------------------------------- |
| gene_id    | Ensembl gene ID                      |
| read_count | stranded read count by featureCounts |
| fpkm       | FPKM                                 |
| fpkm_uq    | FPKM-UQ                              |

Use the gene information ([`gencode.gene.info.v22.tsv`][gene-info-tsv]; md5: `0a3f1d9b0a679e2a426de36d8d74fbf9`) on the [GDC Reference Files] page to get extra information about each gene.

Each batch of execution will produce a TSV `analysis_summary.dat` containing all output with the following columns:

|   Column    |      Description       |
| :---------- | :--------------------- |
| # case      | (same as the BAM map)  |
| disease     | (same as the BAM map)  |
| data_path   | Path to the output TSV |
| file_format | Always `TSV`           |
| sample_name | (same as the BAM map)  |
| sample_uuid | (same as the BAM map)  |

[GDC Reference Files]: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
[GDC's formula]: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#upper-quartile-fpkm
[gene-info-tsv]: https://api.gdc.cancer.gov/data/b011ee3e-14d8-4a97-aed4-e0b10f6bbe82


## Installation
Use the Docker image on Docker Hub [`lbwang/cptac_rna_expression`][docker-image]. See `Dockerfile` for the specific tool requirements.

[docker-image]: https://hub.docker.com/r/lbwang/cptac_rna_expression


## Input preparation
The pipeline requires the following inputs, which are passed through the `snakemake_config.json`:

```json
{
    "sample_list": "List of sample names",
    "bam_map": "Path to the (CPTAC) BAM map",
    "gdc_gtf": "Path to GDC's annotation GTF (e.g. gencode.v22.annotation.gtf)",
    "gdc_gene_info": "Path to GDC's gene info TSV (e.g. gencode.gene.info.v22.tsv)",
    "workflow_root": "Local file path to this repository (e.g. /path/to/cptac_rna_expression)"
}
```


## Pipeline execution
Clone this repository as `pipeline_workflow`:

    git clone https://github.com/ding-lab/cptac_rna_expression pipeline_workflow

The folder structure should be:

    <root>
    - pipeline_workflow
    - <batch1>
    - <batch2>


### Run the pipeline on katmai

Create the conda environment if not available (assuming [bioconda] is configured):

    conda create -n cptac_expression python=3.7 \
        snakemake-minimal=5.10.0 \
        pandas=1.0.1 \
        samtools=1.9 htslib=1.9 \
        subread=1.6.4

[bioconda]: https://bioconda.github.io/index.html

Create a new batch folder, say `cptac2_prospective`:

    mkdir cptac2_prospective

Create `snakemake_config.json`:

```json
{
    "sample_list": "samples.list",
    "bam_map": "bammap.converted.tsv",
    "gdc_gtf": "/diskmnt/Datasets/Reference/GDC/gencode.v22.annotation.gtf",
    "gdc_gene_info": "/diskmnt/Datasets/Reference/GDC/gencode.gene.info.v22.tsv",
    "workflow_root": "/diskmnt/Projects/cptac_scratch/CPTAC_expression/pipeline_workflow"
}
```

Copy the bam map (usually from Matt's CPTAC3.catalog):

    rsync -a --info=progress2 \
        katmai.BamMap.dat \
        katmai:<batch_location>/katmai.BamMap.cptac3_catalog_commit_$(git rev-parse --short HEAD).dat

Run the full pipeline by:

    snakemake \
        --configfile=snakemake_config.json \
        -s ../pipeline_workflow/Snakefile \
        -j 50 --resources io_heavy=4 -- \
        make_analysis_summary


### Run the pipeline on compute1
Refer to the example at `/storage1/fs1/lding/Active/CPTAC3/Analysis/rna_expression_pipeline/2020-02-25_PDA`.

Create a new batch by:

    cd /storage1/fs1/path/to/rna_expression_pipeline/

    # Create a batch name
    export BATCH="2020-02-25_PDA"

    # Create the batch folder
    mkdir $BATCH
    cd $BATCH

    # Create samples.list. For example,
    cut -f2 CPTAC3.catalog/BamMap/compute1.BamMap.dat | tail +2 > samples.list

    # Cache the BAM map, preferably add the current commit hash of CPTAC3.catalog in the file name.
    cp CPTAC3.catalog/BamMap/compute1.BamMap.dat compute1.BamMap.cptac3_catalog_commit_6509b69.dat

    # Create snakemake_config.json
    {
        "sample_list": "samples.list",
        "bam_map": "compute1.BamMap.cptac3_catalog_commit_6509b69.dat",
        "gdc_gtf": "../annotations/gencode.v22.annotation.gtf",
        "gdc_gene_info": "../annotations/gencode.gene.info.v22.tsv",
        "workflow_root": "/storage1/fs1/lding/Active/CPTAC3/Analysis/rna_expression_pipeline/pipeline_workflow"
    }

    # Set up the log structure
    mkdir logs
    mkdir logs/cluster

    # Copy the bash scripts run.sh and run_master_job.sh


Set up the snakemake profile. The default profile is at `/storage1/fs1/lding/Active/CPTAC3/Analysis/rna_expression_pipeline/ris_lsf`. Change the following files if necessary:

- `config.yaml` controls the number of restart times when the same job is failed and the total number of jobs submitted.
- `lsf-submit.py` sets extra LSF options. Note that job group is always required.

    ```python
    def generate_lsf_command(job_properties: dict, cluster: dict) -> str:
        # ...
        # RIS compute1 commands
        queue = "general"
        job_grp = "/liang-bo.wang/cptac_rna"
        docker_img = "lbwang/cptac_rna_expression"
        # ...
    ```

If a new snakemake profile is used, modify `run_master_job.sh` to point to the new profile.

Run the pipeline by:

    bash run.sh


## Apply the pipeline to other datasets
The BAMs should be aligned to hg38 (perferrably same as [GDC Reference Files]). The BAM map given by `bam_map` should mimic CPTAC3 catalog's BAM map columns, which must have the following columns in TSV format:

- `# sample_name`: Note `# ` is part of the column name;
                   sample name must be unique in the BAM map
- `case`
- `disease`
- `UUID`: can be any identifier if UUID is not available
- `data_pth`: the absolute path to the sample BAM;
              file name must be *.rna_seq.genomic.gdc_realn.bam
- `experimental_strategy`: must be `RNA-Seq`
- `data_format`: must be `BAM`
- `reference`: must be `hg38`