# Usage

## Run the pipeline

From the repository root:

```bash
nextflow run main.nf -params-file params.yml
```

This default run uses the local environment and does not require Singularity.
To run with the configured container instead:

```bash
nextflow run main.nf -params-file params.yml -profile singularity
```

Or override parameters on the CLI:

```bash
nextflow run main.nf \
  --samplesheet assets/samplesheet.de.example.csv \
  --counts_matrix assets/counts.example.tsv \
  --outdir results
```

To enable alignment, provide a STAR genome index directory:

```bash
nextflow run main.nf \
  -params-file params.yml \
  --star_index /abs/path/to/star_index
```

This scaffold currently includes:

- input validation
- samplesheet parsing
- FASTQ QC/trimming with `fastp`
- MultiQC aggregation of `fastp` reports
- optional STAR alignment of trimmed reads
- differential expression analysis with DESeq2

## QC parameters

Configured in `params.yml`:

- `fastp_threads`: thread count passed to `fastp`
- `fastp_extra`: optional additional `fastp` CLI arguments
- `star_index`: optional STAR genome index directory; if unset, alignment is skipped
- `star_threads`: thread count passed to STAR and `samtools index`
- `star_extra`: optional additional STAR CLI arguments

QC outputs are written to:

- `results/qc/fastp/`
- per-sample trimmed reads (`*.trimmed_R1.fastq.gz`, `*.trimmed_R2.fastq.gz`)
- per-sample reports (`*.fastp.html`, `*.fastp.json`)
- `results/qc/multiqc/`
- aggregated report (`multiqc_report.html`) and parsed data directory (`multiqc_data/`)

If `star_index` is set, alignment outputs are written to:

- `results/alignment/star/`
- coordinate-sorted BAMs (`*.sorted.bam`) and indexes (`*.sorted.bam.bai`)
- STAR summary logs (`*.Log.final.out`)
- splice junction tables (`*.SJ.out.tab`)

## Differential expression parameters

Configured in `params.yml`:

- `counts_matrix`: TSV with first column `gene` and remaining columns as sample IDs
- `condition_col`: samplesheet column used for group labels
- `control_level`: reference condition
- `treatment_level`: condition compared against control
- `fdr_threshold`: significance cutoff for filtered output

Expected samplesheet columns:

- `sample,fastq_1,fastq_2,strandedness,condition`

## Output files

- `results/differential_expression.tsv`: all genes with `padj` (Benjamini-Hochberg FDR)
- `results/differential_expression_fdr.tsv`: subset where `padj <= fdr_threshold`

## Run a fast validation execution

```bash
nextflow run main.nf -params-file params.yml -stub-run
```

## Build the Singularity/Apptainer container

Using Apptainer (recommended):

```bash
apptainer build containers/yeast-rnaseq.sif containers/Singularity.def
```

Using Singularity:

```bash
singularity build containers/yeast-rnaseq.sif containers/Singularity.def
```

After building, run the pipeline with:

- `-profile singularity`
- `container: file://./containers/yeast-rnaseq.sif` (via `params.yml`)

You can override container path at runtime if needed:

```bash
nextflow run main.nf -params-file params.yml --container file:///abs/path/to/container.sif
```
