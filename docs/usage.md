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
  --outdir results
```

To use a prebuilt STAR genome index instead of generating one:

```bash
nextflow run main.nf \
  -params-file params.yml \
  --star_index /abs/path/to/star_index
```

Reference-related inputs can be declared in `params.yml`:

- `reference_fasta`: reference genome FASTA path
- `star_index`: STAR genome index directory
- `annotation_gtf`: gene annotation GTF path

If `star_index` is omitted, the pipeline builds it from `reference_fasta` and
`annotation_gtf`.

The default pipeline path is:

```bash
nextflow run main.nf \
  -params-file params.yml
```

To force DE to use an external count matrix instead:

```bash
nextflow run main.nf \
  -params-file params.yml \
  --de_counts_source external \
  --counts_matrix /abs/path/to/counts.tsv
```

This scaffold currently includes:

- input validation
- samplesheet parsing
- FASTQ QC/trimming with `fastp`
- MultiQC aggregation of `fastp` reports
- STAR alignment of trimmed reads
- post-alignment QC with `samtools flagstat`
- per-sample `featureCounts` quantification
- differential expression analysis with DESeq2

## QC parameters

Configured in `params.yml`:

- `fastp_threads`: thread count passed to `fastp`
- `fastp_extra`: optional additional `fastp` CLI arguments
- `reference_fasta`: reference genome FASTA path used when `star_index` is not supplied
- `star_index`: optional STAR genome index directory; if unset, the pipeline can build one from `reference_fasta` and `annotation_gtf`
- `star_index_threads`: thread count passed to STAR genome generation
- `star_index_overhang`: STAR `sjdbOverhang` used during genome generation
- `star_index_extra`: optional additional STAR genome generation CLI arguments
- `star_threads`: thread count passed to STAR and `samtools index`
- `star_extra`: optional additional STAR CLI arguments
- `annotation_gtf`: required annotation file for STAR index generation and `featureCounts`
- `featurecounts_threads`: thread count passed to `featureCounts`
- `featurecounts_extra`: optional additional `featureCounts` CLI arguments

QC outputs are written to:

- `results/qc/fastp/`
- per-sample trimmed reads (`*.trimmed_R1.fastq.gz`, `*.trimmed_R2.fastq.gz`)
- per-sample reports (`*.fastp.html`, `*.fastp.json`)
- `results/qc/multiqc/`
- aggregated report (`multiqc_report.html`) and parsed data directory (`multiqc_data/`)
- MultiQC includes `fastp` metrics and, when alignment runs, `samtools flagstat` summaries

Alignment outputs are written to:

- `results/alignment/star/`
- coordinate-sorted BAMs (`*.sorted.bam`) and indexes (`*.sorted.bam.bai`)
- STAR summary logs (`*.Log.final.out`)
- splice junction tables (`*.SJ.out.tab`)

If `star_index` is not set and `reference_fasta` is provided, the generated index
is written to:

- `results/reference/star_index/`

Post-alignment QC outputs are written to:

- `results/qc/alignment/`
- `samtools flagstat` summaries (`*.flagstat.txt`)

Count outputs are written to:

- `results/counts/per_sample/`
- per-sample count tables (`*.featurecounts.tsv`) and summaries (`*.featurecounts.summary`)

Those per-sample count tables are merged internally by the differential
expression step when `de_counts_source: generated` is used. The merged count
matrix is not currently published as a standalone file in `results/counts/`.

## Differential expression parameters

Configured in `params.yml`:

- `counts_matrix`: external TSV with first column `gene` and remaining columns as sample IDs; only used with `de_counts_source: external`
- `de_counts_source`: one of `external` or `generated`
- `condition_col`: samplesheet column used for group labels
- `control_level`: reference condition
- `treatment_level`: condition compared against control
- `fdr_threshold`: significance cutoff for filtered output

Default behavior is `de_counts_source: generated`.
Use `de_counts_source: external` only when you explicitly want to bypass generated counts.

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
