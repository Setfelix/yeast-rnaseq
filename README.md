# yeast-rnaseq-nf

Minimal, reproducible **Nextflow DSL2** RNA-seq pipeline from FASTQ to differential expression.

## Repository layout

- `main.nf`: DSL2 workflow wiring FASTQ QC and DE analysis
- `nextflow.config`: global config and includes for `conf/`
- `params.yml`: example runtime parameters
- `conf/`: base, Singularity, and profile-specific config
- `modules/fastp.nf`: FASTQ QC/trimming module using `fastp`
- `modules/multiqc.nf`: aggregate QC reporting with `MultiQC`
- `modules/star_index.nf`: STAR genome index generation when needed
- `modules/star_align.nf`: STAR alignment module
- `modules/alignment_qc.nf`: post-alignment QC with `samtools flagstat`
- `modules/infer_strandedness.nf`: infer strandedness with `infer_experiment.py` when samplesheet values are missing
- `modules/featurecounts.nf`: per-sample `featureCounts` quantification
- `containers/Singularity.def`: container recipe with core RNA-seq tools and DESeq2
- `assets/samplesheet.example.csv`: example input sheet
- `assets/samplesheet.de.example.csv`: DE-ready example input sheet with condition labels
- `assets/samplesheet.de.infer.example.csv`: DE-ready example input sheet without strandedness values
- `assets/counts.example.tsv`: example gene count matrix
- `.github/workflows/ci.yml`: CI checks (`nextflow config`, `-stub-run`)
- `docs/usage.md`: usage and container build instructions

## Quick start

```bash
nextflow run main.nf -params-file params.yml
```

For containerized execution with Singularity/Apptainer:

```bash
nextflow run main.nf -params-file params.yml -profile singularity
```

## Stub test run

```bash
nextflow run main.nf -params-file params.yml -stub-run
```

## Differential expression output

The pipeline writes:

- `results/differential_expression.tsv`: all genes with BH-FDR corrected `padj`
- `results/differential_expression_fdr.tsv`: genes passing `padj <= fdr_threshold`

By default, the pipeline runs:

- `FASTP`
- `STAR_INDEX` or a supplied `star_index`
- `STAR_ALIGN`
- `ALIGNMENT_QC`
- `FEATURECOUNTS`
- `MULTIQC`
- `DIFFERENTIAL_EXPRESSION`

Differential expression uses generated counts by default
(`--de_counts_source generated`). Use `--de_counts_source external` with
`--counts_matrix` only when you want to bypass the generated count matrix.

Reference-related paths are configured in `params.yml`:

- `reference_fasta`
- `star_index`
- `annotation_gtf`
- `annotation_bed12` (required only when strandedness must be inferred)

If `star_index` is not provided, the pipeline can build a STAR index from
`reference_fasta` and `annotation_gtf`.

## QC output

The pipeline writes per-sample `fastp` outputs to `results/qc/fastp/`:

- `${sample}.trimmed_R1.fastq.gz`
- `${sample}.trimmed_R2.fastq.gz`
- `${sample}.fastp.html`
- `${sample}.fastp.json`

It also writes an aggregated MultiQC report to `results/qc/multiqc/` covering
`fastp` outputs and, when alignment runs, `samtools flagstat` summaries:

- `multiqc_report.html`
- `multiqc_data/`

## Alignment output

The pipeline writes STAR alignment outputs to `results/alignment/star/`:

- `${sample}.sorted.bam`
- `${sample}.sorted.bam.bai`
- `${sample}.Log.final.out`
- `${sample}.SJ.out.tab`

If the STAR index is generated in-pipeline, it is written to
`results/reference/star_index/`.

It also writes post-alignment QC summaries to `results/qc/alignment/`:

- `${sample}.flagstat.txt`

When `strandedness` is omitted from the samplesheet and `annotation_bed12` is
provided, the pipeline infers strandedness after alignment and writes reports to
`results/qc/strandedness/`:

- `${sample}.infer_experiment.txt`

The pipeline writes per-sample count outputs to `results/counts/per_sample/`:

- `${sample}.featurecounts.tsv`
- `${sample}.featurecounts.summary`

These per-sample count tables are merged internally by the differential
expression step. The merged matrix is used as DESeq2 input, but it is not
currently published as a standalone file under `results/counts/`.

## License

This project is licensed under the MIT License. See `LICENSE`.
