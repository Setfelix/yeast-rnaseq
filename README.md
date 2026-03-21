# yeast-rnaseq-nf

Minimal, reproducible scaffold for a **Nextflow DSL2** yeast RNA-seq pipeline.

## Repository layout

- `main.nf`: DSL2 workflow wiring FASTQ QC and DE analysis
- `nextflow.config`: global config and includes for `conf/`
- `params.yml`: example runtime parameters
- `conf/`: base, Singularity, and profile-specific config
- `modules/fastp.nf`: FASTQ QC/trimming module using `fastp`
- `modules/multiqc.nf`: aggregate QC reporting with `MultiQC`
- `modules/star_align.nf`: optional STAR alignment module
- `containers/Singularity.def`: container recipe with core RNA-seq tools and DESeq2
- `assets/samplesheet.example.csv`: example input sheet
- `assets/samplesheet.de.example.csv`: DE-ready example input sheet with condition labels
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

## QC output

The pipeline writes per-sample `fastp` outputs to `results/qc/fastp/`:

- `${sample}.trimmed_R1.fastq.gz`
- `${sample}.trimmed_R2.fastq.gz`
- `${sample}.fastp.html`
- `${sample}.fastp.json`

It also writes an aggregated MultiQC report to `results/qc/multiqc/`:

- `multiqc_report.html`
- `multiqc_data/`

## Optional alignment output

If `--star_index` is provided, the pipeline also writes STAR alignment outputs to
`results/alignment/star/`:

- `${sample}.sorted.bam`
- `${sample}.sorted.bam.bai`
- `${sample}.Log.final.out`
- `${sample}.SJ.out.tab`

## License

This project is licensed under the MIT License. See `LICENSE`.
