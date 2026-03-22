nextflow.enable.dsl = 2

include { FASTP } from './modules/fastp'
include { MULTIQC } from './modules/multiqc'
include { STAR_INDEX } from './modules/star_index'
include { STAR_ALIGN } from './modules/star_align'
include { ALIGNMENT_QC } from './modules/alignment_qc'
include { FEATURECOUNTS } from './modules/featurecounts'

/*
 * Validate required runtime parameters.
 */
def validateParams() {
  if (!params.samplesheet) {
    error "Missing required parameter: --samplesheet"
  }
  def sheet = file(params.samplesheet)
  if (!sheet.exists()) {
    error "Samplesheet not found: ${params.samplesheet}"
  }

  if (!['auto', 'external', 'generated'].contains(params.de_counts_source)) {
    error "Invalid --de_counts_source: ${params.de_counts_source}. Use one of: auto, external, generated"
  }

  if (!params.counts_matrix && !params.annotation_gtf) {
    error "Provide either --counts_matrix or --annotation_gtf"
  }

  if (params.counts_matrix) {
    def counts = file(params.counts_matrix)
    if (!counts.exists()) {
      error "Counts matrix not found: ${params.counts_matrix}"
    }
  }

  if (!params.condition_col) {
    error "Missing required parameter: --condition_col"
  }
  if (!params.control_level || !params.treatment_level) {
    error "Missing required contrast levels: --control_level and --treatment_level"
  }

  if (params.star_index) {
    def starIndex = file(params.star_index)
    if (!starIndex.exists()) {
      error "STAR index directory not found: ${params.star_index}"
    }
  }

  if (params.reference_fasta) {
    def reference = file(params.reference_fasta)
    if (!reference.exists()) {
      error "Reference FASTA not found: ${params.reference_fasta}"
    }
    if (!params.star_index && !params.annotation_gtf) {
      error "--reference_fasta requires --annotation_gtf when STAR index generation is needed"
    }
  }

  if (params.annotation_gtf) {
    def annotation = file(params.annotation_gtf)
    if (!annotation.exists()) {
      error "Annotation GTF not found: ${params.annotation_gtf}"
    }
    if (!params.star_index && !params.reference_fasta) {
      error "--annotation_gtf requires either --star_index or --reference_fasta"
    }
  }

  if (params.de_counts_source == 'external' && !params.counts_matrix) {
    error "--de_counts_source external requires --counts_matrix"
  }

  if (params.de_counts_source == 'generated' && !params.annotation_gtf) {
    error "--de_counts_source generated requires --annotation_gtf and either --star_index or --reference_fasta"
  }
}

process DIFFERENTIAL_EXPRESSION {
  tag "${params.treatment_level}_vs_${params.control_level}"

  publishDir params.outdir, mode: 'copy'

  input:
  path counts_matrix
  path samplesheet

  output:
  path 'differential_expression.tsv'
  path 'differential_expression_fdr.tsv'

  script:
  """
  Rscript -e '
  suppressPackageStartupMessages(library(DESeq2))

  counts <- read.table("${counts_matrix}", header=TRUE, sep="\t", check.names=FALSE)
  if (!("gene" %in% colnames(counts))) {
    stop("Counts matrix must include a gene column")
  }

  rownames(counts) <- counts$gene
  counts$gene <- NULL
  counts <- round(as.matrix(counts))

  meta <- read.csv("${samplesheet}", stringsAsFactors=FALSE)
  required <- c("sample", "${params.condition_col}")
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0) {
    stop(paste("Samplesheet missing required column(s):", paste(missing, collapse=", ")))
  }

  if (!all(meta$sample %in% colnames(counts))) {
    missing_samples <- meta$sample[!(meta$sample %in% colnames(counts))]
    stop(paste("Samplesheet sample(s) not in counts matrix:", paste(missing_samples, collapse=", ")))
  }

  meta <- meta[match(colnames(counts), meta$sample), , drop=FALSE]
  meta[["${params.condition_col}"]] <- factor(meta[["${params.condition_col}"]], levels=c("${params.control_level}", "${params.treatment_level}"))

  if (sum(meta[["${params.condition_col}"]] == "${params.control_level}", na.rm=TRUE) < 2 ||
      sum(meta[["${params.condition_col}"]] == "${params.treatment_level}", na.rm=TRUE) < 2) {
    stop("Need at least two samples per group for DESeq2 in this scaffold")
  }

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = as.formula(paste0("~", "${params.condition_col}"))
  )

  dds <- DESeq(dds)
  res <- results(dds, contrast=c("${params.condition_col}", "${params.treatment_level}", "${params.control_level}"))

  out <- as.data.frame(res)
  out$gene <- rownames(out)
  out <- out[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  out <- out[order(out$padj, na.last=TRUE), ]

  write.table(out, file="differential_expression.tsv", sep="\t", quote=FALSE, row.names=FALSE)

  sig <- subset(out, !is.na(padj) & padj <= ${params.fdr_threshold})
  write.table(sig, file="differential_expression_fdr.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  '
  """

  stub:
  """
  cat <<'EOF' > differential_expression.tsv
  gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
  YAL001C	100	1.20	0.25	4.80	1.6e-06	3.2e-06
  YAL002W	85	-0.90	0.30	-3.00	0.0027	0.0041
  EOF

  cat <<'EOF' > differential_expression_fdr.tsv
  gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
  YAL001C	100	1.20	0.25	4.80	1.6e-06	3.2e-06
  YAL002W	85	-0.90	0.30	-3.00	0.0027	0.0041
  EOF
  """
}

workflow {
  validateParams()

  def samples_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
      def required = ['sample', 'fastq_1', 'fastq_2', 'strandedness', params.condition_col as String]
      def missing = required.findAll { !row.containsKey(it) || !row[it] }
      if (missing) {
        error "Missing required column(s) in samplesheet row: ${missing.join(', ')}"
      }
      tuple(row.sample, row.fastq_1, row.fastq_2, row.strandedness, row[params.condition_col])
    }

  def external_counts_ch = params.counts_matrix ? Channel.value(file(params.counts_matrix)) : null
  def generated_counts_ch = null
  def resolved_star_index_ch = null

  FASTP(samples_ch)

  def multiqc_inputs = FASTP.out.html.mix(FASTP.out.json)

  if (params.star_index || params.reference_fasta) {
    if (params.star_index) {
      resolved_star_index_ch = Channel.value(file(params.star_index))
    } else {
      STAR_INDEX(
        file(params.reference_fasta),
        file(params.annotation_gtf)
      )
      resolved_star_index_ch = STAR_INDEX.out.index
    }

    def star_reads_ch = FASTP.out.reads.combine(resolved_star_index_ch).map { sample, read_1, read_2, strandedness, condition, star_index ->
      tuple(sample, read_1, read_2, strandedness, condition, star_index)
    }

    STAR_ALIGN(star_reads_ch)

    ALIGNMENT_QC(STAR_ALIGN.out.bam)

    multiqc_inputs = multiqc_inputs.mix(ALIGNMENT_QC.out.flagstat)

    if (params.annotation_gtf) {
      def featurecounts_input_ch = STAR_ALIGN.out.bam
        .map { sample, bam, bai, strandedness, condition ->
          [
            sample      : sample as String,
            bam         : bam,
            strandedness: strandedness as String
          ]
        }
        .collect()
        .map { bamRows ->
          def sampleNames = bamRows.collect { row -> row.sample }
          def strandednessValues = bamRows.collect { row -> row.strandedness }.toSet().toList()
          if (strandednessValues.size() != 1) {
            error "FEATURECOUNTS currently requires the same strandedness for all samples; found: ${strandednessValues.join(', ')}"
          }
          def bamPaths = bamRows.collect { row -> row.bam }
        tuple(sampleNames.join('\t'), strandednessValues[0], bamPaths)
        }

      FEATURECOUNTS(
        featurecounts_input_ch,
        file(params.annotation_gtf)
      )
      generated_counts_ch = FEATURECOUNTS.out.counts
    }
  }

  def counts_for_de
  switch (params.de_counts_source) {
    case 'external':
      counts_for_de = external_counts_ch
      break
    case 'generated':
      counts_for_de = generated_counts_ch
      break
    default:
      counts_for_de = generated_counts_ch ?: external_counts_ch
  }

  MULTIQC(
    multiqc_inputs.collect()
  )

  DIFFERENTIAL_EXPRESSION(
    counts_for_de,
    file(params.samplesheet)
  )
}
