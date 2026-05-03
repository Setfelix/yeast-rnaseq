nextflow.enable.dsl = 2

include { FASTP } from './modules/fastp'
include { MULTIQC } from './modules/multiqc'
include { STAR_INDEX } from './modules/star_index'
include { STAR_ALIGN } from './modules/star_align'
include { ALIGNMENT_QC } from './modules/alignment_qc'
include { INFER_STRANDEDNESS } from './modules/infer_strandedness'
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

  if (!['external', 'generated'].contains(params.de_counts_source)) {
    error "Invalid --de_counts_source: ${params.de_counts_source}. Use one of: external, generated"
  }

  if (params.de_counts_source == 'external') {
    if (!params.counts_matrix) {
      error "--de_counts_source external requires --counts_matrix"
    }
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

  if (!params.annotation_gtf) {
    error "Missing required parameter: --annotation_gtf"
  }

  def annotation = file(params.annotation_gtf)
  if (!annotation.exists()) {
    error "Annotation GTF not found: ${params.annotation_gtf}"
  }

  if (params.annotation_bed12) {
    def annotationBed12 = file(params.annotation_bed12)
    if (!annotationBed12.exists()) {
      error "Annotation BED12 not found: ${params.annotation_bed12}"
    }
  }

  if (!params.star_index && !params.reference_fasta) {
    error "Provide either --star_index or --reference_fasta"
  }

  if (params.de_counts_source == 'generated' && !params.annotation_gtf) {
    error "--de_counts_source generated requires --annotation_gtf"
  }
}

process DIFFERENTIAL_EXPRESSION {
  tag "${params.treatment_level}_vs_${params.control_level}"

  publishDir params.outdir, mode: 'copy'

  input:
  path count_inputs
  path samplesheet

  output:
  path 'differential_expression.tsv'
  path 'differential_expression_fdr.tsv'

  script:
  """
  Rscript -e '
  suppressPackageStartupMessages(library(DESeq2))

  meta <- read.csv("${samplesheet}", stringsAsFactors=FALSE)
  required <- c("sample", "${params.condition_col}")
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0) {
    stop(paste("Samplesheet missing required column(s):", paste(missing, collapse=", ")))
  }

  count_paths <- commandArgs(trailingOnly=TRUE)
  if (length(count_paths) == 0) {
    stop("No count inputs provided to differential expression")
  }

  if (length(count_paths) == 1) {
    counts <- read.table(count_paths[[1]], header=TRUE, sep="\t", check.names=FALSE)
    if (!("gene" %in% colnames(counts))) {
      stop("Counts matrix must include a gene column")
    }
  } else {
    count_tables <- lapply(count_paths, function(path) {
      tab <- read.table(path, header=TRUE, sep="\t", check.names=FALSE)
      if (ncol(tab) != 2 || !"gene" %in% colnames(tab)) {
        stop(paste("Per-sample counts file must have columns gene and sample:", path))
      }
      tab
    })
    counts <- Reduce(function(x, y) merge(x, y, by="gene", sort=FALSE), count_tables)
  }

  rownames(counts) <- counts[["gene"]]
  counts[["gene"]] <- NULL
  counts <- round(as.matrix(counts))

  if (!all(meta[["sample"]] %in% colnames(counts))) {
    missing_samples <- meta[["sample"]][!(meta[["sample"]] %in% colnames(counts))]
    stop(paste("Samplesheet sample(s) not in counts matrix:", paste(missing_samples, collapse=", ")))
  }

  meta <- meta[match(colnames(counts), meta[["sample"]]), , drop=FALSE]
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
  out[["gene"]] <- rownames(out)
  out <- out[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  out <- out[order(out[["padj"]], na.last=TRUE), ]

  write.table(out, file="differential_expression.tsv", sep="\t", quote=FALSE, row.names=FALSE)

  sig <- subset(out, !is.na(padj) & padj <= ${params.fdr_threshold})
  write.table(sig, file="differential_expression_fdr.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  ' ${count_inputs}
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
      def required = ['sample', 'fastq_1', 'fastq_2', params.condition_col as String]
      def missing = required.findAll { !row.containsKey(it) || !row[it] }
      if (missing) {
        error "Missing required column(s) in samplesheet row: ${missing.join(', ')}"
      }
      def fastq1 = file(row.fastq_1, checkIfExists: true)
      def fastq2 = file(row.fastq_2, checkIfExists: true)
      def strandedness = row.containsKey('strandedness') && row.strandedness ? row.strandedness.toString().trim().toLowerCase() : ''
      if (!strandedness && !params.annotation_bed12) {
        error "Samplesheet row for '${row.sample}' is missing strandedness. Provide a strandedness value or set --annotation_bed12 for inference."
      }
      tuple(row.sample, fastq1, fastq2, strandedness, row[params.condition_col])
    }

  def external_counts_ch = params.counts_matrix ? Channel.value(file(params.counts_matrix)) : null
  def generated_counts_ch = null
  def resolved_star_index_ch = null

  FASTP(samples_ch)

  def multiqc_inputs = FASTP.out.html.mix(FASTP.out.json)

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

  def counts_input_ch = STAR_ALIGN.out.bam
  if (params.annotation_bed12) {
    INFER_STRANDEDNESS(
      STAR_ALIGN.out.bam,
      file(params.annotation_bed12)
    )
    multiqc_inputs = multiqc_inputs.mix(INFER_STRANDEDNESS.out.report)
    counts_input_ch = INFER_STRANDEDNESS.out.resolved.map { sample, bam, bai, condition, strandedness_file ->
      tuple(sample, bam, bai, strandedness_file.text.trim(), condition)
    }
  }

  FEATURECOUNTS(
    counts_input_ch,
    file(params.annotation_gtf)
  )

  generated_counts_ch = FEATURECOUNTS.out.counts.collect()

  def counts_for_de
  if (params.de_counts_source == 'external') {
    counts_for_de = external_counts_ch
  } else {
    counts_for_de = generated_counts_ch
  }

  MULTIQC(
    multiqc_inputs.collect()
  )

  DIFFERENTIAL_EXPRESSION(
    counts_for_de,
    file(params.samplesheet)
  )
}
