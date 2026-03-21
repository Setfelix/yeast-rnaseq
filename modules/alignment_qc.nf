process ALIGNMENT_QC {
  tag "$sample"

  publishDir "${params.outdir}/qc/alignment", mode: 'copy'

  input:
  tuple val(sample), path(bam), path(bai), val(strandedness), val(condition)

  output:
  path "${sample}.flagstat.txt", emit: flagstat

  script:
  """
  samtools flagstat "${bam}" > "${sample}.flagstat.txt"
  """

  stub:
  """
  cat <<'EOF' > "${sample}.flagstat.txt"
  2 + 0 in total (QC-passed reads + QC-failed reads)
  2 + 0 primary
  2 + 0 mapped (100.00% : N/A)
  2 + 0 properly paired (100.00% : N/A)
  EOF
  """
}
