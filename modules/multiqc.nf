process MULTIQC {
  publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

  input:
  path qc_reports

  output:
  path 'multiqc_report.html', emit: report
  path 'multiqc_data', emit: data

  script:
  """
  multiqc \
    --force \
    --outdir . \
    .
  """

  stub:
  """
  mkdir -p multiqc_data
  cat <<'EOF' > multiqc_report.html
  <html><body><h1>MultiQC stub report</h1></body></html>
  EOF
  cat <<'EOF' > multiqc_data/multiqc_fastp.txt
  sample	total_reads
  stub	1
  EOF
  """
}
