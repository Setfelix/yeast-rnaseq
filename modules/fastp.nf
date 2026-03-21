process FASTP {
  tag "$sample"

  publishDir "${params.outdir}/qc/fastp", mode: 'copy'

  input:
  tuple val(sample), val(fastq_1), val(fastq_2), val(strandedness), val(condition)

  output:
  tuple val(sample), path("${sample}.trimmed_R1.fastq.gz"), path("${sample}.trimmed_R2.fastq.gz"), val(strandedness), val(condition), emit: reads
  path "${sample}.fastp.html", emit: html
  path "${sample}.fastp.json", emit: json

  script:
  """
  fastp \
    --in1 "${fastq_1}" \
    --in2 "${fastq_2}" \
    --out1 "${sample}.trimmed_R1.fastq.gz" \
    --out2 "${sample}.trimmed_R2.fastq.gz" \
    --thread ${params.fastp_threads} \
    --html "${sample}.fastp.html" \
    --json "${sample}.fastp.json" \
    ${params.fastp_extra}
  """

  stub:
  """
  printf "@SEQ\\nACGT\\n+\\n!!!!\\n" | gzip -c > "${sample}.trimmed_R1.fastq.gz"
  printf "@SEQ\\nTGCA\\n+\\n!!!!\\n" | gzip -c > "${sample}.trimmed_R2.fastq.gz"
  cat <<'EOF' > "${sample}.fastp.html"
  <html><body><h1>fastp stub report</h1></body></html>
  EOF
  cat <<'EOF' > "${sample}.fastp.json"
  {"summary":{"before_filtering":{"total_reads":1},"after_filtering":{"total_reads":1}}}
  EOF
  """
}
