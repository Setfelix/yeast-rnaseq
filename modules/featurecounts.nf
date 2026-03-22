process FEATURECOUNTS {
  publishDir "${params.outdir}/counts", mode: 'copy'

  input:
  tuple val(sample_header), val(strandedness), path(bams)
  path annotation_gtf

  output:
  path 'featurecounts.tsv', emit: counts
  path 'featurecounts.summary', emit: summary

  script:
  def strandModeMap = [
    'unstranded': '0',
    'forward'   : '1',
    'reverse'   : '2'
  ]

  def strandMode = strandModeMap[strandedness]
  if (!strandMode) {
    error "Unsupported strandedness '${strandedness}'. Use one of: unstranded, forward, reverse"
  }
  """
  featureCounts \
    -T ${params.featurecounts_threads} \
    -s ${strandMode} \
    -a "${annotation_gtf}" \
    -o featurecounts.raw.tsv \
    ${params.featurecounts_extra} \
    ${bams.join(' ')}

  awk 'BEGIN{FS=OFS="\t"} NR==1 {print "gene\t${sample_header}"; next} NR>2 {printf \$1; for (i=7; i<=NF; i++) printf OFS \$i; printf "\\n"}' \
    featurecounts.raw.tsv > featurecounts.tsv

  awk 'BEGIN{FS=OFS="\t"} NR==1 {print "Status\t${sample_header}"; next} NR>1 {print}' \
    featurecounts.raw.tsv.summary > featurecounts.summary

  """

  stub:
  """
  cat <<'EOF' > featurecounts.tsv
  gene	ctrl_1	ctrl_2	treat_1	treat_2
  YAL001C	51	49	103	101
  YAL002W	60	58	22	20
  EOF
  cat <<'EOF' > featurecounts.summary
  Status	ctrl_1	ctrl_2	treat_1	treat_2
  Assigned	111	107	125	121
  EOF
  """
}
