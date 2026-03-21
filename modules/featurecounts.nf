process FEATURECOUNTS {
  publishDir "${params.outdir}/counts", mode: 'copy'

  input:
  path bams
  path annotation_gtf

  output:
  path 'featurecounts.tsv', emit: counts
  path 'featurecounts.summary', emit: summary

  script:
  """
  featureCounts \
    -T ${params.featurecounts_threads} \
    -a "${annotation_gtf}" \
    -o featurecounts.raw.tsv \
    ${params.featurecounts_extra} \
    ${bams.join(' ')}

  awk 'BEGIN{FS=OFS="\t"} NR==1 {printf "gene"; for (i=7; i<=NF; i++) {n=split(\$i, parts, "/"); printf OFS parts[n]} printf "\\n"; next} NR>2 {printf \$1; for (i=7; i<=NF; i++) printf OFS \$i; printf "\\n"}' \
    featurecounts.raw.tsv > featurecounts.tsv

  cp featurecounts.raw.tsv.summary featurecounts.summary
  """

  stub:
  """
  cat <<'EOF' > featurecounts.tsv
  gene	ctrl_1.sorted.bam	ctrl_2.sorted.bam	treat_1.sorted.bam	treat_2.sorted.bam
  YAL001C	51	49	103	101
  YAL002W	60	58	22	20
  EOF
  cat <<'EOF' > featurecounts.summary
  Status	ctrl_1.sorted.bam	ctrl_2.sorted.bam	treat_1.sorted.bam	treat_2.sorted.bam
  Assigned	111	107	125	121
  EOF
  """
}
