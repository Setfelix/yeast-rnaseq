process FEATURECOUNTS {
  tag "$sample"

  publishDir "${params.outdir}/counts/per_sample", mode: 'copy'

  input:
  tuple val(sample), path(bam), path(bai), val(strandedness), val(condition)
  path annotation_gtf

  output:
  path "${sample}.featurecounts.tsv", emit: counts
  path "${sample}.featurecounts.summary", emit: summary

  script:
  """
  case "${strandedness}" in
    unstranded) strand_mode=0 ;;
    forward) strand_mode=1 ;;
    reverse) strand_mode=2 ;;
    *)
      echo "Unsupported strandedness '${strandedness}'. Use one of: unstranded, forward, reverse" >&2
      exit 1
      ;;
  esac

  featureCounts \
    -T ${params.featurecounts_threads} \
    -s "\$strand_mode" \
    -a "${annotation_gtf}" \
    -o "${sample}.featurecounts.raw.tsv" \
    ${params.featurecounts_extra} \
    "${bam}"

  awk 'BEGIN{FS=OFS="\t"} NR==1 {print "gene\t${sample}"; next} NR>2 {print \$1, \$7}' \
    "${sample}.featurecounts.raw.tsv" > "${sample}.featurecounts.tsv"

  awk 'BEGIN{FS=OFS="\t"} NR==1 {print "Status\t${sample}"; next} NR>1 {print \$1, \$2}' \
    "${sample}.featurecounts.raw.tsv.summary" > "${sample}.featurecounts.summary"
  """

  stub:
  """
  cat <<'EOF' > "${sample}.featurecounts.tsv"
  gene	${sample}
  YAL001C	51
  YAL002W	60
  EOF
  cat <<'EOF' > "${sample}.featurecounts.summary"
  Status	${sample}
  Assigned	111
  EOF
  """
}
