process STAR_ALIGN {
  tag "$sample"

  publishDir "${params.outdir}/alignment/star", mode: 'copy'

  input:
  tuple val(sample), path(read_1), path(read_2), val(strandedness), val(condition), path(star_index)

  output:
  tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), val(strandedness), val(condition), emit: bam
  path "${sample}.Log.final.out", emit: log
  path "${sample}.SJ.out.tab", emit: junctions

  script:
  """
  STAR \
    --genomeDir "${star_index}" \
    --readFilesIn "${read_1}" "${read_2}" \
    --readFilesCommand zcat \
    --runThreadN ${params.star_threads} \
    --outFileNamePrefix "${sample}." \
    --outSAMtype BAM SortedByCoordinate \
    ${params.star_extra}

  mv "${sample}.Aligned.sortedByCoord.out.bam" "${sample}.sorted.bam"
  samtools index -@ ${params.star_threads} "${sample}.sorted.bam"
  """

  stub:
  """
  printf "BAM_STUB\\n" > "${sample}.sorted.bam"
  printf "BAI_STUB\\n" > "${sample}.sorted.bam.bai"
  cat <<'EOF' > "${sample}.Log.final.out"
  Started job on | stub
  Number of input reads | 2
  Uniquely mapped reads number | 2
  EOF
  cat <<'EOF' > "${sample}.SJ.out.tab"
  chrI	10	20	1	1	0	1	0	10
  EOF
  """
}
