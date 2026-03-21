process STAR_INDEX {

  input:
  path reference_fasta
  path annotation_gtf

  output:
  path 'star_index', emit: index

  script:
  def gtfArg = annotation_gtf ? "--sjdbGTFfile \"${annotation_gtf}\"" : ''
  """
  mkdir -p star_index

  STAR \
    --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles "${reference_fasta}" \
    ${gtfArg} \
    --runThreadN ${params.star_index_threads} \
    --sjdbOverhang ${params.star_index_overhang} \
    ${params.star_index_extra}

  mkdir -p "${projectDir}/${params.outdir}/reference/star_index"
  cp -r star_index/. "${projectDir}/${params.outdir}/reference/star_index/"
  """

  stub:
  """
  mkdir -p star_index
  touch star_index/Genome
  touch star_index/SA
  touch star_index/SAindex

  mkdir -p "${projectDir}/${params.outdir}/reference/star_index"
  cp -r star_index/. "${projectDir}/${params.outdir}/reference/star_index/"
  """
}
