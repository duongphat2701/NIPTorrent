process create_ref_sanefalcon {

  input:
  val bamdir

  output:
  path("*"), emit: sanefalcon_ref

  publishDir "${params.outdir}/sanefalcon_ref", mode: 'copy'

  script:
  """
  export SANEFALCON_DIR="${projectDir}/assets/sanefalcon"
  export FF_REF_TSV="${params.ff_ref_tsv}"
  export OUTDIR="${params.outdir}/sanefalcon_ref"

  bash ${projectDir}/bin/sanefalcon.sh \
      --sanefalcon-dir \$SANEFALCON_DIR \
      --bamdir "$bamdir" \
      --outdir \$OUTDIR \
      --threads ${task.cpus} \
      --ff_ref_tsv \$FF_REF_TSV \
      --train
  """
}
