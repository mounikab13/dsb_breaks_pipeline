#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.genome_ref = "~/test_project/genome_ref/"

process bwaIndex{
    publishDir params.genome_ref, mode: 'copy'

    input:
    tuple val(ref_name), path(reference) // Input reference genome file
  
    output:
    tuple val("${ref_name}"), path("${ref_name}.*")// Output BWA index files
  
    """
    bwa index ${reference} -p ${ref_name}
    """

}

workflow {

  Channel
      .fromPath(params.genome_ref + "*.fasta", checkIfExists: true)
      .map { tuple(it.getSimpleName(),it) }
      .set { input_ref_ch }

  Channel
      bwaindex(input_ref_ch)
      .set { indexed_bwa_ref_ch }

}
