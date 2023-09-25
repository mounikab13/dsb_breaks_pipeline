#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.base_dir = "~/test_project/data_tmp/"
params.pub_dir = "~/test_project/data_tmp/dsb_breaks/"
params.raw_reads = "/test_project/raw_reads/coding-test-main-data-fastqs/data/fastqs/"
params.mapped_reads = params.base_dir + "/mapped_reads/"
params.genome_ref = "/test_project/genome_ref/"
params.bed_processed = params.base_dir + "/bed_processed/"
params.ref_bed = "/genome_ref/"
params.ref_prefix = "chr21"

process bwaMapping{
    memory '32 GB'
    cpus 30
    publishDir params.mapped_reads, mode:'copy'

    // Define input and output channels
    input:
    tuple val(sample_id), path(reads), path(genome_dir)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:

    """
    bwa mem -t ${task.cpus} -M -B 7 -w 33 -O 5 -E 2 -T 33 -Y ${genome_dir} ${fastq} |  samtools view -h -bS - | samtools sort > "${sample_id}.bam"
 
    """


}

process bamtobed{
    publishDir params.mapped_reads, mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bed")

    script:
    """

    bedtools bamtobed -i ${bam} > "${sample_id}.bed"
    """

}

process bedProcess {
    publishDir params.bed_processed, mode: 'copy'

    input:
    tuple val(sample_id), path(bed_file)

    output:
    path("${sample_id}_processed.bed")


    script:
    """
    #!/bin/bash
    awk 'BEGIN{OFS="\t"} \$6 == "+" { \$3 = \$2 + 1 } \$6 == "-" { \$2 = \$3 - 1 } { print }' ${bed_file} > "${sample_id}_processed.bed"
    """

}

process intersectBed{
    publishDir params.bed_processed, mode: 'copy'

    input:
    tuple val(sample_id), path(bed_file), path(ref_bed)

    output:
    path("${sample_id}_intersected.bed")

    script:
    """
    bedtools intersect -a ${ref_bed} -b ${bed_file} > "${sample_id}_intersected.bed"

    """

}


workflow{

    Channel
        .fromPath(params.genome_ref + "${params.ref_prefix}", checkIfExists: true) 
        .set { ref_ch }

    Channel
        .fromPath(params.raw_reads + "*.fastq.gz", checkIfExists: true) 
        .map { tuple(it.getSimpleName(),it) }
        .set { raw_reads_ch }

    Channel
        raw_reads_ch
        .combine(ref_ch)
        .set{ raw_reads_with_ref_ch }


    Channel
        bwaMapping(raw_reads_with_ref_ch)
        .set { mapped_reads_ch }

    Channel
        bamtobed(mapped_reads_ch)
        .set{ bed_ch }

    Channel
        bedProcess(bed_ch)
        .map { tuple(it.getSimpleName(),it) }
        .set{ bed_processed_ch }

    Channel
        .fromPath(params.ref_bed + "*.bed", checkIfExists: true)
        .set { ref_bed_ch }

    Channel
        bed_processed_ch
        .combine(ref_bed_ch)
        .set { bed_with_refbed_ch }


    Channel
        intersectBed(bed_with_refbed_ch)
        .map { tuple(it.getSimpleName(),it) }
        .set{ intersect_bed_ch }
        

}
