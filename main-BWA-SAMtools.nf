#!/usr/bin/env nextflow

params.reads = "./data/*_{1,2}.fastq.gz"
params.ref   = "./ref/ref.fa"
params.outdir = "./results"

process index_reference {
    input:
    path ref from file(params.ref)

    output:
    path("${ref}.bwt"), emit: index_files

    script:
    """
    bwa index $ref
    """
}

process bwa_mem_align {
    tag "${sample_id}"
    input:
    tuple val(sample_id), path(reads)
    path ref_index from index_reference.out.index_files.collect()

    output:
    path("${sample_id}.bam")

    script:
    """
    bwa mem -t 16 ${params.ref} ${reads[0]} ${reads[1]} | \
      samtools view -@ 4 -Sb - > ${sample_id}.bam
    """
}

process sort_bam {
    input:
    path bam

    output:
    path("${bam.simpleName}.sorted.bam")

    script:
    """
    samtools sort -@ 4 ${bam} -o ${bam.simpleName}.sorted.bam
    """
}

process index_bam {
    input:
    path sorted_bam

    output:
    path("${sorted_bam}.bai")

    script:
    """
    samtools index ${sorted_bam}
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, flat: true)
        .map { id, reads -> tuple(id.replaceAll(/_1|_2/, ""), reads) }
        .set { sample_reads }

    sample_reads
        | bwa_mem_align
        | sort_bam
        | index_bam
}
