params.reads = "data/sample_R{1,2}.fastq.gz"
params.ref = "data/ref.fa"
params.known_sites = "data/dbsnp.vcf"

process BWA_MEM {
  input:
  tuple val(sample_id), path(reads)
  
  output:
  path("aligned.sorted.bam"), emit: sorted_bam
  
  script:
  """
  bwa mem -t 8 ${params.ref} ${reads[0]} ${reads[1]} | \
    samtools view -Sb - | \
    samtools sort -@ 4 -o aligned.sorted.bam
  samtools index aligned.sorted.bam
  """
}

process MarkDuplicates {
  input:
  path sorted_bam
  
  output:
  path("dedup.bam"), path("metrics.txt")
  
  script:
  """
  gatk MarkDuplicates \
    -I aligned.sorted.bam \
    -O dedup.bam \
    -M metrics.txt
  """
}

process BQSR {
  input:
  path dedup_bam
  
  output:
  path("recal_data.table")
  
  script:
  """
  gatk BaseRecalibrator \
    -I dedup.bam \
    -R ${params.ref} \
    --known-sites ${params.known_sites} \
    -O recal_data.table
  """
}

process ApplyBQSR {
  input:
  path dedup_bam
  path recal_table
  
  output:
  path("recal.bam")
  
  script:
  """
  gatk ApplyBQSR \
    -I dedup.bam \
    -R ${params.ref} \
    --bqsr-recal-file recal_data.table \
    -O recal.bam
  """
}

workflow {
  sample_id = "sample"
  reads = [ file("data/sample_R1.fastq.gz"), file("data/sample_R2.fastq.gz") ]

  BWA_MEM(sample_id, reads)
  MarkDuplicates(BWA_MEM.out.sorted_bam)
  BQSR(MarkDuplicates.out.dedup_bam)
  ApplyBQSR(MarkDuplicates.out.dedup_bam, BQSR.out.recal_data.table)
}
